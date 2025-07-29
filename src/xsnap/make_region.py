#!/usr/bin/env python
"""
Pipeline module to generate DS9-compliant ICRS source and background regions.

This script produces two region files for automated processing:

  1. source.reg
     • Circular region at the specified RA/Dec with radius r_in (″).
  2. bkg.reg
     • First-pass annulus from r_in to r_out (″).
     • If an XIMAGE detection (excluding the source) overlaps, relocates the
       background to a single circle of radius r_out (″) at the first clean spot
       within around 1′ from the source.
Command-line usage:
    make_region <evtfile> <ra> <dec> [r_in] [r_out] [outdir]
                [--ds9 DS9_PATH] [--expimg EXP_IMG]

Positional arguments:
  evtfile    Path to event FITS file.
  ra, dec    Source coordinates in ICRS (decimal degrees).
  r_in       Inner radius of source region (arcsec; defaults telescope-specific).
  r_out      Outer radius for background region (arcsec; defaults telescope-specific).
  outdir     Output directory (default: same as evtfile’s directory).

Options:
  --ds9      Path to DS9 executable (if not on $PATH or in $DS9/_PATH).
  --expimg   Exposure-map FITS file for contamination checks.

Automates region creation and contamination checking for headless, repeatable workflows.
"""

import argparse, subprocess, sys, textwrap
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import os, shutil, subprocess, re, shlex, pathlib


# ------------------------------- CLI -------------------------------------------
def get_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__)
    )
    p.add_argument("evtfile")
    p.add_argument("ra")
    p.add_argument("dec")
    p.add_argument("r_in",  nargs="?", type=float)
    p.add_argument("r_out", nargs="?", type=float)
    p.add_argument("outdir", nargs="?")
    p.add_argument("--ds9")
    p.add_argument("--expimg") 
    return p.parse_args()

# ----------------------- telescope defaults & DS9 lookup ---------------------------------
def default_radii(tel: str) -> tuple[float, float]:
    if tel.lower().startswith("chandra"):
        return (2.0, 45.0)
    elif tel.lower().startswith("nustar"):
        return (45.0, 125.0)
    return (25.0, 125.0)

def _shell_resolve(shell: str) -> str | None:
    """Return ds9 path seen by *shell* or None."""
    if not pathlib.Path(shell).exists():
        return None
    try:
        # command -v prints the target of aliases, funcs, hashes, executables
        cmd = f'{shell} -lc "command -v ds9"'
        path = subprocess.check_output(cmd, shell=True, text=True).strip()
        if not path:
            return None
        # If it's an alias string like:  alias ds9='/Applications/.../ds9'
        m = re.match(r"alias ds9=['\"]?([^'\"]+)['\"]?", path)
        return m.group(1) if m else path       # plain path or extracted target
    except subprocess.CalledProcessError:
        return None

def find_ds9(cli_path: str | None = None) -> str | None:
    """Locate a usable DS9 binary, even if defined only by shell alias."""
    # 1. explicit flag
    if cli_path:
        return cli_path

    # 2. environment variables
    for env in ("DS9", "DS9_PATH"):
        if env in os.environ and os.environ[env]:
            return os.environ[env]

    # 3. normal PATH lookup
    exe = shutil.which("ds9")
    if exe:
        return exe

    # 4. ask common shells to resolve aliases / functions / hashes
    for sh in ("/bin/bash", "/bin/zsh"):
        resolved = _shell_resolve(sh)
        if resolved:
            return resolved

    return None  # no DS9 found, return None

def guess_exposure(evt: pathlib.Path) -> pathlib.Path:
    """
    Return the best-guess exposure image for *evt* using this priority:

    1.   <event>.expmap
    2.   <event>_ex.img
    3.   first *.img file in the same directory (alphabetical)
    4.   fall back to the event file itself
    """
    dir_ = evt.parent
    cand: list[pathlib.Path] = []

    # 1) <stem>.expmap
    p1 = evt.with_suffix(".expmap")
    if p1.is_file():
        cand.append(p1)

    # 2) <stem>_ex.img
    p2 = evt.with_name(evt.stem + "_ex.img")
    if p2.is_file():
        cand.append(p2)

    # 3) any *.img – pick the first alphabetically
    for img in sorted(dir_.glob("*.img")):
        if img not in cand:
            cand.append(img)
            break                       # only need the first one

    return cand[0] if cand else evt

# ------------------------------- region writer -------------------------------------------
HEADER = (
    "# Region file format: DS9 version 4.1\n"
    "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" "
    "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"
    "icrs\n"
)

def write_icrs_regions(outdir: pathlib.Path,
                       src: SkyCoord,
                       r_in: float, r_out: float,
                       bkg_center: SkyCoord | None = None,
                       bkg_annulus: bool = False   # ←★ new flag
                      ) -> tuple[pathlib.Path, pathlib.Path]:
    """Write source.reg & bkg.reg (annulus or circle)."""
    outdir.mkdir(parents=True, exist_ok=True)
    src_f = outdir / "source.reg"
    bkg_f = outdir / "bkg.reg"

    src_f.write_text(HEADER + f"circle({src.ra.deg:.8f},{src.dec.deg:.8f},{r_in}\")\n")

    if bkg_center is None:
        bkg_center = src

    if bkg_annulus:                                      # ←★ annulus first pass
        bkg_f.write_text(
            HEADER + f"annulus({src.ra.deg:.8f},{src.dec.deg:.8f},{r_in}\",{r_out}\")\n"
        )
    else:
        bkg_f.write_text(
            HEADER + f"circle({bkg_center.ra.deg:.8f},{bkg_center.dec.deg:.8f},{r_out}\")\n"
        )
    return src_f, bkg_f

# ------------------------------ XIMAGE detect wrapper --------------------------------------------
def _detect(evt: pathlib.Path, exp: pathlib.Path,
            src_center: SkyCoord,             # ←★ need source for exclusion
            bkg_center: SkyCoord, r_out: float, r_det: float,
            src_exclude: float = 30.0) -> bool:
    """
    Return True if *any* detection, **other than the one closest to the source
    within src_exclude″**, overlaps the background circle.
    """
    if "HEADAS" not in os.environ:
        print("⚠️  HEADAS not initialised – skipping contamination check.")
        return False

    script = f"""
    source $HEADAS/headas-init.sh
    ximage << EOF
      read/ecol=pi/emin=30/emax=1000/size=600 {evt}
      read/exposure/size=600 {exp}
      detect
    exit
    EOF
    """
    out = subprocess.run(["bash", "-lc", script], capture_output=True, text=True).stdout

    det_ra, det_dec = [], []
    for ln in out.splitlines():
        if not re.match(r"\s*\d+\s", ln):
            continue
        p = ln.split()
        if len(p) < 11:
            continue
        rah, ram, ras = map(float, p[5:8])
        dd, dm, ds    = map(float, p[8:11])
        det_ra.append((rah + ram/60 + ras/3600) * 15.0)
        sign = 1 if dd >= 0 else -1
        det_dec.append(sign * (abs(dd) + dm/60 + ds/3600))

    if not det_ra:
        return False

    det = SkyCoord(det_ra, det_dec, unit="deg")
    src_sep = src_center.separation(det).to(u.arcsec).value

    # exclude the closest detection if it's within src_exclude″
    if (src_sep <= src_exclude).any():
        idx_excl = np.argmin(src_sep)
        det = det[np.arange(len(det)) != idx_excl]   # drop that one
        if len(det) == 0:
            return False

    sep = bkg_center.separation(det).to_value(u.arcsec)
    return bool((sep <= (r_out + r_det)).any())

def find_clean_bkg_center(src: SkyCoord,
                          evt: pathlib.Path, exp: pathlib.Path,
                          r_in: float, r_out: float, gap: float,
                          telescope: str,
                          ang_step: int = 30,
                          rad_step: float = 10.0) -> SkyCoord:
    """Try directions every `ang_step`; grow radius by `rad_step` until clean."""
    r_det = 5.0 if "chandra" in telescope.lower() else 30.0
    base  = r_in + r_out + gap
    cosd  = np.cos(src.dec.radian)

    for extra in np.arange(0, 361, rad_step):
        sep_deg = (base + extra) / 3600
        for th in range(0, 360, ang_step):
            th_rad = np.deg2rad(th)
            dra  =  sep_deg * np.cos(th_rad) / cosd
            ddec =  sep_deg * np.sin(th_rad)
            cand = SkyCoord(src.ra + dra * u.deg, src.dec + ddec * u.deg, frame="icrs")
            if not _detect(evt, exp, src, cand, r_out, r_det):
                return cand
    return cand                       # fallback last tried

# --------------------------------------------------------------------------
def main() -> None:
    ns  = get_args()
    evt = pathlib.Path(ns.evtfile).expanduser()
    if not evt.is_file():
        sys.exit(f"{evt} not found")

    if ns.expimg:
        exp = pathlib.Path(ns.expimg).expanduser()
        if not exp.is_file():
            sys.exit(f"❌  exposure image {exp} not found")
    else:
        exp = guess_exposure(evt)

    with fits.open(evt) as hd:
        telescope = hd[0].header.get("TELESCOP", "unknown")

    r_in, r_out = (ns.r_in, ns.r_out) if ns.r_in and ns.r_out else default_radii(telescope)
    if r_in <= 0 or r_out <= r_in:
        sys.exit("Need 0 < r_in < r_out")

    src  = SkyCoord(ns.ra, ns.dec, unit=("deg", "deg"), frame="icrs")
    outd = pathlib.Path(ns.outdir) if ns.outdir else evt.parent

    # first pass – annulus
    src_reg, bkg_reg = write_icrs_regions(outd, src, r_in, r_out, bkg_annulus=True)

    r_det = 5.0 if "chandra" in telescope.lower() else 30.0
    contaminated = _detect(evt, exp, src, src, r_out, r_det)

    if contaminated:
        print("❗ Contaminating detection (excluding the source) – relocating background.")
        new_cen = find_clean_bkg_center(src, evt, exp, r_in, r_out,
                                        gap=60.0, telescope=telescope)
        src_reg, bkg_reg = write_icrs_regions(outd, src, r_in, r_out,
                                              bkg_center=new_cen, bkg_annulus=False)

    # launch DS9
    ds9 = find_ds9(ns.ds9)
    print(ds9)
    if ds9:
        src_p = outd / "source_physical.reg"
        bkg_p = outd / "bkg_physical.reg"
        cmd = [ds9, str(evt),
               "-region", "load", str(src_reg),
               "-region", "load", str(bkg_reg),
               "-scale", "log",
               "-region", "system", "physical",
               "-region", "save", str(src_p),
               "-region", "save", str(bkg_p),]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd)
    else:
        print("DS9 not found; regions written but viewer not launched.")
        
def cli() -> int | None:        
    return main()      

if __name__ == "__main__":
    main()