#!/usr/bin/env python3
"""
Pipeline module to run NuSTAR reduction and spectral extraction in one step.

Automates the standard HEASOFT workflow for NuSTAR:
  1. Execute nupipeline for a given OBSID, producing cleaned FPMA/B events.
  2. Run nuproducts to generate source and background spectra plus ARFs/RMFs.
  3. Optionally reuse pre-cleaned events (skip pipeline).

Usage
-----
    extract-nustar <OBSID> <SRC_REG> <BKG_REG> \
        [--indir DIR] [--outdir DIR] \
        [--ra RA --dec DEC] [--no-pipe]

Positional arguments
--------------------
  OBSID        11-digit NuSTAR observation ID.
  SRC_REG      Source region file or literal DS9 region string.
  BKG_REG      Background region file or literal DS9 region string.

Important options
-----------------
  --indir DIR        Directory containing raw event files (default: ./sources/<OBSID>).
  --outdir DIR       Output root for products (default: ./products/<OBSID>).
  --ra / --dec       Source coordinates (deg) passed to nupipeline.
  --no-pipe          Skip nupipeline; assume cleaned events already exist.

Outputs
-------
`FPMA/` and `FPMB/` sub-directories inside outdir with:
    nu<OBSID>A01_sr.pha   nu<OBSID>B01_sr.pha
which both of them are the spectrum file that will be used in :py:class:`~xsnap.spectrum.SpectrumFit`

Requirements
-------------
  • HEASoft with NuSTAR CALDB files installed

"""
from __future__ import annotations
import argparse, subprocess, shutil, pathlib, re, sys, os, textwrap
from typing import Tuple, Optional

TOOLS = ("nupipeline", "nuproducts")
CALDB = os.environ.get("CALDB")

# ───────────────────────── helpers ──────────────────────────────────────────
def parse_region_center(rfile: pathlib.Path) -> Tuple[float, float]:
    """Return (RA,Dec) in deg from the FIRST circle/annulus in a DS9 reg file."""
    def sex2deg(token: str) -> float:
        if ":" not in token:
            return float(token)
        parts = [float(x) for x in token.split(":")]
        if len(parts) == 3:
            sign = -1 if parts[0] < 0 else 1
            parts[0] = abs(parts[0])
            val = parts[0] + parts[1]/60 + parts[2]/3600
            # assume RA if <=24h
            if token.count(":") == 2 and val <= 24:
                val *= 15
            return sign * val
        raise ValueError("Unrecognized sexagesimal token: %s" % token)

    pat = re.compile(r"(circle|annulus)\s*\(\s*([^)]*)\)")
    for ln in rfile.read_text().splitlines():
        m = pat.search(ln)
        if m:
            inside = m.group(2)
            tokens = [t.strip() for t in inside.split(",")]
            if len(tokens) >= 2:
                return sex2deg(tokens[0]), sex2deg(tokens[1])
    raise RuntimeError(f"No circle/annulus line with RA/Dec in {rfile}")


def sh(cmd: list[str], **kw):
    print("➜", " ".join(cmd))
    subprocess.run(cmd, check=True, **kw)


def find_default_region(obsid: str, indir: pathlib.Path, suffix: str) -> Optional[pathlib.Path]:
    """
    Search for region file matching nu<obsid><suffix>_src.reg under indir/event_cl.
    suffix: 'A01' or 'B01'
    """
    pattern = f"nu{obsid}{suffix}_src.reg"
    for root, dirs, files in os.walk(indir):
        for fname in files:
            if fname.lower() == pattern.lower():
                return pathlib.Path(root) / fname
    return None


def runpipeline(obsid, ra, dec, indir, outdir):
    stem = f"nu{obsid}"
    cmd = [
        "nupipeline",
        f"indir={indir}", f"outdir={outdir}",
        f"steminputs={stem}"
    ]
    if ra is not None and dec is not None:
        cmd += [f"srcra={ra}", f"srcdec={dec}"]
    sh(cmd)


def runproducts(obsid, inst, outdir, prodroot, region):
    stem = f"nu{obsid}"
    prodroot.mkdir(parents=True, exist_ok=True)
    sh([
        "nuproducts",
        f"indir={outdir}", f"outdir={prodroot}",
        f"steminputs={stem}", f"instrument={inst}",
        f"srcregionfile={region}",
        f"bkgregionfile={outdir/'bkg.reg'}",
        "bkgextract=yes", "clobber=yes"
    ])

# ───────────────────────── core routine ─────────────────────────────────────
def extract_nustar(obsid: str,
                   src_reg: Optional[pathlib.Path],
                   bkg_reg: Optional[pathlib.Path],
                   *, indir: pathlib.Path,
                         outdir: pathlib.Path,
                         prodroot: pathlib.Path,
                         ra: Optional[float],
                         dec: Optional[float],
                         run_pipe: bool):
    # pre-check
    if any(shutil.which(t) is None for t in TOOLS):
        sys.exit("❌  NuSTAR FTOOLS not found.")
    if CALDB is None:
        sys.exit("❌  CALDB not set")

    indir = indir.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    prodroot = prodroot.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Determine source region for FPMA/B
    if src_reg is None:
        regA = find_default_region(obsid, indir, 'A01')
        regB = find_default_region(obsid, indir, 'B01')
    else:
        regA = regB = src_reg
    if not regA or not regB:
        sys.exit(f"❌  Could not find default region files for {obsid} in {indir}")

    # copy background region if provided
    if bkg_reg:
        (outdir / 'bkg.reg').write_text(bkg_reg.read_text())
    else:
        # optional: assume a default background
        (outdir / 'bkg.reg').write_text('')

    # copy source regions
    (outdir / 'source_A.reg').write_text(regA.read_text())
    (outdir / 'source_B.reg').write_text(regB.read_text())

    if run_pipe:
        if src_reg is None:
            print("No source.reg supplied or found: running nupipeline with default pointing")
            sh([
                "nupipeline",
                f"indir={indir}", f"outdir={outdir}",
                f"steminputs=nu{obsid}"
            ])
        else:
            # only parse/override RA/Dec if we have a region
            if ra is None or dec is None:
                ra, dec = parse_region_center(regA)
                print(f"(RA,Dec) taken from region: {ra:.6f}, {dec:.6f}")
            else:
                print(f"(RA,Dec) supplied: {ra:.6f}, {dec:.6f}")

            runpipeline(obsid, ra, dec, indir, outdir)

    # extract products for FPMA and FPMB
    runproducts(obsid, 'FPMA', outdir, prodroot / 'FPMA', outdir / 'source_A.reg')
    runproducts(obsid, 'FPMB', outdir, prodroot / 'FPMB', outdir / 'source_B.reg')

    print(f"\n✓ DONE – products under {prodroot}")

# ───────────────────────── CLI entry-point ──────────────────────────────────

def cli() -> None:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__)
    )
    parser.add_argument('obsid', help='NuSTAR observation ID')
    parser.add_argument('--src', type=pathlib.Path,
                    help='Source region file (optional)')
    parser.add_argument('--bkg', type=pathlib.Path, required= True,
                    help='Background region file')
    parser.add_argument('--indir', default=None)
    parser.add_argument('--outdir', default=None)
    parser.add_argument('--prod', default=None)
    parser.add_argument('--ra', type=float)
    parser.add_argument('--dec', type=float)
    parser.add_argument('--no-pipe', action='store_true')
    args = parser.parse_args()

    indir = pathlib.Path(args.indir or f"./sources/{args.obsid}")
    outdir = pathlib.Path(args.outdir or f"./output/{args.obsid}")
    prodroot = pathlib.Path(args.prod or f"./products/{args.obsid}")

    extract_nustar(
        args.obsid,
        args.src,
        args.bkg,
        indir=indir,
        outdir=outdir,
        prodroot=prodroot,
        ra=args.ra,
        dec=args.dec,
        run_pipe=not args.no_pipe
    )

if __name__ == '__main__':
    cli()
