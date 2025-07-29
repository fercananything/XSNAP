"""
Pipeline module to extract 0.3–10 keV Swift/XRT spectra in PC or WT mode.

Automates the standard HEASoft workflow:
  1. (Optionally) run xrtpipeline to calibrate a raw OBSID.
  2. Filter cleaned events with user-supplied source / background regions.
  3. Create exposure-corrected ARFs and grouped PHA files.

Quick usage
-----------
Full calibration + extraction (PC mode):

    extract-swift 00012345001 src.reg bkg.reg \
        --src-ra 150.123 --src-dec -12.345 --mode PC

Skip the pipeline, work on existing events/images:

    extract-swift --no-pipe 00012345001 src.reg bkg.reg \
        --indir ./events --evt sw00012345001xpcw3po_cl.evt

Positional arguments
--------------------
  OBSID        11-digit Swift observation ID.
  SRC_REG      Source region (DS9 “physical” or “sky”).
  BKG_REG      Background region.

Key options
-----------
  --indir DIR        Directory with event/image files (default: ./<OBSID>).
  --outdir DIR       Destination for PHA/ARF products (default: ./).
  --src-ra / --src-dec
                     Only needed when **xrtpipeline** is run.
  --mode {PC,WT}     Extraction mode (Photon-Counting default).
  --no-pipe          Bypass **xrtpipeline**; use existing files.
  --evt / --img      Explicit event / image filenames in --no-pipe mode.

Outputs
--------
PC mode → `pcsou.pha`  `pcbkg.pha`  `pcsou.arf`  `pcsougr1.pha`  
WT mode → `wtsou.pha`  `wtbkg.pha`  `wtsou.arf`  `wtsougr1.pha`

Requirements
-------------
  • HEASoft with Swift-XRT CALDB files installed

"""

from __future__ import annotations
import argparse, subprocess, shutil, pathlib, re, sys, os, textwrap
from typing import List, Tuple

CALDB = os.environ.get("CALDB")
TOOLS = ("xselect", "xrtmkarf", "grppha", "xrtpipeline")

# ───────────────────────── helper utilities ───────────────────────────────

def abort(msg: str):
    sys.stderr.write(f"❌  {msg}\n")
    sys.exit(1)


def check_env():
    if CALDB is None:
        abort("CALDB environment variable not set.")
    if any(shutil.which(t) is None for t in TOOLS[:3]):
        abort("HEADAS tools (xselect/xrtmkarf/grppha) not found in PATH.")


def sex2deg(tok: str) -> float:
    if ":" not in tok:
        return float(tok)
    h, m, s = [float(t) for t in tok.split(":")]
    sign = -1 if h < 0 else 1
    h = abs(h)
    deg = sign * (h + m / 60 + s / 3600)
    return deg * 15 if deg <= 24 else deg


def parse_region_center(rfile: pathlib.Path) -> Tuple[float, float]:
    pat = re.compile(r"(?:circle|annulus)\s*\(\s*([^,]+)\s*,\s*([^,]+)")
    for line in rfile.read_text().splitlines():
        m = pat.search(line)
        if m:
            return sex2deg(m.group(1)), sex2deg(m.group(2))
    abort(f"No RA/Dec in {rfile}")


def sh(cmd: List[str]):
    print("➜", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_xselect(sess: str, evt: str, reg: str, pha: str, cwd: pathlib.Path):
    script = f"""{sess}
read eve
{cwd}
{evt}
yes
filter reg {reg}
set phaname PI
extract spec
save spec {pha}
exit
no
"""
    subprocess.run(["xselect"], input=script, text=True, cwd=cwd, check=True)


def rmf_path(mode: str) -> str:
    base = pathlib.Path(CALDB) / "data/swift/xrt/cpf/rmf"
    return (base / ("swxpc0to12s6_20130101v014.rmf" if mode == "PC" else "swxwt0to2s6_20131212v015.rmf")).as_posix()


def runpipeline(obsid: str, ra: float, dec: float, mode: str, indir: pathlib.Path, outdir: pathlib.Path):
    sh([
        "xrtpipeline", "obsmode=pointing", f"datamode={mode}",
        f"srcra={ra}", f"srcdec={dec}",
        "pntra=POINT", "pntdec=POINT",
        f"indir={indir}", f"outdir={outdir}", f"steminputs={obsid}",
        "cleanup=no", "cleancols=no", "clobber=yes", "createexpomap=yes", "vigflag=yes",
    ])

# ───────────────────────── spectral extraction ────────────────────────────

def extract_spectra(evt: pathlib.Path, img: pathlib.Path, src_reg: pathlib.Path, bkg_reg: pathlib.Path, outdir: pathlib.Path, mode: str):
    evt, img = evt.resolve(), img.resolve()
    src_reg, bkg_reg = src_reg.resolve(), bkg_reg.resolve()

    prefix = "pc" if mode == "PC" else "wt"
    src_pha = outdir / f"{prefix}sou.pha"
    bkg_pha = outdir / f"{prefix}bk.pha"
    src_arf = outdir / f"{prefix}sou.arf"
    grp_pha = outdir / f"{prefix}sougr1.pha"

    run_xselect("srcsess", evt.name, src_reg.as_posix(), src_pha.as_posix(), cwd=outdir)
    run_xselect("bkgsess", evt.name, bkg_reg.as_posix(), bkg_pha.as_posix(), cwd=outdir)

    sh(["xrtmkarf", f"outfile={src_arf}", f"phafile={src_pha}", "srcx=-1", "srcy=-1", "psfflag=yes", "clobber=yes", f"expofile={img}"])

    gpscr = f"""{src_pha}
{grp_pha}
chkey respf {rmf_path(mode)}
chkey ancrf {src_arf}
chkey backf {bkg_pha}
group min 1
exit
"""
    subprocess.run(["grppha"], input=gpscr, text=True, check=True)
    print("\n✓ DONE – products written to", outdir)

# ───────────────────────── CLI  ─────────────────────────────────────────––

def cli():
    check_env()
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__)
    )
    mode_default = "PC"

    # global flags
    p.add_argument("--no-pipe", action="store_true", help="skip xrtpipeline and work from cleaned files")
    p.add_argument("--mode", default=mode_default, choices=["PC", "WT"], help="XRT datamode (default PC)")
    p.add_argument("--indir", required=False, help="input directory (raw data or cleaned event)")
    p.add_argument("--outdir", required=False, help="output directory (default ./output/OBSID or same as --indir in no-pipe mode)")

    # direct-mode specifics
    p.add_argument("--evt", help="cleaned event file (required with --no-pipe if multiple *.evt present)")
    p.add_argument("--img", help="exposure map image (required with --no-pipe if multiple *.img present)")

    # RA/Dec only used with pipeline
    p.add_argument("--src-ra", type=float, help="source RA in deg (pipeline mode)")
    p.add_argument("--src-dec", type=float, help="source Dec in deg (pipeline mode)")

    # required positionals
    p.add_argument("obsid", nargs="?", help="Swift observation ID (required unless --no-pipe)")
    p.add_argument("src_reg", help="DS9 region for source")
    p.add_argument("bkg_reg", help="DS9 region for background")

    args = p.parse_args()

    if not args.no_pipe and args.obsid is None:
        abort("OBSID is required in pipeline mode (omit it only with --no-pipe)")

    mode = args.mode.upper()
    src_reg = pathlib.Path(args.src_reg).expanduser().resolve()
    bkg_reg = pathlib.Path(args.bkg_reg).expanduser().resolve()

    # select indir / outdir defaults
    if args.no_pipe:
        if args.indir is None:
            abort("--indir is required with --no-pipe")
        indir = pathlib.Path(args.indir).expanduser().resolve()
        outdir = pathlib.Path(args.outdir).expanduser().resolve() if args.outdir else indir
    else:
        indir = pathlib.Path(args.indir).expanduser().resolve() if args.indir else pathlib.Path(f"./sources/{args.obsid}").resolve()
        default_out = pathlib.Path(f"./output/{args.obsid}").resolve()
        outdir = pathlib.Path(args.outdir).expanduser().resolve() if args.outdir else default_out

    outdir.mkdir(parents=True, exist_ok=True)

    if args.no_pipe:
        # resolve evt/img
        if args.evt and args.img:
            evt = pathlib.Path(args.evt).expanduser().resolve()
            img = pathlib.Path(args.img).expanduser().resolve()
        else:
            evt_files = sorted(indir.glob("*.evt"))
            img_files = sorted(indir.glob("*.img"))
            if len(evt_files) != 1 or len(img_files) != 1:
                abort("Specify --evt and --img when multiple event/image files exist in --indir")
            evt, img = evt_files[0], img_files[0]
        extract_spectra(evt, img, src_reg, bkg_reg, outdir, mode)
    else:
        ra = args.src_ra
        dec = args.src_dec
        if (ra is None) ^ (dec is None):
            abort("Both --src-ra and --src-dec must be provided together, or neither.")
        if ra is None:
            ra, dec = parse_region_center(src_reg)
        runpipeline(args.obsid, ra, dec, mode, indir, outdir)
        tag = "pcw" if mode == "PC" else "wtw"
        evt_files = sorted(outdir.glob(f"*{args.obsid}x{tag}*_cl.evt"))
        img_files = sorted(outdir.glob(f"*{args.obsid}x{tag}*_ex.img"))
        if len(evt_files) != 1 or len(img_files) != 1:
            abort("Could not uniquely locate pipeline products in --outdir")
        extract_spectra(evt_files[0], img_files[0], src_reg, bkg_reg, outdir, mode)

if __name__ == "__main__":
    cli()
