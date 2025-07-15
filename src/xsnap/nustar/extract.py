#!/usr/bin/env python3
"""
extract_nustar.py – one-command NuSTAR spectrum extractor
=========================================================

Steps performed
---------------
1. **nupipeline**   – calibrate & clean the observation\n
2. **nuproducts**   – extract spectra / ARFs for FPMA and FPMB

Usage
-----
extract-nustar  OBSID  SRC.REG  BKG.REG
               [--indir  DIR]   (default: ./sources/OBSID)
               [--outdir DIR]   (default: ./output/OBSID)
               [--prod   DIR]   (default: ./products/OBSID)
               [--ra RA --dec DEC]

Positional arguments
--------------------
OBSID       11-digit NuSTAR observation ID (e.g. 80902505002)  
SRC.REG     Source region file (WCS, FK5)  
BKG.REG     Background region file (WCS, FK5)

Optional arguments
------------------
--indir DIR     directory containing the unfiltered event files  
--outdir DIR    directory where **nupipeline** writes the cleaned data  
--prod DIR      root directory where **nuproducts** writes FPMA/ FPMB products  
--ra  RA        override RA  (deg) fed to *nupipeline*  
--dec DEC       override Dec (deg) fed to *nupipeline*

If --ra or --dec are not supplied, the centre of `SRC.REG`
(first circle/annulus) is used.
"""
from __future__ import annotations
import argparse, subprocess, shutil, pathlib, re, sys, os, textwrap
from typing import Tuple

TOOLS = ("nupipeline", "nuproducts")
CALDB = os.environ.get("CALDB")

# ───────────────────────── helpers ──────────────────────────────────────────
def parse_region_center(rfile: pathlib.Path) -> Tuple[float, float]:
    """
    Return (RA,Dec) in deg from the FIRST circle/annulus found in a DS9
    region file (FK5/ICRS).  Works for lines like

        circle(210.9106746,54.3116511,0.00346")
        annulus ( 05:34:32.0 , -69:33:33 , 30" , 45" )

    Sexagesimal values are auto-converted to degrees.
    """
    def sex2deg(token: str) -> float:
        if ":" not in token:
            return float(token)
        parts = [float(x) for x in token.split(":")]
        if len(parts) == 3:                       # h:m:s  or  d:m:s
            sign = -1 if parts[0] < 0 else 1
            parts[0] = abs(parts[0])
            val = parts[0] + parts[1]/60 + parts[2]/3600
            return sign * val * (15 if token.count(":") == 2 and
                                 token.strip().lstrip("+-").count(":") == 2
                                 and val <= 24 else 1)
        raise ValueError("Unrecognised sexagesimal token")

    pat = re.compile(r"(circle|annulus)\s*\(\s*([^)]*)\)")
    with open(rfile) as fh:
        for ln in fh:
            m = pat.search(ln)
            if m:
                inside = m.group(2)
                tokens = [t.strip() for t in inside.split(",")]
                if len(tokens) < 2:
                    break
                ra_deg  = sex2deg(tokens[0])
                dec_deg = sex2deg(tokens[1])
                return ra_deg, dec_deg
    raise RuntimeError(f"No circle/annulus line with RA/Dec in {rfile}")

def sh(cmd: list[str], **kw):
    print("➜", " ".join(cmd))
    subprocess.run(cmd, check=True, **kw)

def runpipeline(obsid, ra, dec, indir, outdir):
    stem = f"nu{obsid}"
    sh([
        "nupipeline",
        f"indir={indir}", f"outdir={outdir}",
        f"steminputs={stem}", 
        f"srcra={ra}", f"srcdec={dec}",
    ])

def runproducts(obsid, instrument, outdir, prodroot):
    stem = f"nu{obsid}"
    prodroot.mkdir(parents=True, exist_ok=True)
    sh([
        "nuproducts",
        f"indir={outdir}", f"outdir={prodroot}",
        f"steminputs={stem}", f"instrument={instrument}",
        f"srcregionfile={outdir/'source.reg'}",
        f"bkgregionfile={outdir/'bkg.reg'}",
        "bkgextract=yes", "clobber=yes"
    ])

# ───────────────────────── core routine ─────────────────────────────────────
def extract_nustar(obsid, src_reg, bkg_reg,
                   *, indir, outdir, prodroot,
                   ra=None, dec=None, run_pipe=True):
    if any(shutil.which(t) is None for t in TOOLS):
        sys.exit("❌  NuSTAR FTOOLS not found.")
    if CALDB is None:
        sys.exit("❌  CALDB not set")

    indir, outdir, prodroot = map(lambda p: p.expanduser().resolve(),
                                  (indir, outdir, prodroot))
    outdir.mkdir(parents=True, exist_ok=True)

    # copy regions into outdir with canonical names
    (outdir / "source.reg").write_text(pathlib.Path(src_reg).read_text())
    (outdir / "bkg.reg").write_text(pathlib.Path(bkg_reg).read_text())

    if ra is None or dec is None:
        ra, dec = parse_region_center(src_reg)
        print("(RA,Dec) taken from region:", f"{ra:.6f} {dec:.6f}")
    else:
        print("(RA,Dec) supplied:", f"{ra:.6f} {dec:.6f}")

    if run_pipe:
        runpipeline(obsid, ra, dec, indir, outdir)

    for inst in ("FPMA", "FPMB"):
        runproducts(obsid, inst, outdir, prodroot / inst)

    print("\n✓ DONE – products under", prodroot)

# ───────────────────────── CLI entry-point ──────────────────────────────────
def cli() -> None:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__)
    )
    ap.add_argument("obsid")
    ap.add_argument("src_reg")
    ap.add_argument("bkg_reg")
    ap.add_argument("--indir",  help="input dir  (default ./sources/OBSID)")
    ap.add_argument("--outdir", help="pipeline out dir (default ./output/OBSID)")
    ap.add_argument("--prod",   help="products root dir (default ./products/OBSID)")
    ap.add_argument("--ra",  type=float, help="source RA  deg (override region)")
    ap.add_argument("--dec", type=float, help="source Dec deg (override region)")
    ap.add_argument("--no-pipe", action="store_true",
                    help="skip nupipeline (use existing cleaned data)")
    ns = ap.parse_args()

    obs = ns.obsid
    indir  = pathlib.Path(ns.indir or f"./sources/{obs}")
    outdir = pathlib.Path(ns.outdir or f"./output/{obs}")
    prod   = pathlib.Path(ns.prod  or f"./products/{obs}")

    extract_nustar(obs, pathlib.Path(ns.src_reg), pathlib.Path(ns.bkg_reg),
                   indir=indir, outdir=outdir, prodroot=prod,
                   ra=ns.ra, dec=ns.dec, run_pipe=not ns.no_pipe)

if __name__ == "__main__":
    cli()