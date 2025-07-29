#!/usr/bin/env python3
"""
Pipeline module to generate 500–8000 eV images, PSF maps, and extracted spectra from a Chandra evt2 file.

Streamline the standard CIAO workflow:
  1. Locate the user’s Conda environment.
  2. Locate the appropriate aspect solution (.asol1.fits) and mask (.msk1.fits) files.
  3. Run dmcopy to filter CCD 7 events in 500–8000 eV and bin to a small image.
  4. Create a PSF map at 1.4967 keV with 90% encircled-energy fraction.
  5. Run specextract to build source and background spectra (NUM_CTS grouping), applying correct weights.

Usage:
    extract-chandra <COND ENV> <EVT2_FILE> <SRC_REG> <BKG_REG> [OUTDIR]

Positional arguments:
  COND_ENV     Name of the Conda environment to check against $CONDA_DEFAULT_ENV.
  EVT_FILE    Path to the Chandra level-2 event file.
  SRC_REG      Source region (filename or literal sky region).
  BKG_REG      Background region (filename or literal sky region).
  OUTDIR       Optional directory for outputs (defaults to EVT2_FILE’s directory).

Requirements:
  • CIAO (including dmcopy, mkpsfmap, specextract) installed.
  • CONDA_DEFAULT_ENV set to the requested environment name.

"""
import sys, os, subprocess
from pathlib import Path
from shutil import which

def warn_env(name):
    cur = os.environ.get("CONDA_DEFAULT_ENV","")
    if cur != name:
        print(f"Warning: you asked for env '{name}', but CONDA_DEFAULT_ENV={cur}")

def find_unique(globstr, where):
    L = list(Path(where).glob(globstr))
    if len(L)!=1:
        sys.exit(f"ERROR: expected exactly one {globstr!r} in {where}, found {len(L)}")
    return str(L[0])

def run(cmd, **kw):
    print(">>>", " ".join(cmd))
    subprocess.run(cmd, check=True, **kw)
    
def main():
    argv = sys.argv
    if len(argv) == 2 and argv[1] in ("-h","--help", "-help"):
        print(__doc__)
        sys.exit(0)
    if len(argv) not in (5,6):
        sys.exit("Usage: script.py [-h] ENV EVT2_FILE SRC_REG BKG_REG [OUTDIR]")

    _, env, evt2_path, src_reg, bkg_reg, *rest = argv
    outdir = rest[0] if rest else None

    warn_env(env)

    evt2 = Path(evt2_path).resolve()
    if not evt2.exists():
        sys.exit(f"ERROR: {evt2} not found")

    primary   = evt2.parent
    secondary = primary.parent / "secondary"

    asol = find_unique("*asol1.fits*", primary)
    msk  = find_unique("*msk1.fits*",  secondary)

    # determine outdir
    if outdir:
        outdir = Path(outdir)
    else:
        outdir = primary
    outdir.mkdir(parents=True, exist_ok=True)

    # derive base names
    prefix = "evt_500_8000"
    
    # derive epoch from the parent folder name
    epoch = evt2.parent.parent.name

    # build your new outroot
    outroot = str(outdir / f"spec{epoch}")

    # 1) dmcopy
    img_out = outdir / f"{prefix}_img.fits"
    spec_slice = "[ccd_id=7][energy=500:8000]" + \
                 "[bin x=3848:4360:1.0,y=3860:4372:1.0]"
    run([
        "dmcopy",
        f"{str(evt2)}{spec_slice}",
        str(img_out)
    ])

    # 2) mkpsfmap
    psf_out = outdir / f"{prefix}_psfmap.fits"
    run([
        "mkpsfmap",
        str(img_out),
        f"outfile={psf_out}",
        "energy=1.4967",
        "ecf=0.90"
    ])

    # 3) specextract
    run([
        "specextract",
        f"infile={evt2}[sky=region({src_reg})]",
        f"outroot={outroot}",
        f"bkgfile={evt2}[sky=region({bkg_reg})]",
        "weight=no",
        "correct=yes",
        f"asp={asol}",
        f"mskfile={msk}",
        "grouptype=NUM_CTS",
        "binspec=1"
    ])

    print("\nDONE — products in", outdir)
    print("  ", img_out.name)
    print("  ", psf_out.name)
    print("  ", f"spec{epoch}_grp.pi")
    
def cli() -> int | None:       
    return main() 

if __name__=="__main__":
    main()