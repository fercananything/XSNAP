#!/usr/bin/env python3
"""
Pipeline module to extract and group XMM-Newton EPIC (PN, MOS1, MOS2) spectra.

Automates the standard SAS (Science Analysis System) workflow:
  1. Build CCF, ingest the ODF, and set SAS environment variables.
  2. Locate the PPS directory alongside the ODF and prepare event files.
  3. Extract source and background spectra.
  4. Generate response files (RMF/ARF) and group spectra for PN, MOS1, and MOS2.
  5. Move final grouped spectra to the specified output directory.

Usage:
    extract-xmm <ODF_DIR> <source.reg> <PNbkg.reg> [MOS1bkg.reg] [MOS2bkg.reg] [OUTDIR]

Positional arguments:
  ODF_DIR        Path to the XMM-Newton ODF directory.
  source.reg     Source region in physical coordinates (filename or literal physical “circle(x,y,r)”).
  PNbkg.reg      PN background region in physical coordinates (filename or literal).
  MOS1bkg.reg    Optional MOS1 background region in physical coordinates (defaults to PNbkg.reg if omitted).
  MOS2bkg.reg    Optional MOS2 background region in physical coordinates (defaults to PNbkg.reg if omitted).
  OUTDIR         Destination for grouped spectra (default: PPS directory).

Requirements:
  - SAS (Science Analysis System) must be installed.
  - PPS and ODF directory must have the same parent directory

"""
from __future__ import annotations
import argparse, subprocess, os, re, shutil, sys, textwrap
from pathlib import Path
from contextlib import contextmanager
from typing import Union


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: str, cwd: Path|None = None):
    """Run *cmd* in Bash, streaming output and raising on non‑zero status."""
    print(f"\n>> {cmd}\n", flush=True)
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, cwd=cwd)

@contextmanager
def cd(path: Path):
    """Temporarily ``chdir`` to *path*."""
    prev = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)

def parse_region(arg: Union[str, Path]) -> str:
    """
    Return a DS9 region string for either circle(...) or annulus(...), by
    reading *arg* as a file or literal.

    If *arg* is a path to an existing file, this reads it line by line,
    strips off any “#…” comments, and returns the first line matching
    either “circle(…)” or “annulus(…)”.

    Otherwise if *arg* itself (as str) starts with one of those shapes,
    it’s returned (comments are not allowed in-line for literals).

    Raises:
        ValueError: if no matching shape is found or *arg* is neither a file
                    nor a supported literal.
    """
    # compile once
    shape_re = re.compile(r'^\s*(circle|annulus)\s*\(.*\)', re.IGNORECASE)

    p = Path(arg).expanduser()
    if p.exists():
        # region file mode
        for raw in p.read_text().splitlines():
            line = raw.split('#', 1)[0].strip()
            if shape_re.match(line):
                return line
        raise ValueError(f"No circle(...) or annulus(...) found in {p}")

    # literal mode
    s = str(arg).strip()
    if shape_re.match(s):
        return s

    raise ValueError(f"{arg!r} is neither a region file nor a supported literal")

# ---------------------------------------------------------------------------
# Main recipe
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=textwrap.dedent(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('odf_dir',            type=Path)
    parser.add_argument('src_reg',            type=str)
    parser.add_argument('pn_bkg_reg',         type=str)
    parser.add_argument('mos1_bkg_reg', nargs='?', default=None, type=str)
    parser.add_argument('mos2_bkg_reg', nargs='?', default=None, type=str)
    parser.add_argument('outdir',       nargs='?', default=None, type=Path)
    args = parser.parse_args()

    odf_dir: Path = args.odf_dir.expanduser().resolve()
    if not odf_dir.is_dir():
        sys.exit(f"ODF directory not found: {odf_dir}")

    # region strings -------------------------------------------------------
    src_reg  = parse_region(args.src_reg)
    pn_bkg   = parse_region(args.pn_bkg_reg)
    mos1_bkg = parse_region(args.mos1_bkg_reg or args.pn_bkg_reg)
    mos2_bkg = parse_region(args.mos2_bkg_reg or args.pn_bkg_reg)

    # --------------------------------------------------------------------
    # 1. Build CCF & ingest ODF
    # --------------------------------------------------------------------
    os.environ['SAS_ODF'] = str(odf_dir)
    with cd(odf_dir):
        run('cifbuild')
        os.environ['SAS_CCF'] = str(odf_dir / 'ccf.cif')
        run('odfingest')
        sum_sas = next(Path('.').glob('*SUM.SAS'))
        os.environ['SAS_ODF'] = str(sum_sas.resolve())

    # --------------------------------------------------------------------
    # 2. Locate PPS directory (sibling of ODF)
    # --------------------------------------------------------------------
    pps = (odf_dir.parent / 'pps').resolve()
    if not pps.is_dir():
        sys.exit(f'pps directory not found next to ODF: {pps}')

    pn_evt = next(pps.glob('*PN*EVLI*FTZ'))
    m1_evt = next(pps.glob('*M1*EVLI*FTZ'))
    m2_evt = next(pps.glob('*M2*EVLI*FTZ'))
    # use friendly names for subsequent SAS commands
    shutil.copy2(pn_evt, pps/'PN_evt.FTZ')
    shutil.copy2(m1_evt, pps/'M1_evt.FTZ')
    shutil.copy2(m2_evt, pps/'M2_evt.FTZ')

    # --------------------------------------------------------------------
    # 3. Work inside PPS directory
    # --------------------------------------------------------------------
    with cd(pps):
        # create images (optional)
        run('evselect table=PN_evt.FTZ  imagebinning=binSize imageset=PNimage.fits  withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80')
        run('evselect table=M1_evt.FTZ imagebinning=binSize imageset=MOS1image.fits withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80')
        run('evselect table=M2_evt.FTZ imagebinning=binSize imageset=MOS2image.fits withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80')

        # helper to extract and group spectra
        def extract(inst: str, evt: str, src: str, bkg: str, pattern: str, max_chan: str):
            src_spec = f"{inst}source_spectrum.fits"
            bkg_spec = f"{inst}bkg_spectrum.fits"
            rmf      = f"{inst}.rmf"
            arf      = f"{inst}.arf"
            grp      = f"{inst}_spectrum_grp.fits"
            run(f"evselect table={evt} withspectrumset=yes spectrumset={src_spec} energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={max_chan} expression='(FLAG==0) && ({pattern}) && ((X,Y) IN {src})'")
            run(f"evselect table={evt} withspectrumset=yes spectrumset={bkg_spec} energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={max_chan} expression='(FLAG==0) && ({pattern}) && ((X,Y) IN {bkg})'")
            run(f"backscale spectrumset={src_spec} badpixlocation={evt}")
            run(f"backscale spectrumset={bkg_spec} badpixlocation={evt}")
            run(f"rmfgen spectrumset={src_spec} rmfset={rmf}")
            run(f"arfgen spectrumset={src_spec} arfset={arf} withrmfset=yes rmfset={rmf} badpixlocation={evt} detmaptype=psf")
            run(f"specgroup spectrumset={src_spec} mincounts=1 rmfset={rmf} arfset={arf} backgndset={bkg_spec} groupedset={grp}")
            return grp

        pn_grp = extract('PN', 'PN_evt.FTZ', src_reg, pn_bkg,
                         'PATTERN<=4', '20479')
        m1_grp = extract('M1', 'M1_evt.FTZ', src_reg, mos1_bkg,
                         '#XMMEA_EM && (PATTERN<=12)', '11999')
        m2_grp = extract('M2', 'M2_evt.FTZ', src_reg, mos2_bkg,
                         '#XMMEA_EM && (PATTERN<=12)', '11999')

    # --------------------------------------------------------------------
    # 4. Move grouped spectra to output directory
    # --------------------------------------------------------------------
    outdir = (args.outdir or pps).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    for grp in (pn_grp, m1_grp, m2_grp):
        shutil.move(str(pps/grp), outdir/grp)
    print("\nDone – grouped spectra written to", outdir)

def cli() -> int | None:       
    return main()      

if __name__ == "__main__":
    main()
