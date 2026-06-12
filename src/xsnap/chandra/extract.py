#!/usr/bin/env python3
"""
Pipeline module to generate Chandra images, PSF maps, auto source/background
regions, and extracted spectra from a Chandra evt2 file.

Streamline the standard CIAO workflow:
  1. Locate the user’s Conda environment.
  2. Locate the appropriate aspect solution (.asol1.fits) and mask (.msk1.fits) files.
  3. Convert the requested RA/Dec to SKY coordinates with dmcoords.
  4. Run dmcopy to make a filtered image cutout.
  5. Create a PSF map at the requested energy and encircled-energy fraction.
  6. If regions are not supplied, run wavdetect on the cutout image.
  7. Auto-create physical/SKY regions using the requested RA/Dec unless
     wavdetect finds a source within the match radius.
  8. Run specextract to build source and background spectra.

Auto-region behavior
--------------------
If --src-reg and --bkg-reg are not supplied:

  • Run wavdetect on the image cutout.
  • If a wavdetect source is found within --detect-match-radius of the input
    RA/Dec, use that detection as the center of both regions.
  • Convert the selected center RA/Dec to physical/SKY X,Y with dmcoords.
  • Write a circular physical source region with radius --src-radius.
  • Write a physical annular background region with inner/outer radii
    --bkg-inner and --bkg-outer.
  • If any wavdetect source falls inside that annulus, replace the annulus
    background with a circular physical region of radius --bkg-circle-radius.

Defaults:
  • --detect-match-radius 5 arcsec
  • --src-radius 2 arcsec
  • --bkg-inner defaults to --src-radius
  • --bkg-outer 45 arcsec
  • --bkg-circle-radius defaults to --bkg-outer

Examples
--------
HRC, auto regions:

    extract-chandra ciao-4.18 hrc_evt2.fits \
        --ra 10.83625 \
        --dec 41.30911 \
        --outdir . \
        --no-energy-filter \
        --cutout-size 512 \
        --bin-size 1

HRC, existing regions:

    extract-chandra ciao-4.18 hrc_evt2.fits \
        --src-reg source.reg \
        --bkg-reg bkg.reg \
        --outdir . \
        --no-energy-filter \
        --bin-size 1

ACIS, auto regions:

    extract-chandra ciao-4.18 acis_evt2.fits \
        --ra 10.83625 \
        --dec 41.30911 \
        --outdir . \
        --emin 500 \
        --emax 8000 \
        --cutout-size 512 \
        --bin-size 1

ACIS, existing regions:

    extract-chandra ciao-4.18 acis_evt2.fits \
        --src-reg source.reg \
        --bkg-reg bkg.reg \
        --outdir . \
        --emin 500 \
        --emax 8000 \
        --bin-size 1
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
import astropy.units as u


PHYS_REG_HEADER = (
    "# Region file format: DS9 version 4.1\n"
    "physical\n"
)


CHANDRA_SKY_PIX_ARCSEC = 0.492

def is_physical_region(reg_path):
    """
    Return True if a DS9 region file is already in physical coordinates.
    """
    text = Path(reg_path).read_text().lower()

    for line in text.splitlines():
        line = line.strip()

        if not line or line.startswith("#") or line.startswith("global"):
            continue

        return line == "physical"

    return False


def parse_simple_region(reg_path):
    """
    Parse a simple DS9 circle or annulus region.

    Supported input coordinate systems:
      • physical
      • fk5
      • icrs

    Supported shapes:
      • circle(x,y,r)
      • annulus(x,y,r_in,r_out)
    """
    lines = Path(reg_path).read_text().splitlines()

    coord_sys = None

    for line in lines:
        line = line.strip()

        if not line or line.startswith("#") or line.startswith("global"):
            continue

        low = line.lower()

        if low in {"physical", "fk5", "icrs"}:
            coord_sys = low
            continue

        m = re.match(r"(circle|annulus)\(([^)]*)\)", low)
        if not m:
            continue

        shape = m.group(1)
        vals = [v.strip().replace('"', "") for v in m.group(2).split(",")]
        vals = [float(v) for v in vals]

        return coord_sys, shape, vals

    raise RuntimeError(f"Could not parse simple circle/annulus region: {reg_path}")


def convert_region_to_physical(reg_path, outdir, evt2, asol, suffix):
    """
    Convert a simple celestial DS9 region to physical/SKY coordinates.

    If the input region is already physical, return it unchanged.
    """
    reg_path = Path(reg_path).resolve()

    if is_physical_region(reg_path):
        print(f"Region already physical: {reg_path}")
        return reg_path

    coord_sys, shape, vals = parse_simple_region(reg_path)

    if coord_sys not in {"fk5", "icrs"}:
        raise RuntimeError(
            f"Unsupported region coordinate system in {reg_path}: {coord_sys}"
        )

    ra = vals[0]
    dec = vals[1]

    x, y = get_sky_xy(evt2, asol, ra, dec)

    out_reg = Path(outdir) / f"{reg_path.stem}_physical_{suffix}.reg"

    if shape == "circle":
        radius_pix = arcsec_to_sky_pix(vals[2])

        out_reg.write_text(
            PHYS_REG_HEADER
            + f"circle({x:.3f},{y:.3f},{radius_pix:.3f})\n"
        )

    elif shape == "annulus":
        inner_pix = arcsec_to_sky_pix(vals[2])
        outer_pix = arcsec_to_sky_pix(vals[3])

        out_reg.write_text(
            PHYS_REG_HEADER
            + f"annulus({x:.3f},{y:.3f},{inner_pix:.3f},{outer_pix:.3f})\n"
        )

    else:
        raise RuntimeError(f"Unsupported region shape in {reg_path}: {shape}")

    print(f"Converted region to physical: {out_reg}")
    return out_reg

def get_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )

    p.add_argument("env", help="Expected conda environment name")
    p.add_argument("evt2_file", help="Chandra evt2 file")

    p.add_argument("--src-reg", default=None, help="Existing source region")
    p.add_argument("--bkg-reg", default=None, help="Existing background region")
    p.add_argument("--outdir", default=None, help="Output directory")

    p.add_argument("--ra", type=float, default=None, help="Source RA in degrees")
    p.add_argument("--dec", type=float, default=None, help="Source Dec in degrees")

    p.add_argument("--src-radius", type=float, default=2.0,
                   help="Auto source radius in arcsec")
    p.add_argument("--bkg-inner", type=float, default=None,
                   help="Auto background annulus inner radius in arcsec")
    p.add_argument("--bkg-outer", type=float, default=45.0,
                   help="Auto background annulus outer radius in arcsec")
    p.add_argument("--bkg-circle-radius", type=float, default=None,
               help="Fallback circular background radius in arcsec")
    p.add_argument("--detect-match-radius", type=float, default=5,
                   help="Use wavdetect source as center if within this radius, arcsec")

    p.add_argument("--wav-scales", default="1 2 4 8",
                   help="wavdetect wavelet scales in image pixels")
    p.add_argument("--wav-sigthresh", type=float, default=1e-6,
                   help="wavdetect significance threshold")

    p.add_argument("--cutout-size", type=int, default=512,
                   help="Image cutout size in SKY pixels")
    p.add_argument("--bin-size", type=float, default=1.0,
                   help="Image bin size in SKY pixels")

    p.add_argument("--emin", type=int, default=500)
    p.add_argument("--emax", type=int, default=8000)
    p.add_argument("--ccd-id", type=int, default=None)
    p.add_argument("--no-energy-filter", action="store_true")

    p.add_argument("--psf-energy", type=float, default=1.4967)
    p.add_argument("--ecf", type=float, default=0.90)

    return p.parse_args()


def run(cmd):
    print(">>>", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True)


def warn_env(expected):
    current = os.environ.get("CONDA_DEFAULT_ENV", "")
    if current != expected:
        print(f"Warning: requested env '{expected}', but CONDA_DEFAULT_ENV={current}")


def find_unique(pattern, directories):
    matches = []
    for d in directories:
        d = Path(d)
        if d.exists():
            matches.extend(sorted(d.glob(pattern)))

    if len(matches) != 1:
        msg = "\n".join(str(m) for m in matches)
        sys.exit(
            f"ERROR: expected exactly one {pattern!r}, found {len(matches)}.\n{msg}"
        )

    return matches[0]


def arcsec_to_sky_pix(value_arcsec):
    """
    Convert arcsec to native Chandra SKY pixels.

    Native Chandra SKY pixels are approximately 0.492 arcsec/pixel.
    """
    return value_arcsec / CHANDRA_SKY_PIX_ARCSEC


def get_sky_xy(evt2, asol, ra, dec):
    """
    Run dmcoords and parse SKY(X,Y) from verbose output.

    This is used for both the requested RA/Dec and any wavdetect-selected
    source center, so all auto regions are written in physical/SKY coordinates.
    """
    cmd = [
        "dmcoords",
        f"infile={evt2}",
        f"asolfile={asol}",
        "option=cel",
        "celfmt=deg",
        f"ra={ra}",
        f"dec={dec}",
        "verbose=1",
    ]

    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)

    m = re.search(
        r"SKY\(X,Y\):\s*([0-9.+-]+)\s+([0-9.+-]+)",
        proc.stdout,
    )

    if not m:
        print(proc.stdout)
        print(proc.stderr)
        raise RuntimeError("Could not parse SKY(X,Y) from dmcoords output.")

    x = float(m.group(1))
    y = float(m.group(2))

    print(f"RA={ra:.8f}, Dec={dec:.8f} -> SKY(X,Y)=({x:.2f}, {y:.2f})")
    return x, y


def make_image_slice(args, instrument, x, y):
    half = args.cutout_size / 2

    xmin = int(x - half)
    xmax = int(x + half)
    ymin = int(y - half)
    ymax = int(y + half)

    filters = ""

    use_energy = not args.no_energy_filter
    is_hrc = "hrc" in (instrument or "").lower()

    if is_hrc and use_energy:
        print("HRC detected; disabling energy filter because HRC evt2 has no ENERGY column.")
        use_energy = False

    ccd_id = args.ccd_id

    if ccd_id is None and use_energy:
        ccd_id = 7

    if ccd_id is not None:
        filters += f"[ccd_id={ccd_id}]"

    if use_energy:
        filters += f"[energy={args.emin}:{args.emax}]"

    image_bin = (
        f"[bin x={xmin}:{xmax}:{args.bin_size},"
        f"y={ymin}:{ymax}:{args.bin_size}]"
    )

    return filters + image_bin


def run_wavdetect(img_out, psf_out, outdir, args):
    """
    Run CIAO wavdetect on the cutout image.

    The FITS source list is used only for source-position decisions. The final
    extraction regions are written separately as physical/SKY regions.
    """
    src_list = outdir / "wavdetect_src.fits"
    scell = outdir / "wavdetect_scell.fits"
    recon = outdir / "wavdetect_recon.fits"
    nbkg = outdir / "wavdetect_nbkg.fits"
    reg = outdir / "wavdetect.reg"

    run([
        "wavdetect",
        f"infile={img_out}",
        f"outfile={src_list}",
        f"scellfile={scell}",
        f"imagefile={recon}",
        f"defnbkgfile={nbkg}",
        f"regfile={reg}",
        f"psffile={psf_out}",
        f"scales={args.wav_scales}",
        f"sigthresh={args.wav_sigthresh}",
        "clobber=yes",
    ])

    if not src_list.exists():
        print("No wavdetect source list found; using input RA/Dec for regions.")
        return []

    tab = Table.read(src_list)

    if "RA" not in tab.colnames or "DEC" not in tab.colnames:
        print("wavdetect source list has no RA/DEC columns; using input RA/Dec.")
        return []

    detections = []

    for row in tab:
        coord = SkyCoord(float(row["RA"]), float(row["DEC"]), unit="deg", frame="icrs")
        detections.append(coord)

    print(f"wavdetect detections read: {len(detections)}")
    return detections

def make_full_image_slice(args, instrument):
    filters = ""

    use_energy = not args.no_energy_filter
    is_hrc = "hrc" in (instrument or "").lower()

    if is_hrc and use_energy:
        print("HRC detected; disabling energy filter because HRC evt2 has no ENERGY column.")
        use_energy = False

    ccd_id = args.ccd_id

    if ccd_id is None and use_energy:
        ccd_id = 7

    if ccd_id is not None:
        filters += f"[ccd_id={ccd_id}]"

    if use_energy:
        filters += f"[energy={args.emin}:{args.emax}]"

    filters += f"[bin sky={args.bin_size}]"

    return filters


def choose_region_center(input_coord, detections, match_radius_arcsec):
    """
    Use the closest wavdetect source as the region center if one or more
    detections fall within the requested matching radius. Otherwise keep the
    input RA/Dec.
    """
    if not detections:
        print("No wavdetect detections; using input RA/Dec as region center.")
        return input_coord, False

    matches = []

    for det in detections:
        sep = input_coord.separation(det).arcsec

        if sep <= match_radius_arcsec:
            matches.append((sep, det))

    if matches:
        matches.sort(key=lambda item: item[0])
        best_sep, center = matches[0]

        print(
            "Using closest wavdetect source within match radius as region center: "
            f"RA={center.ra.deg:.8f}, Dec={center.dec.deg:.8f}, "
            f"offset={best_sep:.2f}\", "
            f"matches={len(matches)}"
        )

        return center, True

    nearest_sep = min(input_coord.separation(det).arcsec for det in detections)

    print(
        "No wavdetect source within "
        f"{match_radius_arcsec:.2f}\"; nearest detection is "
        f"{nearest_sep:.2f}\" away. Using input RA/Dec as region center."
    )

    return input_coord, False


def annulus_contains_detection(center_coord, detections, inner_arcsec, outer_arcsec):
    """
    Return True if any wavdetect source falls inside the proposed background
    annulus.
    """
    for det in detections:
        sep = center_coord.separation(det).arcsec

        if inner_arcsec <= sep <= outer_arcsec:
            print(
                "Detected source inside background annulus at "
                f"{sep:.2f}\" from center; using circular background instead."
            )
            return True

    return False


def write_physical_regions(
    outdir,
    x,
    y,
    src_radius_arcsec,
    bkg_inner_arcsec,
    bkg_outer_arcsec,
    bkg_circle_radius_arcsec,
    use_circle_bkg,
):
    """
    Write CIAO-safe physical/SKY source and background regions.

    Region center is in SKY X,Y pixels. Input radii are given in arcsec and
    converted to native Chandra SKY pixels before writing.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    src_radius_pix = arcsec_to_sky_pix(src_radius_arcsec)
    bkg_inner_pix = arcsec_to_sky_pix(bkg_inner_arcsec)
    bkg_outer_pix = arcsec_to_sky_pix(bkg_outer_arcsec)
    bkg_circle_radius_pix = arcsec_to_sky_pix(bkg_circle_radius_arcsec)

    src_reg = outdir / "source.reg"
    bkg_reg = outdir / "bkg.reg"

    src_reg.write_text(
        PHYS_REG_HEADER
        + f"circle({x:.3f},{y:.3f},{src_radius_pix:.3f})\n"
    )

    if use_circle_bkg:
        bkg_reg.write_text(
            PHYS_REG_HEADER
            + f"circle({x:.3f},{y:.3f},{bkg_circle_radius_pix:.3f})\n"
        )
    else:
        bkg_reg.write_text(
            PHYS_REG_HEADER
            + f"annulus({x:.3f},{y:.3f},{bkg_inner_pix:.3f},{bkg_outer_pix:.3f})\n"
        )

    print(f"Wrote physical source region: {src_reg}")
    print(f"Wrote physical background region: {bkg_reg}")

    return src_reg, bkg_reg


def write_auto_regions(outdir, evt2, asol, args, detections):
    """
    Create physical/SKY source and background regions.

    wavdetect is used to decide whether the center should be shifted from the
    input RA/Dec. dmcoords is then used to convert the final center to SKY X,Y.
    """
    input_coord = SkyCoord(args.ra, args.dec, unit="deg", frame="icrs")

    center_coord, used_detection = choose_region_center(
        input_coord,
        detections,
        args.detect_match_radius,
    )

    center_x, center_y = get_sky_xy(
        evt2,
        asol,
        center_coord.ra.deg,
        center_coord.dec.deg,
    )

    use_circle_bkg = annulus_contains_detection(
        center_coord,
        detections,
        args.bkg_inner,
        args.bkg_outer,
    )

    src_reg, bkg_reg = write_physical_regions(
        outdir=outdir,
        x=center_x,
        y=center_y,
        src_radius_arcsec=args.src_radius,
        bkg_inner_arcsec=args.bkg_inner,
        bkg_outer_arcsec=args.bkg_outer,
        bkg_circle_radius_arcsec=args.bkg_circle_radius,
        use_circle_bkg=use_circle_bkg,
    )

    if used_detection:
        print("Auto regions were centered on the nearest matching wavdetect source.")
    else:
        print("Auto regions were centered on the input RA/Dec.")

    return src_reg, bkg_reg


def main():
    args = get_args()
    if args.bkg_inner is None:
        args.bkg_inner = args.src_radius
    if args.bkg_circle_radius is None:
        args.bkg_circle_radius = args.bkg_outer
    warn_env(args.env)
    
    using_auto_regions = args.src_reg is None or args.bkg_reg is None

    if using_auto_regions and (args.ra is None or args.dec is None):
        sys.exit("ERROR: --ra and --dec are required when auto-generating regions.")

    evt2 = Path(args.evt2_file).resolve()
    if not evt2.exists():
        sys.exit(f"ERROR: {evt2} not found")

    primary = evt2.parent
    secondary = primary.parent / "secondary"

    search_dirs = [primary, secondary, primary.parent]

    asol = find_unique("*asol1.fits*", search_dirs)
    msk = find_unique("*msk1.fits*", search_dirs)

    outdir = Path(args.outdir).resolve() if args.outdir else primary
    outdir.mkdir(parents=True, exist_ok=True)

    with fits.open(evt2) as hdul:
        hdr = hdul[0].header
        telescope = hdr.get("TELESCOP", "")
        instrument = hdr.get("INSTRUME", "")

    print(f"Telescope: {telescope}")
    print(f"Instrument: {instrument}")
    print(f"ASOL: {asol}")
    print(f"MASK: {msk}")

    if args.ra is not None and args.dec is not None:
        input_x, input_y = get_sky_xy(evt2, asol, args.ra, args.dec)
        image_slice = make_image_slice(args, instrument, input_x, input_y)
    else:
        image_slice = make_full_image_slice(args, instrument)

    img_out = outdir / "evt_img.fits"
    psf_out = outdir / "evt_psfmap.fits"

    epoch = evt2.parent.parent.name
    outroot = outdir / f"spec{epoch}"

    run([
        "dmcopy",
        f"{evt2}{image_slice}",
        img_out,
        "clobber=yes",
    ])

    run([
        "mkpsfmap",
        img_out,
        f"outfile={psf_out}",
        f"energy={args.psf_energy}",
        f"ecf={args.ecf}",
        "clobber=yes",
    ])

    if args.src_reg is None or args.bkg_reg is None:
        detections = run_wavdetect(img_out, psf_out, outdir, args)

        src_reg, bkg_reg = write_auto_regions(
            outdir,
            evt2,
            asol,
            args,
            detections,
        )
    else:
        src_reg = Path(args.src_reg).resolve()
        bkg_reg = Path(args.bkg_reg).resolve()

        if not src_reg.exists():
            sys.exit(f"ERROR: source region not found: {src_reg}")
        if not bkg_reg.exists():
            sys.exit(f"ERROR: background region not found: {bkg_reg}")

        src_reg = convert_region_to_physical(
            src_reg,
            outdir,
            evt2,
            asol,
            "src",
        )

        bkg_reg = convert_region_to_physical(
            bkg_reg,
            outdir,
            evt2,
            asol,
            "bkg",
        )

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
        "binspec=1",
        "clobber=yes",
    ])

    print("\nDONE — products in", outdir)
    print(" ", img_out.name)
    print(" ", psf_out.name)
    print(" ", "wavdetect_src.fits")
    print(" ", f"spec{epoch}_grp.pi")
    print(" ", src_reg.name)
    print(" ", bkg_reg.name)


def cli() -> int | None:
    return main()


if __name__ == "__main__":
    main()