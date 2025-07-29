#!/usr/bin/env python3
"""
Pipeline module to run Swift/XRT processing in PC Mode across multiple epochs and stack time-binned products.

Automates the stacking workflow:
  1. Executes xrtpipeline for each epoch listed in a text file, producing calibrated event and image files.
  2. Bins the resulting event (.evt) and image (.img) files into time intervals (days since explosion) of user-defined size.
  3. Stacks events and images in each bin into dedicated subdirectories under `stacks/`, named by their first–last epochs.
  4. Updates the stacked event files with MJD-BEG and MJD-END headers.

Usage:
    swift-stack-pc \
      --epochs source_epochs.txt \
      --bin-size 2 \
      --explosion-mjd 59000.5 \
      --src-ra <RA> \
      --src-dec <Dec> \
      [--pn-tra POINT] [--pn-dec POINT]

Arguments:
  --epochs         File listing one epoch/obsid identifier per line.
  --bin-size       Time-bin width in days (default: 1.0).
  --explosion-mjd  Reference MJD of explosion.
  --src-ra         Source RA (deg) for pipeline processing.
  --src-dec        Source Dec (deg) for pipeline processing.
  --pn-tra         (Optional) PN attitude RA override (default: POINT).
  --pn-dec         (Optional) PN attitude Dec override (default: POINT).

Requirements:
  • HEASoft with Swift-XRT CALDB files installed

"""
import os
import subprocess
import argparse, textwrap
from astropy.io import fits
from astropy.time import Time


def run_pipeline(epoch, src_ra, src_dec, pn_tra, pn_dec):
    indir = f"./sources/{epoch}/"
    outdir = f"./output/out_{epoch}"
    stem = epoch
    cmd = [
        "xrtpipeline",
        "obsmode=pointing", "datamode=PC",
        f"srcra={src_ra}", f"srcdec={src_dec}",
        f"pntra={pn_tra}", f"pntdec={pn_dec}",
        f"indir={indir}", f"outdir={outdir}",
        f"steminputs={stem}",
        "cleanup=no", "cleancols=no", "clobber=yes",
        "createexpomap=yes", "vigflag=yes"
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return outdir


def get_date_obs_mjd(evtfile):
    hdr = fits.getheader(evtfile, ext=1)
    date_obs = hdr.get('DATE-OBS')
    timesys = hdr.get('TIMESYS').lower()
    if not date_obs:
        raise ValueError(f"DATE-OBS not found in {evtfile}")
    t = Time(date_obs, format='isot', scale=timesys)
    return t.mjd


def stack_events(evt_list, out_evt):
    txt = "elist.txt"
    with open(txt, 'w') as f:
        for fn in evt_list:
            f.write(fn + "\n")
    cmd = [
        "extractor",
        f"filename=@{txt}",
        f"eventsout={out_evt}",
        "imgfile=NONE", "phafile=NONE", "fitsbinlc=NONE",
        "regionfile=NONE", "timefile=NONE",
        "xcolf=x", "ycolf=y", "tcol=time"
    ]
    print(f"Stacking {len(evt_list)} event files into {out_evt}")
    subprocess.run(cmd, check=True)


def stack_images(img_list, out_img):
    if len(img_list) < 2:
        os.rename(img_list[0], out_img)
        return
    tmp = "tmp_sum_0.img"
    subprocess.run(["farith", img_list[0], img_list[1], tmp, "+"], check=True)
    for i, img in enumerate(img_list[2:], start=2):
        next_tmp = f"tmp_sum_{i-1}.img"
        subprocess.run(["farith", tmp, img, next_tmp, "+"], check=True)
        os.remove(tmp)
        tmp = next_tmp
    os.rename(tmp, out_img)
    
def update_mjd_header(evts, stacked_evt):
    # 1) Gather all MJDs from the component evt files
    mjds = []
    for fn in evts:
        try:
            mjds.append(get_date_obs_mjd(fn))
        except Exception as e:
            print(f"    Warning: could not read MJD-OBS from {fn}: {e}")
    if not mjds:
        return

    mjds.sort()
    mjd_beg, mjd_end = mjds[0], mjds[-1]

    # 2) Open the stacked EVENTS file and write new keywords
    #    We assume the events data are in extension 1
    with fits.open(stacked_evt, mode='update') as hdul:
        hdr = hdul[1].header
        hdr['MJD-BEG'] = (mjd_beg, 'Start MJD of first epoch')
        hdr['MJD-END'] = (mjd_end, 'End   MJD of last  epoch')
        
        hdr2 = hdul[0].header
        hdr2['MJD-BEG'] = (mjd_beg, 'Start MJD of first epoch')
        hdr2['MJD-END'] = (mjd_end, 'End   MJD of last  epoch')
        hdul.flush()
    print(f"Wrote MJD-BEG={mjd_beg:.5f}, MJD-END={mjd_end:.5f} to {stacked_evt}")


def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__)
    )
    p.add_argument('--epochs', required=True, help='File listing epochs')
    p.add_argument('--bin-size', type=float, default=1.0,
                   help='Bin size in days')
    p.add_argument('--explosion-mjd', type=float, required=True,
                   help='MJD of explosion')
    p.add_argument('--src-ra', required=True, help='RA of source')
    p.add_argument('--src-dec', required=True, help='Dec of source')
    p.add_argument('--pn-tra', default='POINT')
    p.add_argument('--pn-dec', default='POINT')
    args = p.parse_args()

    # Run pipeline and collect epoch, evt, img
    files = []
    with open(args.epochs) as f:
        for line in f:
            epoch = line.strip()
            if not epoch or epoch.startswith('#'):
                continue
            outdir = run_pipeline(epoch, args.src_ra, args.src_dec,
                                  args.pn_tra, args.pn_dec)
            evt = os.path.join(outdir, f"sw{epoch}xpcw3po_cl.evt")
            img = os.path.join(outdir, f"sw{epoch}xpcw3po_ex.img")
            # only record if both files exist
            if os.path.exists(evt) and os.path.exists(img):
                files.append({'epoch': epoch, 'evt': evt, 'img': img})
            else:
                print(f"Skipping epoch {epoch}: missing evt or img")

    # Bin by days since explosion (MJD)
    bins = {}
    for item in files:
        try:
            mjd = get_date_obs_mjd(item['evt'])
        except Exception as e:
            print(f"  Skipping {item['epoch']}: cannot read DATE-OBS ({e})")
            continue
        dt = mjd - args.explosion_mjd
        idx = int(dt // args.bin_size)
        bins.setdefault(idx, []).append(item)

    # Stack each bin into its own subdir under ./stacks/
    for idx, items in sorted(bins.items()):
        if not items:
            print(f"Bin {idx} empty, skipping.")
            continue
        epochs = sorted(it['epoch'] for it in items)
        first_ep, last_ep = epochs[0], epochs[-1]

        subdir = os.path.join('stacks', f"{first_ep}_to_{last_ep}")
        os.makedirs(subdir, exist_ok=True)

        # Prepare lists of existing files
        evts = [it['evt'] for it in items if os.path.exists(it['evt'])]
        imgs = [it['img'] for it in items if os.path.exists(it['img'])]
        if not evts or not imgs:
            print(f"Bin {first_ep}_to_{last_ep} has no valid files, skipping.")
            continue

        evt_out = os.path.join(subdir, f"evt_{first_ep}to{last_ep}.evt")
        img_out = os.path.join(subdir, f"img_{first_ep}to{last_ep}.img")

        stack_events(evts, evt_out)
        update_mjd_header(evts, evt_out)
        stack_images(imgs, img_out)
        
def cli() -> int | None:       
    return main()      

if __name__ == '__main__':
    main()