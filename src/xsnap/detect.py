from __future__ import annotations
import os
import subprocess
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import xspec
import matplotlib.pyplot as plt
from pathlib import Path


def detect_snr3(evt_file, exp_img_file=None, pha_file=None,
                src_ra=None, src_dec=None, 
                snr_thresh=3.0, match_radius=25.0, details = False):
    """
    Run XIMAGE SNR-based detection and check for a source near given coordinates.

    Args:
        evt_file (str): Path to the X-ray event file (FITS).
        exp_img_file (str, optional): Path to the exposure map image file (FITS). 
            If None, uses the same file as `evt_file`. Defaults to None.
        pha_file (str, optional): Path to a PHA/PI spectrum file. If provided,
            RA/Dec are read from its header. Defaults to None.
        src_ra (float, optional): Source right ascension in degrees. Used if
            `pha_file` is not provided. Defaults to None.
        src_dec (float, optional): Source declination in degrees. Used if
            `pha_file` is not provided. Defaults to None.
        snr_thresh (float, optional): SNR threshold for XIMAGE’s detect command.
            Defaults to 3.0.
        match_radius (float, optional): Matching radius around the source
            position, in arcseconds. Defaults to 25.0.
        details (bool, optional): If True, print matching source details.
            Defaults to False.

    Returns:
        tuple:
            result (dict): Mapping `{abs_evt_file_path: 'detected' | 'not detected'}`.
            df (pandas.DataFrame or None): DataFrame of the last XIMAGE detect
            table with columns
            `['index','rate','rate_err','x_img','y_img','exp',
               'ra_deg','dec_deg','err','hbox']`,
            or None if no sources were found.

    Raises:
        RuntimeError: If the HEADAS environment is not set or files do not exist.
        KeyError: If RA/Dec cannot be read from the PHA header.
        ValueError: If neither `pha_file` nor both `src_ra`/`src_dec` are provided.
    """
    # 1. Ensure HEADAS is set
    if 'HEADAS' not in os.environ:
        raise RuntimeError("$HEADAS environment variable not set. Please source the HEADAS setup.")

    # 2. Build and run the XIMAGE detection script
    bash_script = f"""
source $HEADAS/headas-init.sh
ximage << EOF
read/ecol=pi/emin=30/emax=1000/size=600 {evt_file}
read/exposure/size=600 {exp_img_file}
detect/snr_thresh={snr_thresh}
exit
EOF
"""
    proc = subprocess.run(["bash", "-lc", bash_script],
                          capture_output=True, text=True)
    out = proc.stdout

    # 3. Parse the last detect table
    lines = out.splitlines()
    table_started = False
    data = []
    for line in lines:
        parts = line.split()
        # detect start of data lines by integer index at start
        if len(parts) > 0 and parts[0].isdigit():
            table_started = True
        if table_started and len(parts) >= 12 and parts[0].isdigit():
            # parse columns
            idx = int(parts[0])
            # rate and error are in "val+/-err" format
            rate_str = parts[1]
            if '+/-' in rate_str:
                rate_val, rate_err = rate_str.split('+/-')
            else:
                rate_val, rate_err = rate_str, ''
            rate = float(rate_val)
            rate_err = float(rate_err)
            x_img = float(parts[2])
            y_img = float(parts[3])
            exp = float(parts[4])
            # RA h m s in cols 5,6,7
            ra_h, ra_m, ra_s = map(float, parts[5:8])
            # Dec d m s in cols 8,9,10
            dec_d, dec_m, dec_s = map(float, parts[8:11])
            # error and hbox
            err = float(parts[11])
            hbox = float(parts[12]) if len(parts) > 12 else None
            # convert to degrees
            ra_deg = (ra_h + ra_m/60 + ra_s/3600) * 15.0
            sign = 1 if dec_d >= 0 else -1
            dec_deg = sign * (abs(dec_d) + dec_m/60 + dec_s/3600)

            data.append({
                'index': idx,
                'rate': rate,
                'rate_err': rate_err,
                'x_img': x_img,
                'y_img': y_img,
                'exp': exp,
                'ra_deg': ra_deg,
                'dec_deg': dec_deg,
                'err': err,
                'hbox': hbox,
            })
    if not data:
        df = None
    else:
        df = pd.DataFrame(data)

    # 4. Determine source coordinate to match
    if pha_file and (src_ra is None or src_dec is None):
        hdr = fits.getheader(pha_file, ext=1)
        if 'RA_OBJ' not in hdr and 'RA_TARG' not in hdr:
            hdr = fits.getheader(pha_file, ext=0)
        ra_hdr = hdr.get('RA_OBJ') or hdr.get('RA_TARG')
        dec_hdr = hdr.get('DEC_OBJ') or hdr.get('DEC_TARG')
        if ra_hdr is None or dec_hdr is None:
            raise KeyError(f"No RA/Dec in header of {pha_file}")
        src_coord = SkyCoord(ra_hdr, dec_hdr, unit=u.deg)
    elif src_ra is not None and src_dec is not None:
        src_coord = SkyCoord(src_ra, src_dec, unit=u.deg)
    else:
        raise ValueError("Either pha_file or src_ra/src_dec must be provided.")

    # 5. Check for any detection within match_radius"
    if df is None or df.empty:
        key = os.path.abspath(evt_file)
        return {key: 'not detected'}, df
    cat_coords = SkyCoord(df['ra_deg'].values, df['dec_deg'].values, unit=u.deg)
    sep = cat_coords.separation(src_coord)
    detected = bool((sep <= match_radius * u.arcsec).any())
    if detected and details:
        print(f"Source coordinates: RA={src_coord.ra.deg:.6f} deg, Dec={src_coord.dec.deg:.6f} deg")
        for i, sep_i in enumerate(sep.arcsec):
            if sep_i <= match_radius:
                print(f"Detected at index {df.iloc[i]['index']} with separation {sep_i:.2f} arcsec")
        
    key = os.path.abspath(evt_file)
    result = {key: 'detected' if detected else 'not detected'}
    return result, df

class SourceDetection:
    """
    A class to detect sources (via `XIMAGE <https://heasarc.gsfc.nasa.gov/docs/xanadu/ximage/ximage.html>`_) with certain SNR (default SNR :math:`\\geq` 3) across multiple observation datasets.

    Attributes
    -----------
        pha_paths : array
            Array of absoulte PHA (spectrum) file paths.
        evt_paths : array
            Array of absoulte EVT (event) file paths.
        exp_paths : array
            Array of absoulte EXP (exposure image) file paths.
        results : dict[str, str]
            Source detection results per event files, mapping ``{evt_file: 'detected' | 'not detected'}``.
        detect_tables : dict[str, pandas.DataFrame]
            Table of detected sources in an event files, mapping ``{evt_file: pandas.DataFrame}``.
        tExplosion : float 
            Supernova time of explosion in MJD
    """

    def __init__(self, evt_paths=None, exp_paths=None,  pha_paths=None):
        """
        Initialization of the :py:class:`~xsnap.detect.SourceDetection` class.

        Parameters
        ----------
        evt_paths : array_like, optional
            Paths to event (EVT) files. Defaults to ``None``.
        exp_paths : array_like, optional
            Paths to exposure image (EXP) files. Defaults to ``None``.
        pha_paths : array_like, optional
            Paths to spectrum (PHA) files. Defaults to ``None``.
        """
        self.pha_paths = []   # list of PHA files (abs paths)
        self.evt_paths = []   # matching list of EVT files
        self.exp_paths = []   # matching list of exposure images
        self.results = {}     # {pha: 'detected' | 'not detected'}
        self.detect_tables = {}
        self.tExplosion = None
        if evt_paths:
            self.load(evt_paths, exp_paths, pha_paths)

    def clear(self):
        """
        Reset all stored paths, results, and tables.
        """
        self.pha_paths = []
        self.evt_paths = [] 
        self.exp_paths = []   
        self.results = {}   
        self.detect_tables = {}
        self.tExplosion = None

    def load(self, evt_paths, exp_paths = None, pha_paths = None):
        """
        Load parallel lists of EVT, EXP, and PHA files.

        Parameters
        ----------
            evt_paths : array_like
                Event (EVT) file paths.
            exp_paths : array_like, optional
                Exposure image (EXP) file paths. Defaults to ``None``.
            pha_paths : array_like, optional
                Spectrum (PHA) file paths. Defaults to ``None``.

        Raises
        -------
            ValueError: If provided lists are not all the same length.
            FileNotFoundError: If any path does not exist.
        """
        if (pha_paths is not None) and (exp_paths is not None) and not (len(pha_paths) == len(evt_paths) == len(exp_paths)):
            raise ValueError("pha_paths, evt_paths, exp_paths must be the same length")

        self.evt_paths = [self._abs_exists(p, "EVT") for p in evt_paths]
        
        if exp_paths is not None:
            self.exp_paths = [self._abs_exists(p, "EXP") for p in exp_paths]
        else:
            self.exp_paths = [self._abs_exists(p, "EXP") for p in evt_paths]
        
        if pha_paths is not None:
            self.pha_paths = [self._abs_exists(p, "PHA") for p in pha_paths]

    def __iadd__(self, files):
        """
        Append a dataset via ``+=`` with an array of length 1–3.

        Parameters
        ----------
        files : array_like
            Files in an array_like with options: ``(evt)``; ``(pha, evt)``; or ``(pha, evt, exp)``

        Returns
        -------
            *self* : SourceDetection
                The same :py:class:`~xsnap.detect.SourceDetection` instance.

        Raises
        ------
        RuntimeError
            If array length is not 1, 2, or 3.
        FileNotFoundError
            If any file does not exist.
        """
        if len(files) > 3:
            raise RuntimeError("Use analyzer += (pha, evt, exp) or (pha, evt) or just evt")
        if len(files) == 3:
            pha, evt, exp = files
        elif len(files) == 2:
            pha, evt = files
            exp = evt
        else:
            evt = files
            exp = evt
            pha = None
            
        if pha is not None:
            self.pha_paths.append(self._abs_exists(pha, "PHA"))
        self.evt_paths.append(self._abs_exists(evt, "EVT"))
        self.exp_paths.append(self._abs_exists(exp, "EXP"))
        
        return self

    @staticmethod
    def _abs_exists(path, kind):
        """
        Convert a path to an absolute path and verify its existence.

        Parameters
        ----------
            path (str): The file path to check.
            kind (str): A label for the file type (e.g., "EVT", "PHA", "EXP").

        Returns
        --------
            str: The absolute path.

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        if not Path(path).exists():
            raise FileNotFoundError(f"{kind} file '{path}' does not exist")
        return os.path.abspath(path)

    def detect_all(self, src_ra=None, src_dec=None, snr_thresh=3.0, match_radius=25.0, details = False):
        """
        Run detection with certain SNR (default SNR :math:`\\geq` 3) on all loaded datasets.

        Parameters
        ----------
            src_ra : float, optional
                Source RA in degrees. Must be parsed if there are no PHA files. Defaults to ``None``.
            src_dec : float, optional
                Source Dec in degrees. Must be parsed if there are no PHA files. Defaults to ``None``.
            snr_thresh : float, optional
                SNR threshold for detection. Defaults to ``3.0``. 
            match_radius : float, optional
                Detection matching radius in arcseconds. Defaults to ``25.0``.
            details : bool, optional
                If True, print detailed matches and tables. Defaults to ``False``.

        Returns
        --------
            *None* - Populates :py:attr:`results` and :py:attr:`detect_tables`.

        Raises
        --------
            RuntimeError: If the HEADAS environment is not set or files do not exist.
            KeyError: If RA/Dec cannot be read from the PHA header.
            ValueError: If neither ``pha_file`` nor both ``src_ra``/``src_dec`` are provided.
        """
        pha_empty = (len(self.pha_paths) == 0)
        
        for i, evt in enumerate(self.evt_paths):
            exp = self.exp_paths[i]
            if src_ra is None and src_dec is None and not pha_empty:
                result, df = detect_snr3(evt, exp, pha_file=self.pha_paths[i],
                                     snr_thresh=snr_thresh,
                                     match_radius=match_radius, details=details)
            else:
                try:
                    result, df = detect_snr3(evt, exp, src_ra=src_ra,
                                     src_dec=src_dec, snr_thresh=snr_thresh,
                                     match_radius=match_radius, details=details)
                except Exception as e:
                    raise Exception(e)
                
            # print human‑friendly summary
            status = list(result.values())[0]
            print(f"{Path(evt).name}: {status}")
            if status == "detected" and details:
                print(df)
            self.results[evt] = status
            self.detect_tables[evt] = df
            
    def show_source(self, pha_files=None, cmap='viridis'):
        """
        Display the source region image through the pha spectrum files.

        Parameters
        ----------
            pha_files : array_like, optional
                Paths to pha spectrum files to display. If provided, replaces ``self.pha_paths``. Defaults to ``None``.
            cmap : str, optional
                Matplotlib colormap name for ``plt.imshow()``. Defaults to ``viridis``.

        Raises
        --------
            RuntimeError: If neither ``pha_files`` nor :py:attr:`pha_paths` contains any file paths.

        Returns
        --------
            *None*
        """
        if len(self.pha_paths) != 0 or pha_files is not None:
            if pha_files is not None:
                self.pha_paths = pha_files
                
            for pha in self.pha_paths:
                hdul = fits.open(pha)
                data = hdul[0].data
                im = plt.imshow(data, cmap=cmap, vmin=data.min(), vmax=data.max()) 
                plt.colorbar()
                plt.title(f"{pha}")
                plt.show()
                
            return None
            
        raise RuntimeError("Please input the pha paths!")