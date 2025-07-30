from __future__ import annotations
import xspec
import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
import pandas as pd
import sys, os, subprocess
from pathlib import Path
from shutil import which
from astropy.io import fits
from astropy.time import Time
from .detect import SourceDetection

# Constants
C_LIGHT_KM_S = 299792.458       # speed of light in km/s
MPC_TO_CM = 3.085677581e24       # 1 Mpc in cm
    

class SpectrumFit:
    """
    A class to handle a single `XSPEC <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/>`_ spectral fit.

    Attributes
    ----------
    fluxes : dict[str, pandas.DataFrame] or None
        Stored absorbed and unabsorbed flux DataFrames, keyed "absorbed" and "unabsorbed", 
        with columns are ``['data', 'model', 'flux', 'lo_flux_err', 'hi_flux_err', 'phot', 'lo_phot_err', 'hi_phot_err']``
    lumin : pandas.DataFrame or None
        Stored luminosity DataFrame with columns ``['data','model','lumin','lo_lumin_err','hi_lumin_err']``.
    params : pandas.DataFrame or None
        Parameter values and uncertainties with columns ``['data','model', '<component>_<param>', 'lo_<component>_<param>_err', 'hi_<component>_<param>_err']``.
    counts : pandas.DataFrame or None
        Count-rate DataFrame with columns ``['data','model','net_rate','net_err','total_rate','model_rate']``.
    models : list[str]
        Model expressions applied.
    obstime : pandas.DataFrame or None
        Observation time DataFrame with columns ``['data','obs_time','obs_time_err']``
        plus ``['time_since_explosion','time_since_explosion_err']`` if ``tExplosion`` 
        given when calculating ``obstime``.
    tExplosion : float or None
        Supernova time of explosion in MJD.
    pha : list[str]
        Loaded PHA file paths.
    detection : :py:class:`~xsnap.detect.SourceDetection` or None
        Associated :py:class:`~xsnap.detect.SourceDetection` object.
    """
    def __init__(self, abund: str = "aspl", seed: int = None, mute: bool = False):
        """Initialize of the :py:class:`~xsnap.spectrum.SpectrumFit` class

        Parameters
        ----------
        abund : str, optional 
            Abundance table for XSPEC. Defaults to ``"aspl"``.
        seed : int, optional
            Random seed for reproducibility. Defaults to ``None``.
        mute : bool, optional
            Mute XSPEC output, i.e. set ``xspec.Xset.chatter = 0``. Defaults to ``False``

        Raises
        -------
            RuntimeError
                If the ``xspec`` executable is not found in ``$PATH``.
        """
        
        if not which("xspec"):
            raise RuntimeError("Xspec is not available in the PATH. Please install it.")
        
        # Initialize Xspec
        if seed is not None:
            xspec.Xset.seed = seed
        if mute:
            xspec.Xset.chatter = 0
        xspec.Xset.abund = abund
        
        self.fluxes = None
        self.lumin = None
        self.params = None
        self.counts = None
        self.models = []
        self.obstime = None
        self.tExplosion = None
        self.pha = []
        self.detection = None
        
    def __iadd__(self, pha):
        """
        Allow ``spec += pha`` syntax to load a new PHA file.

        Parameters
        ----------
            pha : str
                Path to the PHA file to load. The PHA file must be a grouped PHA 
                file or the background and response file must have the same name 
                as the PHA file for PyXspec to automatically get the background and response file.

        Returns
        --------
            *self* : SpectrumFit
                The same :py:class:`~xsnap.spectrum.SpectrumFit` instance, but with the pha file loaded.
        """
        self.load_data(pha)
        return self
    
    def clear(self):
        """
        Clear all loaded data, models, and derived results.
        """
        xspec.AllData.clear()
        xspec.AllModels.clear()
        self.pha = []
        self.models = []
        self.fluxes = {}
        self.obstime = None
        self.tExplosion = None
        self.params = None
        self.lumin = None
        self.counts = None
        self.detection = None
        print("Cleared all data and models.")
    
    def load_data(self, pha, newGroup: bool = False, clear: bool = False, 
                  bad: bool = True, detection: SourceDetection | None = None):
        """
        Load a grouped PHA file into XSPEC.

        Parameters
        ----------
            pha : str
                Path to the PHA file. Path to the PHA file to load. The PHA file must be a grouped PHA 
                file or the background and response file must have the same name 
                as the PHA file for PyXspec to automatically get the background and response file.
            newGroup : bool, optional
                Load into a new group if ``True``. Defaults to ``False``.
            clear : bool, optional
                Clear existing data before loading. Defaults to ``False``.
            bad : bool, optional
                Ignore bad channels if ``True``. Defaults to ``True``.
            detection : SourceDetection, optional
                Detection object to attach. Defaults to ``None``.

        Raises
        --------
            ValueError
                If ``pha`` is None.
            FileNotFoundError
                If the PHA file or its directory does not exist.
        """
        if pha is None:
            raise ValueError("Spectrum file must be specified.")
        
        pha_abspath = os.path.abspath(pha)
        pha_dir = os.path.dirname(pha_abspath)
        if not os.path.exists(pha_abspath):
            raise FileNotFoundError(f"Spectrum file {pha} does not exist.")
        if not os.path.exists(pha_dir):
            raise FileNotFoundError(f"Directory {pha_dir} does not exist.")
        
        if clear:
            xspec.AllData.clear()
            
        Ngrp = xspec.AllData.nGroups
        Nspec = xspec.AllData.nSpectra
        
        # Load the PHA file into Xspec
        original_dir = os.getcwd()
        os.chdir(pha_dir)
        
        if newGroup:
            # Load into a new group
            xspec.AllData(f"{Ngrp+1}:{Nspec+1} {pha_abspath}")
        else:
            # Load into the current group
            xspec.AllData += pha_abspath
        
        if bad:
            # Ignore bad channels
            xspec.AllData.ignore("bad")
            
        self.pha.append(pha_abspath)
        self.detection = detection
        
        os.chdir(original_dir)            
        
    def ignore(self, ignore=None):
        """
        Ignore specified data channels in XSPEC.

        Parameters
        ----------
            ignore : str, optional
                same ignore command as in XSPEC. Defaults to ``None``.
        """
        if ignore is not None:
            xspec.AllData.ignore(ignore)
            
    def show(self):
        """
        Display the currently loaded data and models in XSPEC format.
        Wouldn't work if XSPEC is muted, i.e. ``xspec.Xset.chatter = 0``.
        """
        if xspec.AllData.nGroups == 0:
            print("No data loaded.")
            return
        xspec.AllData.show()
        if xspec.AllModels.nModels > 0:
            xspec.AllModels.show()
        else:
            print("No models defined.")
            
    def set_rebin(self, minSig=None, maxBins=None, groupNum=None, errType=None):
        """
        Set rebinning parameters for XSPEC plots.

        Parameters
        ----------
            minSig : float, optional
                Minimum significance per bin. Defaults to ``None``.
            maxBins : int, optional
                Maximum number of bins. Defaults to ``None``.
            groupNum : int, optional
                Specific plot group number to rebin. Defaults to ``None``.
            errType : str, optional
                Error type for rebinning. Valid entries are ``"quad", "sqrt", "poiss-1","poiss-2", "poiss-3"``. Defaults to ``None``.
        """
        
        xspec.Plot.setRebin(minSig=minSig, maxBins=maxBins, groupNum=groupNum, errType=errType)
        
    def set_plot(self, args: str ="data", device: str ="/null", xAxis: str ="keV", fileName: str ="plot.png"):
        """
        Configure and execute an XSPEC plot.

        Parameters
        ----------
            args : str, optional
                Plot command string (e.g., ``"data", "ldata"``), same command as in XSPEC.
                Defaults to ``"data"``.
            device : str, optional
                XSPEC plot device. If ``"/null"``, plot via matplotlib.
                Defaults to ``"/null"``.
            xAxis : str, optional
                Units for the x-axis. Valid options: ``"channel", "keV", "MeV", "GeV", "Hz", "angstrom", "cm", "micron", "nm"`` (case-insensitive).
                Defaults to ``"keV"``.
            fileName : str, optional
                Base filename for saving the plot PNG.
                Defaults to ``"plot"``.
                
        Raise
        ------
            Exception
                When invalid xAxis or device string is parsed.

        Returns
        --------
            Plot : tuple(matplotlib.figure.Figure, matplotlib.axes.Axes) or None
                If ``device == "/null"``, returns the Figure and Axes so the user
                can customize further. Otherwise, returns None.
        """
        # configure units
        try:
            xspec.Plot.xAxis = xAxis
            print(f"X-axis set to '{xAxis}'")
        except Exception as e:
            print(f"Invalid xAxis '{xAxis}': {e}")
            print('Valid options: "channel", "keV", "MeV", "GeV", "Hz", '
                  '"angstrom", "cm", "micron", "nm" (case-insensitive)')
        # configure device
        try:
            xspec.Plot.device = device
            print(f"Plot device set to '{device}'")
        except Exception as e:
            print(f"Invalid device '{device}': {e}")

        xspec.Plot(args)

        # if device /null
        if xspec.Plot.device == "/null":
            print("Plotting to null device → matplotlib Figure returned.")

            x    = xspec.Plot.x()
            y    = xspec.Plot.y()
            xErr = xspec.Plot.xErr()
            yErr = xspec.Plot.yErr()

            fig, ax = plt.subplots(figsize=(8,6))
            ax.errorbar(x, y, xerr=xErr, yerr=yErr, fmt='o', label='Data')
            if len(xspec.AllModels.sources) > 0:
                y_model = xspec.Plot.model()
                ax.plot(x, y_model, label='Model')

            xl, yl, title = xspec.Plot.labels()
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            ax.set_title(title)
            if "ldata" in args:
                ax.loglog()
            ax.legend()
            ax.grid()

            fig.savefig(f"{fileName}")
            plt.show()
            return fig, ax 

        return None 
       
    def set_model(self, model_string: str, clear: bool = True, data: int = 0, **kwargs):
        """
        Build and configure the XSPEC model for all data groups.

        Parameters
        ----------
            model_string : str
                XSPEC model expression (e.g., "tbabs*pow").
            clear : bool, optional
                If ``True``, clear existing models. Defaults to True.
            data : int, optional
                Index of PHA header to use for RA/Dec in computing nH.
                Defaults to ``0``.
            **kwargs
                Component parameters in ``component_Param`` format, e.g.,
                ``TBabs_nH="0.059 -1"`` or ``powerlaw_PhoIndex=2``.

        Raises
        -------
            KeyError
                If RA/Dec cannot be found in the PHA header when computing nH.
            RuntimeError
                If HEADAS or nH tool invocation fails.
            ValueError
                If a kwarg key is not "Component_Param".
            AttributeError
                If the specified component does not exist in the model.
                
        Example usage
            1. An absorbed power-law model
        
            .. code-block:: python
            
                SpectrumFit.set_model(
                    "tbabs*ztbabs*pow",
                    Tbabs_nH="0.059 -1",         # component TBabs, parameter nH
                    zTBabs_nH=0.5,               # component zTBabs, parameter nH
                    zTBabs_Redshift=0,           # component zTBabs, parameter Redshift
                    powerlaw_PhoIndex=2,         # component powerlaw, parameter PhoIndex
                    powerlaw_norm=1              # component powerlaw, parameter norm
                )
                
            2. An absorbed Thermal-bremsstrahlung model (with ``TBabs.nH`` calculated automatically from `HEASoft <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_)
            
            .. code-block:: python

                SpectrumFit.set_model(
                    "tbabs*ztbabs*pow",
                    data=1,                      # which pha data to get the RA/Dec from. To get the TBabs nH parameter, default is 0.
                    zTBabs_nH=0.5,               # component zTBabs, parameter nH
                    zTBabs_Redshift=0,           # component zTBabs, parameter Redshift
                    powerlaw_PhoIndex=2,         # component powerlaw, parameter PhoIndex
                    powerlaw_norm=1              # component powerlaw, parameter norm
                )
        """
        # If TBabs_nH not provided, compute from header coordinates
        if 'TBabs_nH' not in kwargs and 'tbabs' in model_string.lower():
            pha_file = self.pha[data]
                        # read header: try ext=1, fallback to ext=0
            hdr = fits.getheader(pha_file, ext=1)
            if 'RA_OBJ' not in hdr and 'RA_TARG' not in hdr:
                hdr = fits.getheader(pha_file, ext=0)
            ra = hdr.get('RA_OBJ') or hdr.get('RA_TARG')
            dec = hdr.get('DEC_OBJ') or hdr.get('DEC_TARG')
            if ra is None or dec is None:
                raise KeyError(f"No RA/Dec in header of {pha_file}")
                        # initialize HEASoft environment
            subprocess.run('source $HEADAS/headas-init.sh', shell=True, executable='/bin/bash', check=True)
            # prepare input for nH tool
            nh_input = f"2000\n{ra}\n{dec}\n"            # run nH, feeding input via pipe
            proc = subprocess.run(
                f'printf "{nh_input}" | nH',
                shell=True, executable='/bin/bash',
                capture_output=True, text=True, check=True
            )
            out = proc.stdout
            # parse weighted average nH (cm^-2) (cm^-2)
            for line in out.splitlines():
                if 'Weighted average nH' in line:
                    val = float(line.split()[-1])
                    break
            else:
                raise RuntimeError("Failed to parse nH output")
            # convert to 1e22 units and freeze
            tbabs_val = val / 1e22
            kwargs['TBabs_nH'] = f"{tbabs_val:.5f} -1"
        # clear existing models
        if clear:
            xspec.AllModels.clear()
        # attach base model to group1
        nSources = len(xspec.AllModels.sources)
        xspec.Model(model_string, sourceNum=nSources+1)
        # clone into further groups
        for _ in range(2, xspec.AllData.nGroups + 1):
            xspec.AllModels += model_string
        # set parameters
        for grp in range(1, xspec.AllData.nGroups + 1):
            mod = xspec.AllModels(grp)
            for key, val in kwargs.items():
                try:
                    comp, par = key.split('_', 1)
                except ValueError:
                    raise ValueError(f"Kwarg '{key}' must be 'Component_Param'")
                if not hasattr(mod, comp):
                    raise AttributeError(f"Model has no component '{comp}'")
                setattr(getattr(mod, comp), par, val)
        # show all models
        xspec.AllModels.show()
        # store model for reference
        self.models.append(model_string)
        
    def simulate(self, nIterations: int =1000, statMethod: str ="cstat",
                 low_energy: float = 0.3, high_energy: float = 10.0,
                 unabsorbed: bool =True, plot=True, kwargs_plot="data", device_plot="/svg"):
        """
        Simulate and count flux upper-limits using the ``FakeIt`` command from XSPEC.
        Will populate :py:attr:`fluxes` and :py:attr:`counts`
        
        Parameters
        ----------
            nIterations : int, optional
                Number of fit iterations. Defaults to ``1000``.
            statMethod : str, optional
                Type of fit statistic method. Valid names: ``'chi', 'cstat', 'lstat', 'pgstat', 'pstat', 'whittle'``. Defaults to ``"cstat"``.
            low_energy : float, optional
                Lower energy bound in keV. Defaults to ``0.3``.
            high_energy : float, optional
                Upper energy bound in keV. Defaults to ``10.0``.
            unabsorbed : bool, optional
                If ``True``, compute unabsorbed flux as well. Defaults to ``True``.
            plot : bool, optional
                If ``True``, plot the simulated data. Defaults to ``True``.
            kwargs_plot : str, optional
                Plot command string (e.g., ``"data", "ldata"``). Defaults to ``"data"``.
            device_plot : str, optional
                XSPEC plot device. If ``"/null"``, plot via matplotlib. Defaults to ``"/svg"``.

        Raises
        ------
            RuntimeError
                When there are no model set, 
                When there are no spectrum loaded,
                When cannot compute unabsorbed flux due to no tbabs/ztbabs attribute in the model.

        Returns
        --------
            Flux upper-limits : dict[str, pandas.DataFrame]
                The same absorbed and unabsorbed flux DataFrames as the one in :py:attr`fluxes`
                with columns ``['data', 'model', 'flux', 'lo_flux_err', 'hi_flux_err', 'phot', 'lo_phot_err', 'hi_phot_err']``.
        """
        original_dir = os.getcwd()
        if len(self.models) == 0:
            raise RuntimeError("Please set model first!")
        
        if xspec.AllData.nSpectra == 0:
            raise RuntimeError("Please load a spectrum file first!")
                
        model_str = self.models[len(self.models) - 1]
        
        rows_abs = []
        
        rows_ct = []
        
        rows_unabs = []
        
        for i, pha in enumerate(self.pha):
            
            data_name = os.path.basename(pha)
            pha_dir = os.path.dirname(pha)
            ct_rate = xspec.AllData(i+1).rate[2]
            
            rows_ct.append({
                'data': data_name,
                'model': self.models[len(self.models)-1] if self.models else '',
                'net_rate': xspec.AllData(i+1).rate[0],
                'net_err': xspec.AllData(i+1).rate[1],
                'total_rate': xspec.AllData(i+1).rate[2],
                'model_rate': xspec.AllData(i+1).rate[3]
            })
            
            os.chdir(pha_dir)
            
            fs = xspec.FakeitSettings()
            xspec.AllData.fakeit(1, fs) 
            
            xspec.AllData.ignore(f"**-{low_energy} {high_energy}-**")
            xspec.Fit.statMethod = statMethod
            xspec.Fit.nIterations = nIterations
            xspec.Fit.perform()
            
            xspec.AllModels.calcFlux(f"{low_energy} {high_energy}")
            
            flux_to_ct = xspec.AllData(1).flux[0]/xspec.AllData(1).rate[0]
            ul_flux_abs = flux_to_ct * ct_rate
            
            phot_to_ct = xspec.AllData(1).flux[3]/xspec.AllData(1).rate[0]
            ul_phot_abs = phot_to_ct * ct_rate
            
            if plot:
                self.set_plot(kwargs_plot, device_plot)
            
            if unabsorbed:
            
                first = xspec.AllModels(1)
                if not (hasattr(first,'TBabs') or hasattr(first,'zTBabs')):
                    raise RuntimeError("Cannot compute unabsorbed flux: missing TBabs/zTBabs.")
                
                m = xspec.AllModels(1)
                # Zero Galactic absorber if it’s there
                if hasattr(m, "TBabs"):
                    m.TBabs.nH = 0
                    
                # Zero red-shifted absorber if it’s there
                if hasattr(m, "zTBabs"):
                    m.zTBabs.nH = 0
                    
                xspec.AllModels.calcFlux(f"{low_energy} {high_energy}")
                
                flux_to_ct_unabs = xspec.AllData(1).flux[0]/xspec.AllData(1).rate[0]
                ul_flux_unabs = flux_to_ct_unabs * ct_rate
                
                phot_to_ct_unabs = xspec.AllData(1).flux[3]/xspec.AllData(1).rate[0]
                ul_phot_unabs = phot_to_ct_unabs * ct_rate
            
            rows_abs.append({
                'data': data_name,
                'model': model_str,
                'flux': ul_flux_abs,
                'lo_flux_err': 0,
                'hi_flux_err': 0,
                'phot': ul_phot_abs,
                'lo_phot_err': 0,
                'hi_phot_err': 0
            })
            
            if unabsorbed:
                rows_unabs.append({
                    'data': data_name,
                    'model': model_str,
                    'flux': ul_flux_unabs,
                    'lo_flux_err': 0,
                    'hi_flux_err': 0,
                    'phot': ul_phot_unabs,
                    'lo_phot_err': 0,
                    'hi_phot_err': 0
                })
                
            os.chdir(original_dir)
            
        df_ul_abs = pd.DataFrame(rows_abs)
        df_ul_unabs = pd.DataFrame(rows_unabs)
        
        result = {'absorbed': df_ul_abs, 'unabsorbed': df_ul_unabs}
        
        self.fluxes = result
        self.counts = pd.DataFrame(rows_ct)
        
        return self.fluxes
               
    
    def fit(self, nIterations: int =1000, statMethod: str ="cstat"):
        """
        Fit the model to the loaded data through XSPEC.

        Parameters
        ----------
            nIterations : int, optional
                Number of fit iterations. Defaults to ``1000``.
            statMethod : str, optional
                Type of fit statistic method. Valid names: ``'chi', 'cstat', 'lstat', 'pgstat', 'pstat', 'whittle'``. Defaults to ``"cstat"``.
        """
        xspec.Fit.renorm()
        xspec.Fit.statMethod = statMethod
        xspec.Fit.nIterations = nIterations
        xspec.Fit.perform()
        
    def _getObsTime(self, file):
        """
        Internal: read FITS header to get observation MJD and uncertainty.

        Args:
            file (str): Path to the FITS file.

        Returns:
            tuple:
                mjd_obs (float): Mid-point MJD of observation.
                mjd_err (float): Uncertainty in MJD.

        Raises:
            RuntimeError: If file does not exist.
            KeyError: If no DATE-OBS or MJD-OBS keywords found.
        """
        
        file = os.path.abspath(file)
        if not Path(file).exists():
            raise RuntimeError(f"The file {file} does not exist!")
        
        hdr = fits.getheader(file, ext=1)
        
        if 'MJD-BEG' in hdr and 'MJD-END' in hdr:
            mjd_beg = float(hdr['MJD-BEG'])
            mjd_end = float(hdr['MJD-END'])
            
            mjd_obs = 0.5 * (mjd_beg + mjd_end)
            mjd_err = np.abs(mjd_end - mjd_obs)
        else:
            mjd_err = 0
            
            # If time keywords missing, try primary HDU
            if 'MJD-OBS' not in hdr and 'DATE-OBS' not in hdr:
                hdr = fits.getheader(file, ext=0)

            # --- derive observation MJD
            if 'MJD-OBS' in hdr:
                mjd_obs = float(hdr['MJD-OBS'])
            else:
                date_obs = hdr.get('DATE-OBS')
                if date_obs is None:
                    raise KeyError(f"No DATE-OBS or MJD-OBS in {file}")
                timesys = hdr.get('TIMESYS', 'UTC').lower()
                mjd_obs = Time(date_obs, format='isot', scale=timesys).mjd
        
        return mjd_obs, mjd_err
        
    def get_time(self, tExplosion: float =None):
        """
        Build a DataFrame of observation times and, if provided, time since explosion.

        Parameters
        ----------
            tExplosion : float, optional
                Supernova time of explosion in MJD. Defaults to ``None``.

        Returns
        ---------
            Observation times : pandas.DataFrame
                Observation times wrapped in DataFrame table with columns ``['data','obs_time','obs_time_err']``
                plus ``['time_since_explosion','time_since_explosion_err']`` if ``tExplosion`` given.

        Raises
        --------
            ValueError
                If computed ``time_since_explosion`` is negative.
        """
        cols = ['data', 'obs_time', 'obs_time_err'] + ([] if tExplosion is None else ['time_since_explosion', 'time_since_explosion_err'])

        if not self.pha:          # no files loaded
            return pd.DataFrame(columns=cols)

        records = []
        for pha_file in self.pha:
            mjd_obs, mjd_err = self._getObsTime(pha_file)
            
            # --- build row
            row = {'data': os.path.basename(pha_file), 'obs_time': mjd_obs, 'obs_time_err': mjd_err}
            
            if tExplosion is not None:
                t_days = mjd_obs - self.tExplosion
                
                if t_days < 0:
                    raise ValueError("Observation date is before explosion date.")
                
                t_err = mjd_err 
                row['time_since_explosion'] = t_days
                row['time_since_explosion_err'] = t_err

            records.append(row)

        obstime_df = pd.DataFrame.from_records(records, columns=cols)
        self.obstime = obstime_df
        self.tExplosion = tExplosion
        return obstime_df
    
    def get_counts(self):
        """
        Return a DataFrame of count rates for each loaded spectrum.

        Returns
        --------
            Count rates : pandas.DataFrame
                Count rates wrapped in DataFrame table with columns ``['data','model','net_rate','net_err','total_rate','model_rate']``.

        Raises
        -------
            RuntimeError
                If no spectra are loaded.
        """
        nspec = xspec.AllData.nSpectra
        if nspec == 0:
            raise RuntimeError("No spectra loaded.")

        rows = []
        for i in range(1, nspec+1):
            spec      = xspec.AllData(i)
            data_name = os.path.basename(self.pha[i-1])

            # use the spec.rate tuple
            net_rate, net_err, total_rate, model_rate = spec.rate

            rows.append({
                'data': data_name,
                'model': self.models[len(self.models)-1] if self.models else '',
                'net_rate': net_rate,
                'net_err': net_err,
                'total_rate': total_rate,
                'model_rate': model_rate
            })
        self.counts = pd.DataFrame(rows)
        return self.counts
    
    def get_params(self, error_args, replace=False):
        """
        Compute best-fit parameters and their uncertainties and store them.

        Parameters
        ----------
            error_args : str or list
                The same arguments passed to :func:`xspec.Fit.error`. 
                Details on the arguments can be found `here <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/fitmanager.html>`_.
            replace : bool, optional
                If True, overwrite existing ``params``.  Defaults to ``False``.

        Returns
        -------
            Parameters : pandas.DataFrame
                Best-fit parameters and uncertainties wrapped in a DataFrame table with columns
                ``['data','model', '<component>_<param>', 'lo_<component>_<param>_err', 'hi_<component>_<param>_err']``.
        """
        # compute errors in XSPEC
        xspec.Fit.error(error_args)
        ngrp = xspec.AllData.nGroups
        records = []
        for idx in range(1, ngrp+1):
            row = {
                'data':  os.path.basename(self.pha[idx-1]),
                'model': self.models[len(self.models)-1]
            }
            mod = xspec.AllModels(idx)
            for comp in mod.componentNames:
                comp_obj = getattr(mod, comp)
                for pname in comp_obj.parameterNames:
                    p    = getattr(comp_obj, pname)
                    val  = p.values[0] if p.values else np.nan
                    errs = p.error or []
                    if p.frozen or not errs or len(errs) < 2 or (errs[0] == 0 and errs[1] == 0):
                        lo = hi = 0.0
                    else:
                        lo = val - errs[0]
                        hi = errs[1] - val
                    key  = f"{comp}_{pname}"
                    row[key]           = val
                    row[f"lo_{key}_err"] = lo
                    row[f"hi_{key}_err"] = hi
            records.append(row)

        df = pd.DataFrame(records)
        if self.params is None or replace:
            self.params = df
        else:
            self.params = pd.concat([self.params, df], ignore_index=True)
        return self.params
                
    def get_fluxes(self, low_energy: float = 0.3, high_energy: float = 10.0,
              errMCMC: float =1000, CI: float =68.0, unabsorbed: bool =True):
        """
        Calculate absorbed (and optionally unabsorbed) fluxes for each spectrum.

        Parameters
        ----------
            low_energy : float, optional
                Lower energy bound in keV. Defaults to ``0.3``.
            high_energy : float, optional
                Upper energy bound in keV. Defaults to ``10.0``.
            errMCMC : int, optional
                Number of MCMC error iterations. Defaults to ``1000``.
            CI : float, optional
                Confidence interval percentage. Defaults to ``68.0``.
            unabsorbed : bool, optional
                If True, compute unabsorbed flux as well. Defaults to ``True``.

        Returns
        ---------
            Fluxes : dict[str, pandas.DataFrame]
                The same absorbed and unabsorbed flux DataFrames as the one in :py:attr`fluxes`
                with columns ``['data', 'model', 'flux', 'lo_flux_err', 'hi_flux_err', 'phot', 'lo_phot_err', 'hi_phot_err']``.

        Raises
        -------
            RuntimeError
                If no spectra loaded or missing absorber for unabsorbed flux.
        """
        nspec = xspec.AllData.nSpectra
        ngrp  = xspec.AllData.nGroups
        if nspec == 0:
            raise RuntimeError("No spectra loaded.")
        if len(self.models) < ngrp:
            raise RuntimeError("Need one model entry per DataGroup in self.models.")

        # --- absorbed fluxes ---
        cmd = f"{low_energy} {high_energy} err {errMCMC} {CI}"
        xspec.AllModels.calcFlux(cmd)

        abs_rows, rel_errs = [], []
        for i in range(1, nspec+1):
            spec      = xspec.AllData(i)
            data_name = os.path.basename(self.pha[i-1])
            model_str = self.models[len(self.models)-1]
            flux, lo_f, hi_f, phot, lo_p, hi_p = spec.flux

            lo_e  = max(0.0, flux - lo_f)
            hi_e  = max(0.0, hi_f  - flux)
            plo_e = max(0.0, phot - lo_p)
            phi_e = max(0.0, hi_p  - phot)

            abs_rows.append({
                'data': data_name,
                'model': model_str,
                'flux': flux,
                'lo_flux_err': lo_e,
                'hi_flux_err': hi_e,
                'phot': phot,
                'lo_phot_err': plo_e,
                'hi_phot_err': phi_e
            })
            rel_errs.append((lo_e/flux if flux else 0.0,
                            hi_e/flux if flux else 0.0))

        df_abs = pd.DataFrame(abs_rows)
        result = {'absorbed': df_abs}

        # --- unabsorbed fluxes ---
        if unabsorbed:
            # --- zero every “nH”-type parameter in every DataGroup -------------
            for g in range(1, ngrp + 1):
                mod = xspec.AllModels(g)

                for comp_name in mod.componentNames:    
                    comp = getattr(mod, comp_name)

                    # look through that component’s parameters
                    for pname in comp.parameterNames:         
                        if pname.lower() == "nh":              # matches “nH” / “NH” / “nh”
                            setattr(comp, pname, 0)            # comp.nH = 0

            xspec.AllModels.calcFlux(f"{low_energy} {high_energy}")

            un_rows = []
            for i, (rlo, rhi) in enumerate(rel_errs, start=1):
                spec      = xspec.AllData(i)
                data_name = os.path.basename(self.pha[i-1])
                model_str = self.models[len(self.models)-1]
                flux_u, _, _, phot_u, _, _ = spec.flux

                un_rows.append({
                    'data': data_name,
                    'model': model_str,
                    'flux': flux_u,
                    'lo_flux_err': flux_u * rlo,
                    'hi_flux_err': flux_u * rhi,
                    'phot': phot_u,
                    'lo_phot_err': phot_u * rlo,
                    'hi_phot_err': phot_u * rhi
                })

            df_unabs = pd.DataFrame(un_rows)
            result['unabsorbed'] = df_unabs

        # --- accumulate into self.fluxes ---
        if self.fluxes is None:
            self.fluxes = result
        else:
            # append absorbed
            self.fluxes['absorbed'] = pd.concat(
                [self.fluxes['absorbed'], df_abs], ignore_index=True
            )
            # append unabsorbed
            if 'unabsorbed' in result:
                prev = self.fluxes.get('unabsorbed')
                if prev is None:
                    self.fluxes['unabsorbed'] = result['unabsorbed']
                else:
                    self.fluxes['unabsorbed'] = pd.concat(
                        [prev, result['unabsorbed']], ignore_index=True
                    )
        return self.fluxes
    
    def get_lumin(self, fluxes: pd.DataFrame, model_name: str = None,
                 distance: float = None, redshift: float = None,
                 lo_dist_err: float = None, hi_dist_err: float = None,
                 lo_z_err: float = None, hi_z_err: float = None,
                 H0: float = 70.0, replace: bool = False):
        """
        Compute luminosities from given fluxes.

        Parameters
        ----------
            fluxes : pandas.DataFrame
                Must contain ``['data','model','flux','lo_flux_err','hi_flux_err']``.
            model_name : str, optional
                Override per-row model names. Defaults to ``None``.
            distance : float, optional
                Distance in Mpc. Required if redshift is None. Defaults to ``None``.
            redshift : float, optional
                Redshift. Required if distance is ``None``. Defaults to ``None``.
            lo_dist_err : float, optional
                Lower uncertainty of distance in Mpc. Defaults to ``None``.
            hi_dist_err : float, optional
                Upper uncertainty of distance in Mpc. Defaults to ``None``.
            lo_z_err : float, optional
                Lower uncertainty of redshift. Defaults to ``None``.
            hi_z_err : float, optional
                Upper uncertainty of redshift. Defaults to ``None``.
            H0 : float, optional
                Hubble constant in km/s/Mpc. Defaults to ``70.0``.
            replace : bool, optional
                If True, overwrite existing :py:attr:`lumin`. Defaults to ``False``.

        Returns
        --------
            Luminosity : pandas.DataFrame
                Luminosity wrapped in DataFrame table with columns ``['data','model','lumin','lo_lumin_err','hi_lumin_err']``.

        Raises
        -------
            ValueError
                If neither distance nor redshift is provided.
        
        .. note::
            If ``distance`` is not supplied, it is inferred from ``redshift`` via the Doppler relation:
            
            .. math::
                v = 
                \\begin{cases}
                    c\,z, & z \\ll 1, \\\\
                    c\,\\frac{(1 + z)^2 - 1}{(1 + z)^2 + 1}, & \\text{otherwise.}
                \\end{cases}
                
            and then:

            .. math::
                d = \\frac{v}{H_0} \ (\\mathrm{Mpc}).
            
            References:
                Hogg, D. W. (1999). *Distance measures in cosmology*. `arXiv:astro-ph/9905116 <https://doi.org/10.48550/arXiv.astro-ph/9905116>`_
        """
        # ── distance (Mpc) & its asymmetric errors ───────────────────────────
        if distance is None and redshift is None:
            raise ValueError("Either distance or redshift must be provided.")

        if distance is None:
            z = redshift
            # simple approximation: v = cz   for z <~ 0.1
            if z < 0.1:
                v = z * C_LIGHT_KM_S
            else:                            # special-relativistic formula
                v = C_LIGHT_KM_S * ((1 + z)**2 - 1) / ((1 + z)**2 + 1)
            distance = v / H0               # Mpc

            # propagate z-uncertainties → distance errors if supplied
            if lo_z_err is not None:
                distance_lo = (z - lo_z_err) * C_LIGHT_KM_S / H0 if z - lo_z_err >= 0 else 0
                lo_dist_err = distance - distance_lo
            if hi_z_err is not None:
                distance_hi = (z + hi_z_err) * C_LIGHT_KM_S / H0
                hi_dist_err = distance_hi - distance
        else:
            # distance given explicitly; propagate any supplied distance errors
            lo_dist_err = lo_dist_err or 0.0
            hi_dist_err = hi_dist_err or 0.0

        distance = float(distance)
        lo_dist_err = float(lo_dist_err or 0.0)
        hi_dist_err = float(hi_dist_err or 0.0)

        # convert to cm
        R_cm  = distance * MPC_TO_CM
        fac   = 4.0 * np.pi * R_cm**2

        # pre-compute D_low / D_high for distance error propagation
        D_lo_cm = max(distance - lo_dist_err, 0) * MPC_TO_CM
        D_hi_cm = (distance + hi_dist_err) * MPC_TO_CM
        fac_lo  = 4.0 * np.pi * D_lo_cm**2
        fac_hi  = 4.0 * np.pi * D_hi_cm**2

        # ── build output records row-by-row ──────────────────────────────────
        records = []
        for _, row in fluxes.iterrows():
            f   = row['flux']
            f_lo  = row.get('lo_flux_err', 0.0)
            f_hi  = row.get('hi_flux_err', 0.0)

            L      = f * fac
            # flux-component of the errors
            dL_lo_f = f_lo * fac
            dL_hi_f = f_hi * fac
            # distance-component  (asymmetric)
            dL_lo_d = max(0.0, L - f * fac_lo)
            dL_hi_d = max(0.0, f * fac_hi - L)

            lo_tot = np.hypot(dL_lo_f, dL_lo_d) 
            hi_tot = np.hypot(dL_hi_f, dL_hi_d)

            records.append({
                'data':  row['data'],
                'model': model_name if model_name is not None else row['model'],
                'lumin': L,
                'lo_lumin_err': lo_tot,
                'hi_lumin_err': hi_tot
            })

        df_lum = pd.DataFrame.from_records(
            records,
            columns=['data', 'model', 'lumin', 'lo_lumin_err', 'hi_lumin_err']
        )

        if replace or not hasattr(self, "lumin"):
            self.lumin = df_lum
        else:
            self.lumin = pd.concat([self.lumin, df_lum], ignore_index=True)

        return self.lumin
    
class SpectrumManager:
    """
    Manage a collection of :py:class:`~xsnap.spectrum.SpectrumFit` objects for batch analysis and plotting.

    Attributes
    -----------
        specs : list[dict]
            Each entry is ``{'spec': SpectrumFit, 'instr': str}``.
        tExplosion : float or None
            Supernova time of explosion in MJD.
        fluxes : dict or None
            Combined flux DataFrames, keyed by 'absorbed'/'unabsorbed'.
        counts : pandas.DataFrame or None
            Combined count-rate DataFrame.
        lumin : pandas.DataFrame or None
            Combined luminosity DataFrame.
        params : pandas.DataFrame or None
            Combined parameters DataFrame.
    """
    def __init__(self, specs=None, tExplosion=None):
        """
        Initialize of the :py:class:`~xsnap.spectrum.SpectrumManager` class

        Parameters
        -----------
        specs : SpectrumFit or array[tuple(SpectrumFit, str)], optional
            Collection of :py:class:`~xsnap.spectrum.SpectrumFit` objects, can be parsed singularly
            or an iterable of (:py:class:`~xsnap.spectrum.SpectrumFit`, instrument) pairs. Defaults to ``None``.
        tExplosion : float, optional
            Supernova time of explosion in MJD to override or set. Defaults to ``None``.
        """
        self.specs = []  # list of {'spec': SpectrumFit, 'instr': str}
        self.tExplosion = tExplosion
        self.fluxes = None
        self.counts = None
        self.lumin = None
        self.params = None
        if specs is not None:
            self.load(specs)
            
    def _read_instr_from_header(self, pha_file: str):
        """Best-effort to guess the instrument label from a PHA file."""
        for ext in (1, 0):             # ext-1 first, then ext-0 fallback
            try:
                hdr = fits.getheader(pha_file, ext=ext)
                tele, instr = (hdr.get(k, '').strip() for k in ('TELESCOP', 'INSTRUME'))
                if tele.upper() == 'SWIFT' and instr:
                    return f'{tele} {instr}'
                return tele or instr or None
            except Exception:
                continue        # keep looping; on final failure we return None
        return None
    
    def _add_common_cols(self, df: pd.DataFrame, obs_time, obs_err,
                         instr: str, t_since_exp, t_err):
        """
        Append timing & instrument columns (including errors)
        """
        return df.assign(obs_time=obs_time,
                        obs_time_err=obs_err,
                        time_since_explosion=t_since_exp,
                        time_since_explosion_err=t_err,
                        instrument=instr)
        
    def _merge_df(self, target: pd.DataFrame | None, incoming: pd.DataFrame) -> pd.DataFrame:
        return incoming if target is None else pd.concat([target, incoming], ignore_index=True)

    def load(self, specs, instrument=None):
        """
        Load one or many :py:class:`~xsnap.spectrum.SpectrumFit` instances, with optional instrument labels.

        Parameters
        -----------
            specs : SpectrumFit or array[tuple(SpectrumFit, str)], optional
                Collection of :py:class:`~xsnap.spectrum.SpectrumFit` objects, can be parsed singularly
                or an iterable of (:py:class:`~xsnap.spectrum.SpectrumFit`, instrument) pairs. Defaults to ``None``.
            instrument : str, optional
                Instrument label to use if not provided per-:py:class:`~xsnap.spectrum.SpectrumFit`.

        Returns
        --------
            *self* : SpectrumManager
                The same :py:class:`~xsnap.spectrum.SpectrumManager` with new specs appended.

        Raises
        ------
            TypeError
                If ``specs`` is not a :py:class:`~xsnap.spectrum.SpectrumFit` or valid iterable.
            ValueError
                If an instrument label cannot be determined.
        """
        from collections.abc import Iterable
        items = []
        if isinstance(specs, SpectrumFit):
            items = [(specs, instrument)]
        elif isinstance(specs, Iterable) and not isinstance(specs, (str, SpectrumFit)):
            for entry in specs:
                if isinstance(entry, SpectrumFit):
                    items.append((entry, instrument))
                elif isinstance(entry, tuple) and len(entry) == 2 and isinstance(entry[0], SpectrumFit):
                    items.append((entry[0], entry[1]))
                else:
                    raise TypeError(
                        "Iterable entries must be SpectrumFit or (SpectrumFit, instr) tuples"
                    )
        else:
            raise TypeError("load() requires a SpectrumFit, or iterable thereof")
        
        for spec, instr in items:
            # Get instrument label from header if not provided
            if instr is None:
                instr = self._read_instr_from_header(spec.pha[0])
            if instr is None:
                raise ValueError("Instrument label required for SpectrumFit instance.")
            
            # ---------- explosion time ----------
            if self.tExplosion is None and spec.tExplosion is not None:
                self.tExplosion = spec.tExplosion

            # ---------- common timing values ----------
            obstime        = spec.obstime['obs_time'].iloc[0]
            obs_err        = spec.obstime['obs_time_err'].iloc[0]
            if 'time_since_explosion' in spec.obstime.columns:
                t_since_exp = spec.obstime['time_since_explosion'].iloc[0]
                t_err       = spec.obstime['time_since_explosion_err'].iloc[0]
            else:
                t_since_exp = (obstime - self.tExplosion) if self.tExplosion else None
                t_err       = obs_err if self.tExplosion else None

            # ---------- data frames ----------
            if spec.fluxes is not None:
                for kind, df_flux in spec.fluxes.items():
                    df_flux     = self._add_common_cols(df_flux.copy(), obstime, obs_err, instr, t_since_exp, t_err)
                    self.fluxes = self.fluxes or {}                  # lazy-init dict
                    self.fluxes[kind] = self._merge_df(self.fluxes.get(kind), df_flux)

            if spec.counts is not None:
                self.counts = self._merge_df(self.counts, self._add_common_cols(spec.counts.copy(), obstime, obs_err, instr, t_since_exp, t_err))
            
            if spec.lumin is not None:
                self.lumin  = self._merge_df(self.lumin,  self._add_common_cols(spec.lumin.copy(),  obstime, obs_err, instr, t_since_exp, t_err))
                
            if spec.params is not None:
                self.params = self._merge_df(self.params, self._add_common_cols(spec.params.copy(), obstime, obs_err, instr, t_since_exp, t_err))

        
            # Book keeping
            self.specs.append({'spec': spec, 'instr': instr})    
                
        return self
    
    def clear(self):
        """Clear all loaded SpectrumFit objects."""
        self.specs = []
        self.tExplosion = None
        self.fluxes = None
        self.counts = None
        self.lumin = None
        self.params = None
        return self

    def plot_flux(self, scatter: bool = True, log: bool = True):
        """
        Plot flux light curves for each model and kind ('absorbed', 'unabsorbed').

        Parameters
        -----------
            scatter : bool, optional
                If ``True``, use scatter markers; otherwise lines. Defaults to ``True``.
            log : bool, optional
                If ``True``, set y-axis to log scale. Defaults to ``True``.

        Returns:
            Flux light curve plots: dict[tuple(str, str), matplotlib.figure.Figure]
                Mapping (model: str, kind: 'absorbed' | 'unabsorbed') → matplotlib.figure.Figure.
                
                Data are grouped and labeled by instruments.
        """
        rows = []
        for entry in self.specs:
            spec = entry['spec']; instr = entry['instr']
            df_time = spec.obstime if spec.obstime is not None else spec.getTime(self.tExplosion or spec.tExplosion)
            for kind, df_flux in spec.fluxes.items():
                tmp = df_flux.merge(df_time, on='data')
                tmp['type'] = kind; tmp['instrument'] = instr
                rows.append(tmp)
        df_all = pd.concat(rows, ignore_index=True)

        figs = {}
        for kind in ['absorbed', 'unabsorbed']:
            for model, dfkm in df_all[df_all['type'] == kind].groupby('model'):
                fig, ax = plt.subplots()
                for instr, grp in dfkm.groupby('instrument'):
                    # determine x
                    if 'time_since_explosion' in grp.columns:
                        x = grp['time_since_explosion']
                        xerr  = grp['time_since_explosion_err'] if 'time_since_explosion_err' in grp else None
                    elif self.tExplosion is not None:
                        x = grp['obs_time'] - self.tExplosion
                        xerr = grp['obs_time_err'] 
                    else:
                        x = grp['obs_time']
                        xerr = grp['obs_time_err']
                    y = grp['flux']
                    if scatter:
                        yerr = (grp['lo_flux_err'].values, grp['hi_flux_err'].values)
                        ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='o', label=instr)
                    else:
                        ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='-o-', label=instr)
                fig.suptitle(f"{kind.capitalize()} {model} flux")
                xlabel = 'Time since explosion (days)' if ('time_since_explosion' in df_all.columns or self.tExplosion is not None) else 'MJD'
                ax.set_xlabel(xlabel)
                if log: ax.set_yscale('log')
                ax.set_ylabel(r'Flux $(\text{erg } \text{s}^{-1} \text{ cm}^{-2})$')
                ax.legend(); ax.grid(True, which="both", ls=":")
                plt.show()
                figs[(model, kind)] = fig
        return figs

    def plot_lumin(self, scatter: bool = True, log: bool = True):
        """
        Plot luminosity light curves for each model.

        Parameters
        -----------
            scatter : bool, optional
                If ``True``, use scatter markers; otherwise lines. Defaults to ``True``.
            log : bool, optional
                If ``True``, set y-axis to log scale. Defaults to ``True``.

        Returns
        --------
            Luminosity light curve plots : dict[str, matplotlib.figure.Figure]
                Mapping model: str → matplotlib.figure.Figure.
                
                Data are grouped and labeled by instruments.
        """

        rows = []
        for entry in self.specs:
            spec = entry['spec']; instr = entry['instr']
            df_time = spec.obstime if spec.obstime is not None else spec.getTime(self.tExplosion or spec.tExplosion)
            df_lum = spec.lumin.copy()
            if 'data' not in df_lum.columns:
                df_lum['data'] = df_time['data'].values
            tmp = df_lum.merge(df_time, on='data')
            tmp['instrument'] = instr; rows.append(tmp)
        df_all = pd.concat(rows, ignore_index=True)

        figs = {}
        for model, grp in df_all.groupby('model'):
            fig, ax = plt.subplots()
            for instr, g2 in grp.groupby('instrument'):
                if 'time_since_explosion' in g2.columns:
                    x = g2['time_since_explosion']
                    xerr  = g2['time_since_explosion_err'] if 'time_since_explosion_err' in grp else None
                elif self.tExplosion is not None:
                    x = g2['obs_time'] - self.tExplosion
                    xerr = g2['obs_time_err']
                else:
                    x = g2['obs_time']
                    xerr = g2['obs_time_err']
                y = g2['lumin']
                if scatter:
                    yerr = (g2['lo_lumin_err'].values, g2['hi_lumin_err'].values)
                    ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='o', label=instr)
                else:
                   ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='-o-', label=instr)
            fig.suptitle(f"Luminosity {model} light curve")
            xlabel = 'Time since explosion (days)' if ('time_since_explosion' in df_all.columns or self.tExplosion is not None) else 'MJD'
            ax.set_xlabel(xlabel)
            if log: ax.set_yscale('log')
            ax.set_ylabel(r'Luminosity $(\text{erg } \text{s}^{-1})$')
            ax.legend(); ax.grid(True, which="both", ls=":")
            plt.show()
            figs[model] = fig
        return figs
    
    def plot_phot(self, scatter: bool = True, log: bool = True):
        """
        Plot photon-flux light curves for each model and kind.

        Parameters
        -----------
            scatter : bool, optional
                If ``True``, use scatter markers; otherwise lines. Defaults to ``True``.
            log : bool, optional
                If ``True``, set y-axis to log scale. Defaults to ``True``.

        Returns
        -------
            Photon flux light curve plots: dict[tuple(str, str), matplotlib.figure.Figure]
                Mapping (model: str, kind: 'absorbed' | 'unabsorbed') → matplotlib.figure.Figure.
                
                Data are grouped and labeled by instruments.
        """
        rows = []
        for entry in self.specs:
            spec = entry['spec']; instr = entry['instr']
            df_time = spec.obstime if spec.obstime is not None else spec.getTime(self.tExplosion or spec.tExplosion)
            for kind, df_flux in spec.fluxes.items():
                tmp = df_flux.merge(df_time, on='data')
                tmp['type'] = kind; tmp['instrument'] = instr
                rows.append(tmp)
        df_all = pd.concat(rows, ignore_index=True)

        figs = {}
        for kind in ['absorbed', 'unabsorbed']:
            for model, dfkm in df_all[df_all['type'] == kind].groupby('model'):
                fig, ax = plt.subplots()
                for instr, grp in dfkm.groupby('instrument'):
                    # determine x
                    if 'time_since_explosion' in grp.columns:
                        x = grp['time_since_explosion']
                        xerr  = grp['time_since_explosion_err'] if 'time_since_explosion_err' in grp else None
                    elif self.tExplosion is not None:
                        x = grp['obs_time'] - self.tExplosion
                        xerr = grp['obs_time_err'] 
                    else:
                        x = grp['obs_time']
                        x = grp['obs_time_err']
                    y = grp['phot']
                    if scatter:
                        yerr = (grp['lo_phot_err'].values, grp['hi_phot_err'].values)
                        ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='o', label=instr)
                    else:
                        ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='-o-', label=instr)
                fig.suptitle(f"{kind.capitalize()} {model} phot")
                xlabel = 'Time since explosion (days)' if ('time_since_explosion' in df_all.columns or self.tExplosion is not None) else 'MJD'
                ax.set_xlabel(xlabel)
                if log: ax.set_yscale('log')
                ax.set_ylabel(r'Photon Flux $(\text{phot } \text{s}^{-1} \text{ cm}^{-2})$')
                ax.legend(); ax.grid(True, which="both", ls=":")
                plt.show()
                figs[(model, kind)] = fig
        return figs
        
    def plot_counts(self, scatter: bool = True, log: bool = True):
        """
        Plot count-rate light curves for each model; one trace per instrument.

        Parameters
        -----------
            scatter : bool, optional
                If ``True``, use scatter markers; otherwise lines. Defaults to ``True``.
            log : bool, optional
                If ``True``, set y-axis to log scale. Defaults to ``True``.

        Returns
        --------
            Count rate light curve plots : dict[str, matplotlib.figure.Figure]
                Mapping model: str → matplotlib.figure.Figure.
                
                Data are grouped and labeled by instruments.
        """
        # ---------- collect & merge ----------
        rows = []
        for entry in self.specs:
            spec  = entry['spec']
            instr = entry['instr']

            df_time = (spec.obstime if spec.obstime is not None else
                    spec.getTime(self.tExplosion or spec.tExplosion))

            df_cnt  = spec.counts.copy()
            df_cnt  = df_cnt.merge(df_time, on="data")
            df_cnt["instrument"] = instr
            rows.append(df_cnt)

        df_all = pd.concat(rows, ignore_index=True)

        # ---------- plot ----------
        figs = {}
        for model, grp_model in df_all.groupby("model"):
            fig, ax = plt.subplots()

            for instr, g2 in grp_model.groupby("instrument"):
                # --- x + x-error ---
                if "time_since_explosion" in g2.columns:
                    x     = g2["time_since_explosion"]
                    xerr  = g2.get("time_since_explosion_err", None)
                elif self.tExplosion is not None:
                    x     = g2["obs_time"] - self.tExplosion
                    xerr  = g2.get("obs_time_err", None)
                else:
                    x     = g2["obs_time"]
                    xerr  = g2.get("obs_time_err", None)

                # --- y + y-error ---
                y    = g2["net_rate"]
                yerr = g2["net_err"]

                # --- draw ---
                if scatter:
                    ax.errorbar(x, y, yerr=yerr, xerr=xerr,
                                fmt="o", capsize=2, label=instr)
                else:
                    ax.errorbar(x, y, yerr=yerr, xerr=xerr,
                                fmt="-o", capsize=2, label=instr)

            # cosmetics
            ax.set_title(f"Count-rate light curve ({model})")
            xlabel = ("Time since explosion (days)"
                    if ("time_since_explosion" in df_all.columns or
                        self.tExplosion is not None) else "MJD")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r"Count rate (cts s$^{-1}$)")
            if log:
                ax.set_yscale("log")
            ax.grid(True, which="both", ls=":")
            ax.legend()

            plt.show()
            figs[model] = fig

        return figs

    def plot_params(self, scatter: bool = True, log: bool = True):
        """
        Plot fit parameter evolution vs. time for each model and parameter.

        Parameters
        -----------
            scatter : bool, optional
                If ``True``, use scatter markers; otherwise lines. Defaults to ``True``.
            log : bool, optional
                If ``True``, set y-axis to log scale. Defaults to ``True``.

        Returns:
            Parameter evolution plots: dict[tuple(str, str), matplotlib.figure.Figure]
                Mapping (model: str, parameter: str) → matplotlib.figure.Figure.
                
                Data are grouped and labeled by instruments.
        """
        import numpy as np

        # 1) Merge every spec.params with its time table
        rows = []
        for entry in self.specs:
            spec  = entry['spec']
            instr = entry['instr']
            dfp   = spec.params.copy()
            df_time = (
                spec.obstime
                if spec.obstime is not None
                else spec.getTime(self.tExplosion or spec.tExplosion)
            )
            dfp = dfp.merge(df_time, on='data')
            dfp['instrument'] = instr
            rows.append(dfp)
        df_all = pd.concat(rows, ignore_index=True)

        figs = {}

        # identify parameter columns (exclude the id‐vars and error columns)
        error_like = lambda c: (
                c.endswith('_err') or            # timing errors
                c.startswith('lo_') or           # our own lower errors
                c.startswith('hi_')              # our own upper errors
        )
        id_vars = {'data','model','obs_time','time_since_explosion',
                'obs_time_err','time_since_explosion_err','instrument'}

        param_cols = [c for c in df_all.columns
                    if c not in id_vars and not error_like(c)]

        # now loop by model → then by parameter
        for model in df_all['model'].unique():
            grp_model = df_all[df_all['model'] == model]
            # only keep parameters that appear (non‐all‐NaN) in this model
            model_params = [p for p in param_cols
                            if grp_model[p].notna().any()]

            for param in model_params:
                # choose which columns actually exist
                cols = ['instrument','data','obs_time']
                if 'time_since_explosion' in grp_model.columns:
                    cols.append('time_since_explosion')
                cols += [param, f'lo_{param}_err', f'hi_{param}_err']

                grp = grp_model[cols].copy()

                fig, ax = plt.subplots()
                # one curve per instrument
                for instr, g2 in grp.groupby('instrument'):
                    # x-axis: prefer time_since_explosion if available
                    if 'time_since_explosion' in grp.columns:
                        x = g2['time_since_explosion']
                    elif self.tExplosion is not None:
                        x = g2['obs_time'] - self.tExplosion
                    else:
                        x = g2['obs_time']

                    y = g2[param]
                    if scatter:
                        lo = g2[f'lo_{param}_err']
                        hi = g2[f'hi_{param}_err']
                        yerr = np.vstack([lo.values, hi.values])
                        ax.errorbar(x, y, yerr=yerr, fmt='o', label=instr)
                    else:
                        ax.plot(x, y, 'o-', label=instr)

                # title & labels
                fig.suptitle(f"{param} ({model}) parameter evolution")
                xlabel = ('Time since explosion (days)'
                          if 'time_since_explosion' in grp.columns or self.tExplosion is not None
                          else 'MJD')
                ax.set_xlabel(xlabel)

                # only log‐scale if _every_ y > 0
                if log and (grp[param] > 0).all():
                    ax.set_yscale('log')

                ax.set_ylabel(param)
                ax.legend()
                ax.grid(True, which="both", ls=":")
                plt.show()

                figs[(model, param)] = fig

        return figs