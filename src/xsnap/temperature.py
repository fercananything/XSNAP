from __future__ import annotations
import numpy as np
from astropy.io import fits 
from astropy.time import Time
import matplotlib.pyplot as plt
import pandas as pd
import emcee, os
from pathlib import Path
from ._fitting import fit_powerlaw_asymmetric, compute_chi2_powerlaw, predict_with_errors

class TemperatureEstimator:
    """
    Estimate temperature evolution from a supernova explosion 
    given a set of temperature evolution data or power-law parameters.
    
    If given a set of temperature evolution data, the :py:class:`~xsnap.temperature.TemperatureEstimator` will
    fit the data in a power-law model and get the best-fit parameters.
        
    Attributes
    ------------
        tExplosion : float
            Supernova time of explosion in MJD.
        norm : float
            Best-fit normalization for the power-law model.
        lo_norm_err : float
            Lower uncertainty of the normalization.
        hi_norm_err : float
            Upper uncertainty of the normalization.
        exp : float
            Best-fit exponent for the power-law model.
        lo_exp_err : float
            Lower uncertainty of the exponent.
        hi_exp_err : float
            Upper uncertainty of the exponent.
        chi2_red : float
            The reduced chi-squared of the fit.
        temperatures : pandas.DataFrame
            Estimated temperatures wrapped in a DataFrame table with columns
            ``['time_since_explosion', 'lo_time_err', 'hi_time_err', 'temperature', 'lo_temp_err', 'hi_temp_err']``
    """
    
    def __init__(self, tExplosion: float | None = None, norm: float | None = None, 
                 exponent: float | None = None):
        """Initialization of the :py:class:`~xspec.temperature.TemperatureEstimator` class

        Parameters
        ----------
            tExplosion : float, optional
                Supernova time of explosion in MJD. If ``None``, must be set before estimating. Defaults to ``None``.
            norm : float, optional
                Normalization for the power-law model (if known). Defaults to ``None``.
            exponent : float, optional
                Exponent for the power-law model (if known). Defaults to ``None``.
        """
        
        self.tExplosion = tExplosion
        self.norm = norm
        self.lo_norm_err = 0
        self.hi_norm_err = 0
        self.exp = exponent
        self.lo_exp_err = 0
        self.hi_exp_err = 0
        self.chi2_red = None
        self.temperatures = pd.DataFrame({'time_since_explosion': [], 'lo_time_err': [], 'hi_time_err': [],
                                          'temperature': [], 'lo_temp_err': [], 'hi_temp_err': []})
        
    def clear(self):
        """
        Reset all fitted parameters and clear the temperature table.
        """
        self.tExplosion = None
        self.norm = None
        self.lo_norm_err = 0
        self.hi_norm_err = 0
        self.exp = None
        self.lo_exp_err = 0
        self.hi_exp_err = 0
        self.chi2_red = None
        self.temperatures = pd.DataFrame({'time_since_explosion': [], 'lo_time_err': [], 'hi_time_err': [],
                                          'temperature': [], 'lo_temp_err': [], 'hi_temp_err': []})
        
    def compute_pl_fit(self, time_since_explosion, temperature, temp_err_lo=None, temp_err_hi=None,
                        time_err_lo=None, time_err_hi=None,
                        nwalkers=200, nsteps=6000, nburn=1000,
                        show_plots=True):
        
        """
        Fit a power-law :math:`T(t) = \\mathrm{norm} \\times t^{\mathrm{exp}}` to the data, 
        via Markov chain Monte-Carlo (MCMC).

        Parameters
        ----------
            time_since_explosion : array_like
                Times since explosion in days.
            temperature : array_like
                Observed/fitted temperatures data..
            temp_err_lo : array_like, optional
                Lower uncertainty on temperature. Defaults to ``None``.
            temp_err_hi : array_like, optional
                Upper uncertainty on temperature. Defaults to ``None``.
            time_err_lo : array_like, optional
                Lower uncertainty on time. Defaults to ``None``.
            time_err_hi : array_like, optional
                Upper uncertainty on time. Defaults to ``None``.
            nwalkers : int, optional
                Number of MCMC walkers. Defaults to ``200``.
            nsteps : int, optional
                Number of MCMC steps per walker. Defaults to ``6000``.
            nburn : int, optional
                Number of burn-in steps to discard. Defaults to ``1000``.
            show_plots (bool, optional
                If ``True``, display fitted and residual plots. Defaults to ``True``.

        Returns
        ---------
            *self* : TemperatureEstimator
                The same :py:class:`~xsnap.temperature.TemperatureEstimator` class with 
                :py:attr:`norm`, :py:attr:`exp`, :py:attr:`lo_norm_err`, :py:attr:`hi_norm_err`, 
                :py:attr:`lo_exp_err`, :py:attr:`hi_exp_err`, and :py:attr:`chi2_red` all set.
        """
        
        norm, lo_norm_err, hi_norm_err, exp, lo_exp_err, hi_exp_err, chain = fit_powerlaw_asymmetric(time_since_explosion, temperature, 
                                                                                                temp_err_lo, temp_err_hi,
                                                                                                time_err_lo, time_err_hi,
                                                                                                nwalkers, nsteps, nburn,
                                                                                                show_plots)
            
        self.norm = norm
        self.lo_norm_err = lo_norm_err
        self.hi_norm_err = hi_norm_err
        self.exp = exp
        self.lo_exp_err = lo_exp_err
        self.hi_exp_err = hi_exp_err
            
        chi2, chi2_red = compute_chi2_powerlaw(time_since_explosion, temperature, 
                                                    temp_err_lo, temp_err_hi,
                                                    A=norm, A_err_lo=lo_norm_err, A_err_hi=hi_norm_err,
                                                    k=exp, k_err_lo=lo_exp_err, k_err_hi=hi_exp_err,
                                                    xlo=time_err_lo, xhi=time_err_hi,
                                                    dof=None, plot_resid=show_plots
                                                )
        self.chi2_red = chi2_red
            
        return self
    
    def _getTime(self, file):
        """
        Read a FITS file header to compute the observation time relative to explosion.

        Args:
            file (str): Path to the FITS file.

        Returns:
            tuple:
                t_days (float): Time since explosion in days.
                t_err (float): Uncertainty on that time.

        Raises:
            RuntimeError: If `tExplosion` is not set or file does not exist.
            KeyError: If no valid DATE-OBS or MJD-OBS keywords are found.
            ValueError: If the observation precedes the explosion time.
        """
        if self.tExplosion is None:
            raise RuntimeError("Please input the time of explosion in MJD!")
        
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
            
        # --- compute time since explosion
        t_days = mjd_obs - self.tExplosion
        
        if t_days < 0:
            raise ValueError("Observation date is before explosion date.")
        
        t_err = mjd_err
        
        return t_days, t_err
        
    def estimate(self, time_since_explosion = None, t_err_lo = 0, t_err_hi = 0, 
                 file: str | None = None, files = None, tExplosion = None):
        """
        Estimate temperature(s) at given time(s) or from FITS file(s) and append to :py:attr:`temperatures`.

        Parameters
        -----------
            time_since_explosion : float or array_like, optional
                Time(s) since explosion in days.
            t_err_lo : float or array_like, optional
                Lower time error(s). Defaults to 0.
            t_err_hi : float or array_like, optional
                Upper time error(s). Defaults to 0.
            file : str, optional
                Single FITS file to read time from.
            files : list of str, optional
                List of FITS files to read times from.
            tExplosion : float, optional
                Explosion time (MJD), if not already set.

        Returns
        ---------
            Temperatures : pandas.DataFrame
                The :py:attr:`temperatures` table wrapped in DataFrame with columns
                ``['time_since_explosion', 'lo_time_err', 'hi_time_err', 'temperature', 'lo_temp_err', 'hi_temp_err']``.

        Raises
        -------
            RuntimeError
                If neither ``time_since_explosion`` nor ``file`` or ``files`` provided,
                or if model parameters have not been fitted,
                or if ``tExplosion`` is still unset when reading FITS.
        """
        
        if time_since_explosion is None and file is None and files is None:
            raise RuntimeError("Please enter something to estimate!")
        
        if self.norm is None or self.exp is None:
            raise RuntimeError("Please compute_pl_fit first!")
        
        A_m = self.norm
        A_err_lo = self.lo_norm_err
        A_err_hi = self.hi_norm_err
        k_m = self.exp
        k_err_lo = self.lo_exp_err 
        k_err_hi = self.hi_exp_err
        # A_err_lo, A_err_hi, k_err_lo, k_err_hi = 0,0,0,0
        
        if self.tExplosion is None:
            self.tExplosion = tExplosion
        
        if file is None and files is None:
            x_m = time_since_explosion, x_err_lo = t_err_lo, x_err_hi = t_err_hi
            y_m, y_err_lo, y_err_hi = predict_with_errors(
                A_m, A_err_lo, A_err_hi,
                k_m, k_err_lo, k_err_hi,
                x_m, x_err_lo = x_err_lo, x_err_hi = x_err_hi
            )
            row = {
                'time_since_explosion': x_m, 
                'lo_time_err': x_err_lo, 
                'hi_time_err': x_err_hi,
                'temperature': y_m,
                'lo_temp_err': y_err_lo, 
                'hi_temp_err': y_err_hi
            }
            self.temperatures = pd.concat(
                [self.temperatures, pd.DataFrame([row])],
                ignore_index=True
            )
            return self.temperatures
        
        if self.tExplosion is None:
            raise RuntimeError("Please input the time of explosion in MJD!")
        
        def get_temp(file):
            file = os.path.abspath(file)
            t_days, t_err = self._getTime(file)
            y_m, y_err_lo, y_err_hi = predict_with_errors(
                A_m, A_err_lo, A_err_hi,
                k_m, k_err_lo, k_err_hi,
                t_days, 
                x_err_lo = t_err, x_err_hi = t_err
            )
            
            row = {
                'time_since_explosion': t_days, 
                'lo_time_err': t_err, 
                'hi_time_err': t_err,
                'temperature': y_m,
                'lo_temp_err': y_err_lo, 
                'hi_temp_err': y_err_hi
            }
            self.temperatures = pd.concat(
                [self.temperatures, pd.DataFrame([row])],
                ignore_index=True
            )
            
        
        if time_since_explosion is None and files is None:
            get_temp(file)
            
            return self.temperatures
        
        if time_since_explosion is None and file is None:
            for file in files:
                get_temp(file)
                
            return self.temperatures