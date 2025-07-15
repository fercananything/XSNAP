#!/usr/bin/env python3
"""
analysis.py - to analysis CSM shock interactions

Provides CSMAnalysis, which loads SpectrumManager class and
get the density & mass-loss rate and plot them vs radius
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .spectrum import SpectrumFit, SpectrumManager
import astropy.units as u 
from ._fitting import fit_powerlaw_asymmetric, compute_chi2_powerlaw, fit_powerlaw_simple, compute_chi2_powerlaw_simple
from scipy.optimize import curve_fit

# Constants
C_LIGHT_KM_S = 299792.458 * u.km / u.s
MPC_TO_CM = (1 * u.Mpc).to(u.cm).value
v_wind = 20 * u.km / u.s
v_shock = 10000 * u.km / u.s
m_p = 1.67262192e-27 * u.kg

class CSMAnalysis:
    """
    Physics-side for Supernova CSM analysis 
    Fitting L(t) ∝ t^x, deriving CSM densities,
    and mass-loss rates from a SpectrumManager.lumin table.
    
    Args:
        manager (SpectrumManager, optional): Pre-populated SpectrumManager instance.

    Attributes:
        manager (SpectrumManager or None): Input manager spectra analysis.
        distance (float or None): Source distance in Mpc.
        r_shock (pandas.DataFrame or None): Shock radius.
        times (pandas.DataFrame or None): Times since explosion.
        fit_lumin_params (pandas.DataFrame or None): Luminosity fit parameters.
        fit_temp_params (pandas.DataFrame or None): Temperature fit parameters.
        fit_density_params (pandas.DataFrame or None): Density fit parameters.
        densities (pandas.DataFrame or None): Computed shock densities.
        mass_loss_rate (pandas.DataFrame): Computed mass-loss rates.
    """
    def __init__(self, manager=None):
        self.manager = None
        self.distance = None
        self.r_shock = None
        self.times = None
        self.fit_lumin_params = None
        self.fit_temp_params = None
        self.fit_density_params = None
        self.densities = None
        self.mass_loss_rate = pd.DataFrame({'fit': [], 'm_dot': [], 'lo_m_dot_err': [], 'hi_m_dot_err': []})
        
        if manager is not None:
            self.load(manager)
        
    def _dist_from_z(self, z, H0):
        """Approximate luminosity distance D(z) [Mpc] (non-rel relativistic)."""
        if z < 0.1:
            v = z * C_LIGHT_KM_S
        else:
            v = C_LIGHT_KM_S * ((1 + z) ** 2 - 1) / ((1 + z) ** 2 + 1)
        return (v / H0).to_value(u.Mpc)
    
    def _distance_and_errors(self, distance: float | None, lo_dist_err: float | None, hi_dist_err: float | None,
                              z: float | None, lo_z_err: float | None, hi_z_err: float | None, H0):
        """
        Resolve distance [Mpc] plus asymmetric errors.

        If only (z, z_err) are given → propagate σ_D = σ_z * dD/dz ≈ σ_z * c/H0.
        """
        if distance is None and z is None:
            raise ValueError("Either distance or redshift must be provided.")

        if distance is None:                                  # use redshift
            distance = self._dist_from_z(z, H0)
            lo_dist_err = (
                self._dist_from_z(max(z - (lo_z_err or 0), 0), H0) - distance
                if lo_z_err is not None
                else 0.0
            )
            lo_dist_err = abs(lo_dist_err)
            hi_dist_err = (
                self._dist_from_z(z + (hi_z_err or 0), H0) - distance
                if hi_z_err is not None
                else 0.0
            )
        else:                                                 # explicit distance
            lo_dist_err = float(lo_dist_err or 0.0)
            hi_dist_err = float(hi_dist_err or 0.0)

        return float(distance), float(lo_dist_err), float(hi_dist_err)

    def load(self, manager: SpectrumManager, distance=None, lo_dist_err=None, hi_dist_err=None,
             z=None, lo_z_err=None, hi_z_err=None, v_shock: float = v_shock.value, H0 = 70 * u.km / u.s / u.Mpc):
        """
        Load manager luminosity & parameter tables and compute shock radii.

        Args:
            manager (SpectrumManager): Must contain one 'bremss' model in `.lumin` & `.params`.
            distance (float, optional): Distance in Mpc. Required if `z` is None.
            z (float, optional): Redshift to compute distance. Required if `distance` is None.
            lo_dist_err (float, optional): Lower uncertainty of distance in Mpc.
            hi_dist_err (float, optional): Upper uncertainty of distance in Mpc.
            lo_z_err (float, optional): Lower uncertainty of redshift.
            hi_z_err (float, optional): Upper uncertainty of redshift.
            H0 (float, optional): Hubble constant (km/s/Mpc). Defaults to 70.0 km/s/Mpc.
            v_shock (float, optional): Shock velocity in km/s. Defaults to 10000 km/s.

        Returns:
            CSMAnalysis: self, with `times` and `r_shock` populated.

        Raises:
            RuntimeError: If more than one bremss model is present.
            ValueError: If required tables/columns are missing.
        """
        self.manager = manager
        
        m_bremss = manager.lumin['model'].str.contains('bremss', case=False, na=False)
        df       = manager.lumin.loc[m_bremss]
        if len(df['model'].unique()) != 1:
            raise RuntimeError("Please to use only one bremss model")
        
        if distance is not None and z is not None:
            D_Mpc, d_lo, d_hi = self._distance_and_errors(
                distance, lo_dist_err, hi_dist_err,
                z, lo_z_err, hi_z_err, H0
            )
            self.distance = [D_Mpc, d_lo, d_hi]
        
        # Check if params and lumin have been populated and whether they have bremss model
        params_ready = (
            (manager.params is not None) and
            manager.params['model'].str.contains('bremss', case=False, na=False).any()
        )

        lumin_ready = (
            (manager.lumin is not None) and
            manager.lumin['model'].str.contains('bremss', case=False, na=False).any()
        )

        if params_ready and lumin_ready:
            temp_df = manager.lumin[['time_since_explosion', 'time_since_explosion_err']].copy()
            temp_df = temp_df.drop_duplicates().sort_values('time_since_explosion').reset_index(drop=True)  
            self.times = temp_df
            temp_time = (np.array(self.times['time_since_explosion']) * u.day).to(u.s).value
            temp_time_err = (np.array(self.times['time_since_explosion_err']) * u.day).to(u.s).value
            r_shock = temp_time * (v_shock * u.km / u.s).to(u.cm / u.s).value
            r_shock_err = temp_time_err * (v_shock * u.km / u.s).to(u.cm / u.s).value
            self.r_shock = pd.DataFrame({
                    'r_shock': r_shock,
                    'r_shock_err': r_shock_err
                })
            return self
                    
        raise ValueError("Please populate the luminosity and parameters DataFrame first and make sure they have 'bremss' in their 'model' column") 
        
    def clear(self):
        """
        Reset all loaded data and computed results.
        """
        self.manager = None
        self.distance = None
        self.r_shock = None
        self.times = None
        self.fit_lumin_params = None
        self.fit_temp_params = None
        self.fit_density_params = None
        self.densities = None
        self.mass_loss_rate = pd.DataFrame({'fit': [], 'm_dot': [], 'lo_m_dot_err': [], 'hi_m_dot_err': []})

    def fit_lumin_mcmc(self, nwalkers=500, nsteps=10000, nburn=2000, show_plots=True):
        """
        Fit luminosity L ∝ t^x using MCMC.

        Args:
            nwalkers (int, optional): Number of MCMC walkers. Defaults to 500.
            nsteps (int, optional): Number of steps per walker. Defaults to 10000.
            nburn (int, optional): Number of burn-in steps. Defaults to 2000.
            show_plots (bool, optional): If True, display diagnostic plots.

        Returns:
            pandas.DataFrame: columns ['model','fit','norm','lo_norm_err','hi_norm_err',
                                       'exp','lo_exp_err','hi_exp_err','ndata'].
        """
        df = self.manager.lumin.copy()
        m_bremss = df['model'].str.contains('bremss', case=False, na=False)
        df       = df.loc[m_bremss]
        models = df['model'].unique().tolist()

        records = []
        figs = {}

        for mod in models:
            sub = df[df['model'] == mod]
            # require positive times & luminosities
            mask = (sub['time_since_explosion'] > 0) & (sub['lumin'] > 0)
            sub = sub.loc[mask]
            if len(sub) < 2:
                continue
            
            t = sub['time_since_explosion']
            t_err = sub['time_since_explosion_err']
            L = sub['lumin']
            lo_L_err = sub['lo_lumin_err']
            hi_L_err = sub['hi_lumin_err']

            norm, lo_norm_err, hi_norm_err, exp, lo_exp_err, hi_exp_err, chain = fit_powerlaw_asymmetric(
                                                                                    t, L, lo_L_err, hi_L_err,
                                                                                    xerr_lo=t_err, xerr_hi=t_err,
                                                                                    nwalkers=nwalkers, nsteps=nsteps, 
                                                                                    nburn=nburn, show_plots=show_plots
                                                                                 )
            
            compute_chi2_powerlaw(
                t, L, lo_L_err, hi_L_err,
                norm, lo_norm_err, hi_norm_err,
                exp, lo_exp_err, hi_exp_err,
                xlo=t_err, xhi=t_err,
                dof=None, plot_resid=show_plots
            )

            records.append({
                'model':     mod,
                'fit':    'mcmc',
                'norm':       norm,
                'lo_norm_err': lo_norm_err,
                'hi_norm_err': hi_norm_err,
                'exp':        exp,
                'lo_exp_err':   lo_exp_err,
                'hi_exp_err':   hi_exp_err,
                'ndata':     len(sub)
            })

        df_fits = pd.DataFrame.from_records(
            records,
            columns=['model', 'fit', 'norm', 'lo_norm_err', 'hi_norm_err', 'exp',
                     'lo_exp_err', 'hi_exp_err', 'ndata']
        )

        if self.fit_lumin_params is None:
            self.fit_lumin_params = df_fits
        else:
            self.fit_lumin_params = pd.concat(
                [self.fit_lumin_params, df_fits],
                axis=0,
                ignore_index=True
            )
        return df_fits
    
    def fit_lumin_simple(self, plot_resid=True):
        """
        Fit luminosity L ∝ t^x using scipy.optimize.curve_fit.

        Args:
            plot_resid (bool, optional): If True, show residual plot. Defaults to True.

        Returns:
            pandas.DataFrame: same columns as `fit_lumin_mcmc` with fit='curve_fit'.
        """
        df = self.manager.lumin.copy()
        m_bremss = df['model'].str.contains('bremss', case=False, na=False)
        df       = df.loc[m_bremss]
        models = df['model'].unique().tolist()

        records = []
        figs = {}

        for mod in models:
            sub = df[df['model'] == mod]
            # require positive times & luminosities
            mask = (sub['time_since_explosion'] > 0) & (sub['lumin'] > 0)
            sub = sub.loc[mask]
            if len(sub) < 2:
                continue
            
            t = sub['time_since_explosion']
            t_err = sub['time_since_explosion_err']
            L = sub['lumin']
            lo_L_err = sub['lo_lumin_err']
            hi_L_err = sub['hi_lumin_err']

            norm, norm_err, exp, exp_err = fit_powerlaw_simple(t, L, yerr_lo=lo_L_err, yerr_hi=hi_L_err)
            
            chi2, chi2red = compute_chi2_powerlaw_simple(t, L, yerr_lo=lo_L_err, yerr_hi=hi_L_err, 
                                                         norm=norm, exp=exp, plot_resid=plot_resid)

            records.append({
                'model':     mod,
                'fit':    'curve_fit',
                'norm':       norm,
                'lo_norm_err': norm_err,
                'hi_norm_err': norm_err,
                'exp':        exp,
                'lo_exp_err':   exp_err,
                'hi_exp_err':   exp_err,
                'ndata':     len(sub)
            })

        df_fits = pd.DataFrame.from_records(
            records,
            columns=['model', 'fit', 'norm', 'lo_norm_err', 'hi_norm_err', 'exp',
                     'lo_exp_err', 'hi_exp_err', 'ndata']
        )

        if self.fit_lumin_params is None:
            self.fit_lumin_params = df_fits
        else:
            self.fit_lumin_params = pd.concat(
                [self.fit_lumin_params, df_fits],
                axis=0,
                ignore_index=True
            )
        return df_fits

    
    def fit_temp_mcmc(self, nwalkers=500, nsteps=10000, nburn=2000, show_plots=True):
        """
        Fit temperature kT ∝ t^x using MCMC from data in `.params`.

        Args:
            nwalkers (int, optional): Number of MCMC walkers. Defaults to 500.
            nsteps (int, optional): Number of steps per walker. Defaults to 10000.
            nburn (int, optional): Number of burn-in steps. Defaults to 2000.
            show_plots (bool, optional): If True, display diagnostic plots.

        Returns:
            pandas.DataFrame: columns ['model','fit','norm','lo_norm_err','hi_norm_err',
                                       'exp','lo_exp_err','hi_exp_err','ndata'].
        """
        df = self.manager.params.copy()

        # ─────────────────────────────────────────────────────────────
        #  keep only rows whose model string contains “bremss”
        # ─────────────────────────────────────────────────────────────
        m_bremss = df['model'].str.contains('bremss', case=False, na=False)
        df       = df.loc[m_bremss]

        if df.empty:
            raise ValueError("manager.params contains no rows with a 'bremss' model.")

        models = [df['model'].unique()[0]]
        
        records = []

        for mod in models:
            sub = df[df['model'] == mod]

            # require positive times & temperatures
            m_valid = (sub['time_since_explosion'] > 0) & (sub['bremss_kT'] > 0)
            sub     = sub.loc[m_valid]
            
            if len(sub) < 2:
                continue                       # not enough points to fit
            
            t = sub['time_since_explosion']
            t_err = sub['time_since_explosion_err']
            T = sub['bremss_kT']
            lo_T_err = sub['lo_bremss_kT_err']
            hi_T_err = sub['hi_bremss_kT_err']

            norm, lo_norm_err, hi_norm_err, exp, lo_exp_err, hi_exp_err, chain = fit_powerlaw_asymmetric(
                                                                                    t, T, lo_T_err, hi_T_err,
                                                                                    xerr_lo=t_err, xerr_hi=t_err,
                                                                                    nwalkers=nwalkers, nsteps=nsteps, 
                                                                                    nburn=nburn, show_plots=show_plots
                                                                                 )
            
            compute_chi2_powerlaw(
                t, T, lo_T_err, hi_T_err,
                norm, lo_norm_err, hi_norm_err,
                exp, lo_exp_err, hi_exp_err,
                xlo=t_err, xhi=t_err,
                dof=None, plot_resid=show_plots
            )

            records.append({
                'model':     mod,
                'fit':    'mcmc',
                'norm':       norm,
                'lo_norm_err': lo_norm_err,
                'hi_norm_err': hi_norm_err,
                'exp':        exp,
                'lo_exp_err':   lo_exp_err,
                'hi_exp_err':   hi_exp_err,
                'ndata':     len(sub)
            })

        if not records:
            raise ValueError("No valid ‘bremss’ rows with positive t and kT to fit.")

        df_fits = pd.DataFrame.from_records(
            records,
            columns=['model', 'fit', 'norm', 'lo_norm_err', 'hi_norm_err', 'exp',
                     'lo_exp_err', 'hi_exp_err', 'ndata']
        )

        if self.fit_temp_params is None:
            self.fit_temp_params = df_fits
        else:
            self.fit_temp_params = pd.concat(
                [self.fit_temp_params, df_fits],
                axis=0,
                ignore_index=True
            )

        return df_fits
    
    def fit_temp_simple(self, plot_resid=True):
        """
        Fit temperature kT ∝ t^x using scipy.optimize.curve_fit from data in `.params`.

        Args:
            plot_resid (bool, optional): If True, show residual plot. Defaults to True.

        Returns:
            pandas.DataFrame: same columns as `fit_temp_mcmc` with fit='curve_fit'.
        """
        df = self.manager.params.copy()

        # ─────────────────────────────────────────────────────────────
        #  keep only rows whose model string contains “bremss”
        # ─────────────────────────────────────────────────────────────
        m_bremss = df['model'].str.contains('bremss', case=False, na=False)
        df       = df.loc[m_bremss]

        if df.empty:
            raise ValueError("manager.params contains no rows with a 'bremss' model.")

        models = [df['model'].unique()[0]]
        
        records = []

        for mod in models:
            sub = df[df['model'] == mod]

            # require positive times & temperatures
            m_valid = (sub['time_since_explosion'] > 0) & (sub['bremss_kT'] > 0)
            sub     = sub.loc[m_valid]
            
            if len(sub) < 2:
                continue                       # not enough points to fit
            
            t = sub['time_since_explosion']
            t_err = sub['time_since_explosion_err']
            T = sub['bremss_kT']
            lo_T_err = sub['lo_bremss_kT_err']
            hi_T_err = sub['hi_bremss_kT_err']

            norm, norm_err, exp, exp_err = fit_powerlaw_simple(t, T, yerr_lo=lo_T_err, yerr_hi=hi_T_err)
            
            chi2, chi2red = compute_chi2_powerlaw_simple(t, T, yerr_lo=lo_T_err, yerr_hi=hi_T_err, 
                                                         norm=norm, exp=exp, plot_resid=plot_resid)

            records.append({
                'model':     mod,
                'fit':    'curve_fit',
                'norm':       norm,
                'lo_norm_err': norm_err,
                'hi_norm_err': norm_err,
                'exp':        exp,
                'lo_exp_err':   exp_err,
                'hi_exp_err':   exp_err,
                'ndata':     len(sub)
            })

        if not records:
            raise ValueError("No valid ‘bremss’ rows with positive t and kT to fit.")

        df_fits = pd.DataFrame.from_records(
            records,
            columns=['model', 'fit', 'norm', 'lo_norm_err', 'hi_norm_err', 'exp',
                     'lo_exp_err', 'hi_exp_err', 'ndata']
        )

        if self.fit_temp_params is None:
            self.fit_temp_params = df_fits
        else:
            self.fit_temp_params = pd.concat(
                [self.fit_temp_params, df_fits],
                axis=0,
                ignore_index=True
            )

        return df_fits
    
    def calc_density_mcmc(self, distance: float = None, lo_dist_err=None, hi_dist_err=None,
                          z: float = None, lo_z_err=None, hi_z_err=None, mu_e=1.14, mu_ion=1.24, 
                          radius_ratio=1.2, f = 1, H0=70 * u.km / u.s / u.Mpc, freeze_norm=None, freeze_exp=None,
                          nwalkers=500, nsteps=10000, nburn=2000, show_plots=True):
        """
        Calculate CSM density ρ(r) and fit its power-law using MCMC based on 'bremss_norm'.
        
        Formula: 
        rho_CSM = m_p/4 (2*EM*mu_e*mu_ion/V_FS)^(1/2)
        bremss_norm = 3.02e-15/(4 pi d^2) * EM
        V_FS = 4/3 pi f (R_out^3 - R_in^3)

        Args:
            distance (float, optional): distance to source in Mpc.
            z (float, optional): redshift to source.
            mu_e (float, optional): mean molecular weight per electron. Defaults to 1.14. from E. A. Zimmerman et al. 2024 for the environment of SN 2023ixf
            mu_ion (float, optional): mean molecular weight per ion. Defaults to 1.24.
            radius_ratio (float, optional): shock Rout/Rin ratio. Defaults to 1.2.
            f (int, optional): filling factor. Defaults to 1.
            H0 (float, optional): Hubble constant in km/s/Mpc. Defaults to 70 km/s/Mpc.

        Raises:
            RuntimeError: When manager.params do not have these columns: ['bremss_norm', 'lo_bremss_norm_err', 'hi_bremss_norm_err', 'time_since_explosion']
            ValueError: When distance or redshift is not entered

        Returns:
            pd.DataFrame : DataFrame of shock density in g cm-3
        """
        df_norms = None
        try:
            df_norms = self.manager.params.loc[  # rows where ...
                self.manager.params['model'].str.contains('bremss', case=False, na=False), 
                ['bremss_norm', 'lo_bremss_norm_err', 'hi_bremss_norm_err', 'time_since_explosion']
                ]
            df_norms = df_norms.dropna(ignore_index=True)
        except Exception as e:
            raise RuntimeError(e)
        
        
        if distance is None and z is None:
            if self.distance is None:
                raise ValueError("Please enter distance/redshift")
            else:
                distance, lo_dist_err, hi_dist_err = self.distance
        
        distance, lo_dist_err, hi_dist_err = self._distance_and_errors(
                distance, lo_dist_err, hi_dist_err,
                z, lo_z_err, hi_z_err, H0
            )
        
        # Calculate density
            
        distance = distance * MPC_TO_CM
        lo_dist_err = lo_dist_err * MPC_TO_CM
        hi_dist_err = hi_dist_err * MPC_TO_CM
        
        r_in = self.r_shock['r_shock']
        r_in_err = self.r_shock['r_shock_err']
        
        EM_const = 4 * np.pi * distance**2 / 3.02e-15

        EM = np.array(df_norms['bremss_norm'].to_list()) * EM_const
        lo_EM_err = np.array(df_norms['lo_bremss_norm_err'].to_list()) * EM_const
        hi_EM_err = np.array(df_norms['hi_bremss_norm_err'].to_list()) * EM_const
        
        V_FS = 4/3 * np.pi * f * (radius_ratio**3 - 1) * r_in**3
        
        densities = (m_p.to(u.g).value)/4 * (2 * EM * mu_e * mu_ion / V_FS)**(1/2)
        
        rel_D_lo = 2 * lo_dist_err / distance
        rel_D_hi = 2 * hi_dist_err / distance
        lo_EM_err = np.hypot(lo_EM_err, EM * rel_D_lo)
        hi_EM_err = np.hypot(hi_EM_err, EM * rel_D_hi)
        
        rel_EM_lo  = 0.5 * lo_EM_err  / EM
        rel_EM_hi  = 0.5 * hi_EM_err  / EM
        rel_r      = 1.5 * r_in_err   / r_in

        lo_rho_err = densities * np.sqrt(rel_EM_lo**2 + rel_r**2)
        hi_rho_err = densities * np.sqrt(rel_EM_hi**2 + rel_r**2)
        
            
        norm, lo_norm_err, hi_norm_err, exp, lo_exp_err, hi_exp_err, chain = fit_powerlaw_asymmetric(
                                                                                    r_in, densities, lo_rho_err, hi_rho_err,
                                                                                    xerr_lo=r_in_err, xerr_hi=r_in_err,
                                                                                    nwalkers=nwalkers, nsteps=nsteps, 
                                                                                    nburn=nburn, show_plots=show_plots,
                                                                                    freeze_norm=freeze_norm, freeze_exp=freeze_exp
                                                                                 )
            
        compute_chi2_powerlaw(
                r_in, densities, lo_rho_err, hi_rho_err,
                norm, lo_norm_err, hi_norm_err,
                exp, lo_exp_err, hi_exp_err,
                xlo=r_in_err, xhi=r_in_err,
                dof=None, plot_resid=show_plots
            )
        
        df_densities = pd.DataFrame({
            'fit': 'mcmc',
            'time_since_explosion': self.times['time_since_explosion'],
            'time_since_explosion_err': self.times['time_since_explosion_err'],
            'rho': densities,
            'lo_rho_err': lo_rho_err,
            'hi_rho_err': hi_rho_err,
        })
        
        df_dens_params = pd.DataFrame({
            'fit':    'mcmc',
            'norm':       norm,
            'lo_norm_err': lo_norm_err,
            'hi_norm_err': hi_norm_err,
            'exp':        exp,
            'lo_exp_err':   lo_exp_err,
            'hi_exp_err':   hi_exp_err,
        }, index=[0])
        
        if self.densities is None:
            self.densities = df_densities
        else:
            self.densities = pd.concat(
                [self.densities, df_densities],
                axis=0,
                ignore_index=True
            )
            
        if self.fit_density_params is None:
            self.fit_density_params = df_dens_params
        else:
            self.fit_density_params = pd.concat(
                [self.fit_density_params, df_dens_params],
                axis=0,
                ignore_index=True
            )
                    
        return df_densities
    
    def calc_density_simple(self, with_error=False, distance: float = None, lo_dist_err=None, hi_dist_err=None,
                            z: float = None, lo_z_err=None, hi_z_err=None, mu_e=1.14, mu_ion=1.24, 
                            radius_ratio=1.2, f = 1, H0=70 * u.km / u.s / u.Mpc, plot_resid=True):
        """
        Calculate CSM density ρ(r) and fit its power-law using scipy.optimize.curve_fit based on 'bremss_norm'.
        
        Formula: 
        rho_CSM = m_p/4 (2*EM*mu_e*mu_ion/V_FS)^(1/2)
        bremss_norm = 3.02e-15/(4 pi d^2) * EM
        V_FS = 4/3 pi f (R_out^3 - R_in^3)

        Args:
            with_error (bool, optional): fit the density with its error or not. Defaults to False.
            distance (float, optional): distance to source in Mpc.
            z (float, optional): redshift to source.
            mu_e (float, optional): mean molecular weight per electron. Defaults to 1.14. from E. A. Zimmerman et al. 2024 for the environment of SN 2023ixf
            mu_ion (float, optional): mean molecular weight per ion. Defaults to 1.24.
            radius_ratio (float, optional): shock Rout/Rin ratio. Defaults to 1.2.
            f (int, optional): filling factor. Defaults to 1.
            H0 (float, optional): Hubble constant in km/s/Mpc. Defaults to 70 km/s/Mpc.

        Raises:
            RuntimeError: When manager.params do not have these columns: ['bremss_norm', 'lo_bremss_norm_err', 'hi_bremss_norm_err', 'time_since_explosion']
            ValueError: When distance or redshift is not entered

        Returns:
            pd.DataFrame : DataFrame of shock density in g cm-3
        """
        df_norms = None
        try:
            df_norms = self.manager.params.loc[  # rows where ...
                self.manager.params['model'].str.contains('bremss', case=False, na=False), 
                ['bremss_norm', 'lo_bremss_norm_err', 'hi_bremss_norm_err', 'time_since_explosion']
                ]
            df_norms = df_norms.dropna(ignore_index=True)
        except Exception as e:
            raise RuntimeError(e)
        
        
        if distance is None and z is None:
            if self.distance is None:
                raise ValueError("Please enter distance/redshift")
            else:
                distance, lo_dist_err, hi_dist_err = self.distance
        
        distance, lo_dist_err, hi_dist_err = self._distance_and_errors(
                distance, lo_dist_err, hi_dist_err,
                z, lo_z_err, hi_z_err, H0
            )
        
        # Calculate density
            
        distance = distance * MPC_TO_CM
        lo_dist_err = lo_dist_err * MPC_TO_CM
        hi_dist_err = hi_dist_err * MPC_TO_CM
        
        r_in = self.r_shock['r_shock']
        r_in_err = self.r_shock['r_shock_err']
        
        EM_const = 4 * np.pi * distance**2 / 3.02e-15

        EM = np.array(df_norms['bremss_norm'].to_list()) * EM_const
        lo_EM_err = np.array(df_norms['lo_bremss_norm_err'].to_list()) * EM_const
        hi_EM_err = np.array(df_norms['hi_bremss_norm_err'].to_list()) * EM_const
        
        V_FS = 4/3 * np.pi * f * (radius_ratio**3 - 1) * r_in**3
        
        densities = (m_p.to(u.g).value)/4 * (2 * EM * mu_e * mu_ion / V_FS)**(1/2)
        
        rel_D_lo = 2 * lo_dist_err / distance
        rel_D_hi = 2 * hi_dist_err / distance
        lo_EM_err = np.hypot(lo_EM_err, EM * rel_D_lo)
        hi_EM_err = np.hypot(hi_EM_err, EM * rel_D_hi)
        
        rel_EM_lo  = 0.5 * lo_EM_err  / EM
        rel_EM_hi  = 0.5 * hi_EM_err  / EM
        rel_r      = 1.5 * r_in_err   / r_in

        lo_rho_err = densities * np.sqrt(rel_EM_lo**2 + rel_r**2)
        hi_rho_err = densities * np.sqrt(rel_EM_hi**2 + rel_r**2)
        
        if with_error:
            norm, norm_err, exp, exp_err = fit_powerlaw_simple(r_in, densities, yerr_lo=lo_rho_err, yerr_hi=hi_rho_err)
                
            chi2, chi2red = compute_chi2_powerlaw_simple(r_in, densities, yerr_lo=lo_rho_err, yerr_hi=hi_rho_err, 
                                                         norm=norm, exp=exp, plot_resid=plot_resid)
        else:
            norm, norm_err, exp, exp_err = fit_powerlaw_simple(r_in, densities)
                
            chi2, chi2red = compute_chi2_powerlaw_simple(r_in, densities, 
                                                         norm=norm, exp=exp, plot_resid=plot_resid)
        
        df_densities = pd.DataFrame({
            'fit': 'curve_fit',
            'time_since_explosion': self.times['time_since_explosion'],
            'time_since_explosion_err': self.times['time_since_explosion_err'],
            'rho': densities,
            'lo_rho_err': lo_rho_err,
            'hi_rho_err': hi_rho_err,
        })
        
        df_dens_params = pd.DataFrame({
            'fit':    'curve_fit',
            'norm':       norm,
            'lo_norm_err': norm_err,
            'hi_norm_err': norm_err,
            'exp':        exp,
            'lo_exp_err':   exp_err,
            'hi_exp_err':   exp_err,
        }, index=[0])
        
        if self.densities is None:
            self.densities = df_densities
        else:
            self.densities = pd.concat(
                [self.densities, df_densities],
                axis=0,
                ignore_index=True
            )
            
        if self.fit_density_params is None:
            self.fit_density_params = df_dens_params
        else:
            self.fit_density_params = pd.concat(
                [self.fit_density_params, df_dens_params],
                axis=0,
                ignore_index=True
            )
                    
        return df_densities
    
    def calc_mass_loss_rate(self, v_wind: float = v_wind.value):
        """
        Calculate mass loss rate in Msolar/yr
        
        Formula:
        M_dot = 4 pi R_sh^2 rho v_wind
        
        Since we know that rho = norm * r^2
        M_dot = 4 pi norm v_wind
        
        Args:
            v_wind (float, optional): stellar wind velocity in km/s. Defaults to 20 km/s.

        Returns:
            Mass-loss rate and its errors in a list: (m_dot, lo_m_dot_err, hi_m_dot_err)
        """
        try:
            dens_params = self.fit_density_params
        except Exception as e:
            raise RuntimeError("Please fit density first! (Hint: use the calc_density_simple() or calc_density_mcmc() function)")

        v_wind = (v_wind * u.km / u.s).to(u.cm / u.s).value
        K = np.array(dens_params['norm'])
        Kerr_lo = np.array(dens_params['lo_norm_err'])
        Kerr_hi = np.array(dens_params['hi_norm_err'])
        
        m_dot = 4 * np.pi * v_wind * K
        lo_m_dot_err = 4 * np.pi * v_wind * Kerr_lo
        hi_m_dot_err = 4 * np.pi * v_wind * Kerr_hi
        
        m_dot = (m_dot * u.g / u.s).to(u.M_sun / u.year).value
        lo_m_dot_err = (lo_m_dot_err * u.g / u.s).to(u.M_sun / u.year).value
        hi_m_dot_err = (hi_m_dot_err * u.g / u.s).to(u.M_sun / u.year).value
        
        df_m_dot = pd.DataFrame({'fit': dens_params['fit'], 'm_dot': m_dot, 'lo_m_dot_err': lo_m_dot_err, 'hi_m_dot_err': hi_m_dot_err})
        self.mass_loss_rate = pd.concat(
                [self.mass_loss_rate, df_m_dot],
                axis=0,
                ignore_index=True
            )
        return self.mass_loss_rate
    
    def plot_lumin(self, mcmc: bool = True, model_color: str = 'red'):
        """
        Plot luminosity light curve with best-fit power-law.

        Args:
            mcmc (bool, optional): If True, use MCMC fit; else curve_fit. Defaults to True.
            model_color (str, optional): Color for the fit line. Defaults to 'red'.

        Returns:
            matplotlib.figure.Figure: Figure with data and fit.
        """
        if self.fit_lumin_params is None:
            raise RuntimeError("Run `fit_lumin()` first.")

        # choose the fit type
        fit_type = 'mcmc' if mcmc else 'curve_fit'
        df_fit = self.fit_lumin_params[self.fit_lumin_params['fit'] == fit_type]
        if df_fit.empty:
            raise ValueError(f"No entries with fit=='{fit_type}' found.")

        # assume only one model row
        row = df_fit.iloc[0]
        model_name = row['model']
        norm = row['norm']
        exp  = row['exp']

        # data for this model
        df_lum = self.manager.lumin.copy()
        sub = df_lum[df_lum['model'] == model_name]
        if sub.empty:
            raise ValueError(f"No luminosity data for model '{model_name}'.")

        # set up figure
        fig, ax = plt.subplots()
        # plot each instrument
        for instr, g in sub.groupby('instrument'):
            t = g['time_since_explosion'].values
            L = g['lumin'].values
            yerr = np.vstack([g['lo_lumin_err'].values,
                              g['hi_lumin_err'].values])
            xerr = g['time_since_explosion_err'].values
            ax.errorbar(t, L, yerr=yerr, xerr=xerr,
                        fmt='o', label=instr)

        # fitted power-law
        t_min, t_max = sub['time_since_explosion'].min(), sub['time_since_explosion'].max()
        t_line = np.logspace(np.log10(t_min), np.log10(t_max), 200)
        L_line = norm * t_line**exp
        ax.plot(t_line, L_line, color=model_color, ls='--',
                label=fr'Fit: $L\propto t^{{{exp:.2f}}}$')

        # cosmetics
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Time since explosion (days)')
        ax.set_ylabel(r'Luminosity (erg s$^{-1}$)')
        ax.set_title(f"Luminosity light curve (fit type: {fit_type})")
        ax.legend()
        ax.grid(True, which='both', ls=':')

        return fig

    def plot_temp(self, mcmc: bool = True, model_color: str = 'red'):
        """
        Plot temperature evolution with best-fit power-law.

        Args:
            mcmc (bool, optional): If True, use MCMC fit; else curve_fit. Defaults to True.
            model_color (str, optional): Color for the fit line. Defaults to 'red'.

        Returns:
            matplotlib.figure.Figure: Figure with data and fit.
        """
        if self.fit_temp_params is None:
            raise RuntimeError("Run `fit_temp()` first to populate self.fit_temp_params.")

        # choose the fit type
        fit_type = 'mcmc' if mcmc else 'curve_fit'
        df_fit = self.fit_temp_params[self.fit_temp_params['fit'] == fit_type]
        if df_fit.empty:
            raise ValueError(f"No entries with fit=='{fit_type}' found.")

        # assume only one model row
        row = df_fit.iloc[0]
        model_name = row['model']
        norm = row['norm']
        exp  = row['exp']

        # data for this model
        df_params = self.manager.params.copy()
        sub = df_params[df_params['model'] == model_name]
        if sub.empty:
            raise ValueError(f"No parameter data for model '{model_name}'.")

        # set up figure
        fig, ax = plt.subplots()
        # plot each instrument
        for instr, g in sub.groupby('instrument'):
            t = g['time_since_explosion'].values
            T = g['bremss_kT'].values
            yerr = np.vstack([g['lo_bremss_kT_err'].values,
                              g['hi_bremss_kT_err'].values])
            xerr = g['time_since_explosion_err'].values
            ax.errorbar(t, T, yerr=yerr, xerr=xerr,
                        fmt='o', label=instr)

        # fitted power-law
        t_min, t_max = sub['time_since_explosion'].min(), sub['time_since_explosion'].max()
        t_line = np.logspace(np.log10(t_min), np.log10(t_max), 200)
        L_line = norm * t_line**exp
        ax.plot(t_line, L_line, color=model_color, ls='--',
                label=fr'Fit: $kT \propto t^{{{exp:.2f}}}$')

        # cosmetics
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Time since explosion (days)')
        ax.set_ylabel(r'Bremsstrahlung $kT$ (keV)')
        ax.set_title(f"Temperature evolution (fit type: {fit_type})")
        ax.grid(True, which='both', ls=':')
        ax.legend()

        return fig
                  
    def plot_density(self, mcmc: bool=True, model_color: str = "red"):
        """
        Plot CSM density profile ρ(r) with best-fit power-law.

        Args:
            mcmc (bool, optional): If True, use MCMC fit; else curve_fit. Defaults to True.
            model_color (str, optional): Color for the fit line. Defaults to 'red'.

        Returns:
            matplotlib.figure.Figure: Figure with data and fit.

        Raises:
            RuntimeError: If densities, r_shock or fit parameters are missing.
            ValueError: If no data for the selected fit type.
        """
        # sanity checks
        if self.densities is None or self.fit_density_params is None or self.r_shock is None:
            raise RuntimeError("Run density calc, r_shock and its fit first.")

        # pick which fit
        fit_type = 'mcmc' if mcmc else 'curve_fit'

        # select only that fit’s points
        df_rho   = self.densities[self.densities['fit'] == fit_type]
        df_shock = self.r_shock     
        df_fit   = self.fit_density_params[self.fit_density_params['fit'] == fit_type]

        if df_rho.empty or df_shock.empty or df_fit.empty:
            raise ValueError(f"No data for fit='{fit_type}'")

        # extract fit parameters
        norm = df_fit["norm"].iloc[0]
        exp  = df_fit["exp"].iloc[0]

        # build figure
        fig, ax = plt.subplots()

        # plot data with asymmetric y-errors
        yerr = np.vstack([df_rho["lo_rho_err"], df_rho["hi_rho_err"]])
        ax.errorbar(
            df_shock["r_shock"], df_rho["rho"],
            yerr=yerr, xerr=df_shock["r_shock_err"],
            fmt="o", label="data"
        )

        # overlay power-law fit
        rmin, rmax = df_shock["r_shock"].min(), df_shock["r_shock"].max()
        r_line   = np.logspace(np.log10(rmin), np.log10(rmax), 200)
        rho_line = norm * r_line**exp
        ax.plot(
            r_line, rho_line,
            color=model_color, ls="--",
            label=fr"Fit: $\rho\propto r^{{{exp:.2f}}}$"
        )

        # log–log formatting
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Shock radius (cm)")
        ax.set_ylabel(r"Density $\rho$ (g cm$^{-3}$)")
        ax.set_title(f"CSM density profile (fit type: {fit_type})")
        ax.grid(True, which="both", ls=":")
        ax.legend()

        return fig