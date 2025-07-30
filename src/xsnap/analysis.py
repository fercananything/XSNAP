from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .spectrum import SpectrumFit, SpectrumManager
import astropy.units as u 
from ._fitting import fit_powerlaw_asymmetric, compute_chi2_powerlaw, fit_mdot_mcmc
from scipy.optimize import curve_fit
import emcee, corner

# Constants
C_LIGHT_KM_S = 299792.458 * u.km / u.s
MPC_TO_CM = (1 * u.Mpc).to(u.cm).value
v_wind = 20 * u.km / u.s
v_shock = 10000 * u.km / u.s
m_p = 1.67262192e-27 * u.kg

class CSMAnalysis:
    """
    Supernova circumstellar medium (CSM) analysis, 
    including fitting :math:`L(t) \\propto t^x` to deriving unshocked CSM densities and 
    mass-loss rates of the supernova progenitor.
    
    .. note::
        As our emission measure (EM) comes from the normalization of the thermal-bremsstrahlung,
        the :py:class:`~xsnap.analysis.CSMAnalysis` class can only be used 
        if the :py:class:`~xsnap.spectrum.SpectrumManager` class contains **exactly one** bremss model.
    
    Attributes
    -----------
        manager : SpectrumManager or None
            Input manager spectra analysis.
        distance : float or None
            Source distance in Mpc.
        r_shock : pandas.DataFrame or None
            Shock radius.
        times : pandas.DataFrame or None
            Times since explosion.
        fit_lumin_params : pandas.DataFrame or None
            Best-fit luminosity parameters wrapped in DataFrame table with 
                columns ``['model', 'norm','lo_norm_err','hi_norm_err', 'exp','lo_exp_err','hi_exp_err','ndata']``.
        fit_temp_params : pandas.DataFrame or None
            Best-fit temperature parameters wrapped in DataFrame table with 
                columns ``['model', 'norm','lo_norm_err','hi_norm_err', 'exp','lo_exp_err','hi_exp_err','ndata']``.
        fit_density_params : pandas.DataFrame or None
            Best-fit density parameters, i.e. the fitted mass-loss rates in :math:`\\rm g \ s^{-1}`, wrapped in 
            DataFrame table with columns ``['mdot', 'lo_mdot_err', 'hi_mdot_err']``
        densities : pandas.DataFrame or None
            Computed unshocked CSM density in :math:`\\rm g \ {cm}^{-3}` wrapped in DataFrame table with columns
            ``['time_since_explosion', 'rho', 'lo_rho_err', 'hi_rho_err']``
        mass_loss_rate : pandas.DataFrame
            Computed mass-loss rates of the supernova progenitor in :math:`\\rm M_{\\bigodot} \ {yr}^{-1}`,
            wrapped in a DataFrame with columns ``['m_dot', 'lo_m_dot_err', 'hi_m_dot_err']``
    """
    def __init__(self, manager=None):
        """
        Initialization of the :py:class:`~xsnap.analysis.CSMAnalysis` class

        Parameters
        ----------
            manager : SpectrumManager, optional
                Pre-populated :py:class:`~xsnap.spectrum.SpectrumManager` instance.
        """
        self.manager = None
        self.distance = None
        self.r_shock = None
        self.times = None
        self.fit_lumin_params = None
        self.fit_temp_params = None
        self.fit_density_params = None
        self.densities = None
        self.mass_loss_rate = pd.DataFrame({'m_dot': [], 'lo_m_dot_err': [], 'hi_m_dot_err': []})
        
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
             z=None, lo_z_err=None, hi_z_err=None, v_shock: float = v_shock.value, H0: float = 70):
        """
        Load manager luminosity & parameter tables and compute shock radius.

        Parameters
        ----------
            manager : SpectrumManager
                Must contain one ``'bremss'`` model in :py:attr:`~xsnap.spectrum.SpectrumManager.lumin` & :py:attr:`~xsnap.spectrum.SpectrumManager.params`.
            distance : float, optional
                Distance in Mpc. Required if ``z`` is ``None``.
            z : float, optional
                Redshift to compute distance. Required if `distance` is ``None``. Defaults to ``None``.
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
            v_shock : float, optional
                Shock velocity in km/s. Defaults to ``10000``.

        Returns
        --------
            *self* : CSMAnalysis
                The same :py:class:`~xsnap.analysis.CSMAnalysis` with 
                :py:attr:`~xsnap.analysis.CSMAnalysis.times` and :py:attr:`~xsnap.analysis.CSMAnalysis.r_shock` populated.

        Raises
        --------
            RuntimeError
                If more than one bremss model is present.
            ValueError
                If required tables/columns are missing.
                
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
        H0 = H0 * u.km / u.s / u.Mpc
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
        self.mass_loss_rate = pd.DataFrame({'m_dot': [], 'lo_m_dot_err': [], 'hi_m_dot_err': []})

    def fit_lumin(self, nwalkers=500, nsteps=10000, nburn=2000, show_plots=True):
        """
        Fit power-law luminosity :math:`L(t) \\propto t^x` using Markov chain Monte-Carlo (MCMC).

        Parameters
        ----------
            nwalkers : int, optional
                Number of MCMC walkers. Defaults to ``500``.
            nsteps : int, optional
                Number of steps per walker. Defaults to ``10000``.
            nburn : int, optional
                Number of burn-in steps. Defaults to ``2000``.
            show_plots : bool, optional
                If ``True``, display diagnostic plots. Defaults to ``True``

        Returns
        -------
            Best-fit luminosity parameters : pandas.DataFrame
                Best-fit luminosity parameters wrapped in DataFrame table with 
                columns ``['model', 'norm','lo_norm_err','hi_norm_err', 'exp','lo_exp_err','hi_exp_err','ndata']``.
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
            columns=['model',   'norm', 'lo_norm_err', 'hi_norm_err', 'exp',
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

    
    def fit_temp(self, nwalkers=500, nsteps=10000, nburn=2000, show_plots=True):
        """
        Fit power-law temperature :math:`T(t) \\propto t^x` using Markov chain Monte-Carlo (MCMC).

        Parameters
        ----------
            nwalkers : int, optional
                Number of MCMC walkers. Defaults to ``500``.
            nsteps : int, optional
                Number of steps per walker. Defaults to ``10000``.
            nburn : int, optional
                Number of burn-in steps. Defaults to ``2000``.
            show_plots : bool, optional
                If ``True``, display diagnostic plots. Defaults to ``True``

        Returns
        --------
            Best-fit temperature parameters : pandas.DataFrame
                Best-fit temperature parameters wrapped in DataFrame table with 
                columns ``['model', 'norm','lo_norm_err','hi_norm_err', 'exp','lo_exp_err','hi_exp_err','ndata']``.
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
            columns=['model',   'norm', 'lo_norm_err', 'hi_norm_err', 'exp',
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
    

    def calc_density(self, distance=None, lo_dist_err=None, hi_dist_err=None,
                    z=None, lo_z_err=None, hi_z_err=None, mu_e=1.14, mu_ion=1.24,
                    radius_ratio=1.2, f=1, v_wind=v_wind.value, H0: float =70,
                    nwalkers=500, nsteps=10000, nburn=2000, show_plots=True):
        """
        Calculate unshocked CSM density :math:`\\rho_{\\mathrm{CSM}}(r)` and fit using Markov chain Monte-Carlo (MCMC) 
        based on ``'bremss_norm'`` in :py:attr:`~xsnap.spectrum.SpectrumManager.params`.
        
        In calculating the density, we use the model of `D. Brethauer et al. (2022) <https://doi.org/10.3847/1538-4357/ac8b14>`_,
        where the unshocked CSM density is expressed by: 
        
        .. math::
            \\rho_{\\mathrm{CSM}}(r) = \\frac{m_p}{4} \left(\\frac{2 \\times \\mathrm{EM}(r)\\mu_e \\mu_I}{V_{\\mathrm{FS}}(r)}\\right)^{1/2}
        
        where we parameterized:
        
        .. math::
            \\mathrm{{bremss} \ {norm}} = \\frac{3.02 \\cdot 10^{-15}}{4 \\pi d^2} \\times \\mathrm{EM} 
            
            \\mathrm{V}_{\\mathrm{FS}} = \\frac{4 \\pi}{3} f \\left(R_{\\mathrm{out}}^3 - R_{\\mathrm{in}}^3\\right)
        
        From there, we fit :math:`\\rho_{\\mathrm{CSM}}(r)` and :math:`\\dot{M}` with:
        
        .. math::
            \\rho_{\\mathrm{CSM}} = \\frac{\\dot{M}}{4 \\pi r_{\\mathrm{shock}}^2 v_{\\mathrm{wind}}}

        Parameters
        ----------
            distance : float, optional
                distance to source in Mpc.
            z : float, optional
                redshift to source.
            mu_e : float, optional
                mean molecular weight per electron. Defaults to ``1.14`` (solar value).
            mu_ion : float, optional
                mean molecular weight per ion. Defaults to ``1.24`` (solar value).
            radius_ratio : float, optional
                shock Rout/Rin ratio. Defaults to ``1.2``.
            f : int, optional
                filling factor. Defaults to ``1``.
            v_wind : float, optional
                wind velocity in km/s. Defaults to ``20``.
            H0 : float, optional
                Hubble constant in km/s/Mpc. Defaults to ``70``.
            nwalkers : int, optional
                Number of MCMC walkers. Defaults to ``500``.
            nsteps : int, optional
                Number of steps per walker. Defaults to ``10000``.
            nburn : int, optional
                Number of burn-in steps. Defaults to ``2000``.
            show_plots : bool, optional
                If ``True``, display diagnostic plots. Defaults to ``True``.
            fit : bool, optional
                If ``True``, will fit density and get mass-loss rate in :math:`\\rm g \ s^{-1}`. Defaults to ``True``

        Returns
        ---------
            Density : pandas.DataFrame 
                Computed unshocked CSM density in :math:`\\rm g \ {cm}^{-3}` wrapped in DataFrame table with columns
                ``['time_since_explosion', 'rho', 'lo_rho_err', 'hi_rho_err']``
                
        Raises
        ---------
            RuntimeError
                When :py:attr:`~xsnap.spectrum.SpectrumManager.params` do not have these columns: 
                ``['bremss_norm', 'lo_bremss_norm_err', 'hi_bremss_norm_err', 'time_since_explosion']``
            ValueError
                When distance or redshift is not entered
                
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
        H0 = H0 * u.km / u.s / u.Mpc
        v_wind = v_wind * u.km / u.s
        df_norms = None
        try:
            df_norms = self.manager.params.loc[ 
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
        D  = distance * MPC_TO_CM
        sD_lo = lo_dist_err * MPC_TO_CM
        sD_hi = hi_dist_err * MPC_TO_CM


        r_in     = self.r_shock["r_shock"]
        r_in_err = self.r_shock["r_shock_err"]


        EM_const  = 4*np.pi*D**2/3.02e-15
        EM        = df_norms["bremss_norm"].values        * EM_const
        sEM_lo    = df_norms["lo_bremss_norm_err"].values * EM_const
        sEM_hi    = df_norms["hi_bremss_norm_err"].values * EM_const

        rel_D_lo  = 2*sD_lo/D
        rel_D_hi  = 2*sD_hi/D
        sEM_lo    = np.hypot(sEM_lo, EM*rel_D_lo)
        sEM_hi    = np.hypot(sEM_hi, EM*rel_D_hi)


        V_FS = 4/3 * np.pi * f * (radius_ratio**3 - 1) * r_in**3
        rho = (m_p.to(u.g).value/4) * np.sqrt( 2*EM*mu_e*mu_ion / V_FS )

        rel_EM_lo = 0.5*sEM_lo/EM
        rel_EM_hi = 0.5*sEM_hi/EM
        rel_r     = 1.5*r_in_err/r_in

        lo_rho_err = rho * np.sqrt(rel_EM_lo**2 + rel_r**2)
        hi_rho_err = rho * np.sqrt(rel_EM_hi**2 + rel_r**2)

        if len(np.array(rho)) >= 2:
            mdot, mdot_lo, mdot_hi = fit_mdot_mcmc(
                r_in, rho, lo_rho_err, hi_rho_err,
                v_wind=v_wind,
                nwalkers=nwalkers, nsteps=nsteps, nburn=nburn,
                show_corner=show_plots
            )


            if show_plots:
                compute_chi2_powerlaw(      # your existing routine
                    r_in, rho, lo_rho_err, hi_rho_err,
                    mdot/(4*np.pi*v_wind.to(u.cm/u.s).value), 0, 0,   # norm (not used)
                    -2, 0, 0,                                       # slope fixed −2
                    xlo=r_in_err, xhi=r_in_err, plot_resid=True
                )

            df_par = pd.DataFrame({
                "mdot": mdot,
                "lo_mdot_err": mdot_lo, "hi_mdot_err": mdot_hi
            }, index=[0])
            
            self.fit_density_params = (df_par if self.fit_density_params is None
                                else pd.concat([self.fit_density_params, df_par],
                                                ignore_index=True))

        df_out = pd.DataFrame({
            "time_since_explosion": df_norms["time_since_explosion"],
            "rho": rho, "lo_rho_err": lo_rho_err, "hi_rho_err": hi_rho_err
        })

        self.densities = (df_out if self.densities is None
                        else pd.concat([self.densities, df_out], ignore_index=True))
        

        return df_out
    
    
    def get_mdot(self):
        """
        Get mass loss rate in :math:`\\rm M_{\\bigodot} \ {yr}^{-1}`
            
        Returns
        ---------
            Mass-loss rate and its errors in :math:`\\rm M_{\\bigodot} \ {yr}^{-1}`, wrapped in
            a DataFrame with columns ``['m_dot', 'lo_m_dot_err', 'hi_m_dot_err']``
        """
        try:
            dens_params = self.fit_density_params
            m_dot = np.array(dens_params['mdot'])
            lo_m_dot_err = np.array(dens_params['lo_mdot_err'])
            hi_m_dot_err = np.array(dens_params['hi_mdot_err'])
        except Exception as e:
            raise RuntimeError("Please fit density first! (Hint: use the calc_density()function)")
        
        m_dot = (m_dot * u.g / u.s).to(u.M_sun / u.year).value
        lo_m_dot_err = (lo_m_dot_err * u.g / u.s).to(u.M_sun / u.year).value
        hi_m_dot_err = (hi_m_dot_err * u.g / u.s).to(u.M_sun / u.year).value
        
        df_m_dot = pd.DataFrame({'m_dot': m_dot, 'lo_m_dot_err': lo_m_dot_err, 'hi_m_dot_err': hi_m_dot_err})
        self.mass_loss_rate = pd.concat(
                [self.mass_loss_rate, df_m_dot],
                axis=0,
                ignore_index=True
            )
        return self.mass_loss_rate
    
    def plot_lumin(self, model_color: str = 'red'):
        """
        Plot luminosity light curve with best-fit power-law.

        Parameters
        -----------
            model_color : str, optional
                Color for the fit line. Defaults to ``'red'``.

        Returns
        ---------
            Plot : matplotlib.figure.Figure
                The fitted luminosity light curve plot with data and fit.
                
                Data are grouped and labeled by instruments.
        """
        if self.fit_lumin_params is None:
            raise RuntimeError("Run `fit_lumin()` first.")

        # choose the fit type
        df_fit = self.fit_lumin_params
        if df_fit.empty:
            raise ValueError(f"Please fit the luminosity with fit_lumin() first.")

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
        ax.set_title(f"Luminosity light curve")
        ax.legend()
        ax.grid(True, which='both', ls=':')

        return fig

    def plot_temp(self, model_color: str = 'red'):
        """
        Plot temperature evolution with best-fit power-law.

        Parameters
        -----------
            model_color : str, optional
                Color for the fit line. Defaults to ``'red'``.

        Returns
        ---------
            Plot : matplotlib.figure.Figure
                The fitted temperature evolution plot with data and fit.
                
                Data are grouped and labeled by instruments.
        """
        if self.fit_temp_params is None:
            raise RuntimeError("Run `fit_temp()` first to populate self.fit_temp_params.")

        # choose the fit type
        df_fit = self.fit_temp_params
        if df_fit.empty:
            raise ValueError(f"Please fit the temperature with fit_temp() first.")

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
        ax.set_title(f"Temperature evolution")
        ax.grid(True, which='both', ls=':')
        ax.legend()

        return fig
    
    def plot_density(self, model_color: str = "red", v_wind: float = 20):
        """
        Plot the unshocked CSM density profile with its best-fit mass-loss rate.

        Parameters
        -----------
            model_color : str, optional
                Color for the fit line. Defaults to ``'red'``.
            v_wind      : float
                Wind velocity (needed as the best-fit parameter is in terms of the mass-loss rate).

        Returns
        ---------
            Plot : matplotlib.figure.Figure
                The density profile plot with data and fit
        """
        v_wind = v_wind * u.km / u.s

        # ─────────── basic checks ──────────────────────────────────────────
        if self.densities is None or self.fit_density_params is None or self.r_shock is None:
            raise RuntimeError("Run density calculation, r_shock and its fit first.")


        df_rho   = self.densities
        df_shock = self.r_shock
        df_fit   = self.fit_density_params

        if df_rho.empty or df_fit.empty or df_shock.empty:
            raise ValueError(f"Please fit density with calc_density() first.")

        # ─────────── decide which parametrisation the fit used ─────────────
        if "mdot" in df_fit.columns:               # ← new Ṁ fit
            mdot   = df_fit["mdot"].iloc[0]                 # [Msol yr⁻¹]
            v_cms  = v_wind.to(u.cm/u.s).value
            def rho_model(r):                       # r in cm
                return mdot / (4*np.pi*r**2*v_cms)
            fit_label = fr"$\rho \,\propto\, r^{{-2}}$"     # slope fixed −2

        elif {"norm", "exp"}.issubset(df_fit.columns):      # ← legacy norm/exp
            norm = df_fit["norm"].iloc[0]
            exp  = df_fit["exp"].iloc[0]
            def rho_model(r):
                return norm * r**exp
            fit_label = fr"$\rho \propto r^{{{exp:.2f}}}$"

        else:
            raise RuntimeError("Fit table has neither (norm,exp) nor mdot column.")

        # ─────────── build the figure ──────────────────────────────────────
        fig, ax = plt.subplots()

        # data with asymmetric errors
        yerr = np.vstack([df_rho["lo_rho_err"], df_rho["hi_rho_err"]])
        ax.errorbar(
            df_shock["r_shock"], df_rho["rho"],
            yerr=yerr, xerr=df_shock["r_shock_err"],
            fmt="o", label="data"
        )

        # fitted curve
        rmin, rmax = df_shock["r_shock"].min(), df_shock["r_shock"].max()
        r_line   = np.logspace(np.log10(rmin), np.log10(rmax), 200)
        ax.plot(r_line, rho_model(r_line), ls="--", color=model_color, label=fit_label)

        # aesthetics
        ax.set_xscale("log");  ax.set_yscale("log")
        ax.set_xlabel("Shock radius (cm)")
        ax.set_ylabel(r"Density $\rho\;(\mathrm{g\,cm^{-3}})$")
        ax.set_title("CSM density profile")
        ax.grid(True, which="both", ls=":")
        ax.legend()

        return fig