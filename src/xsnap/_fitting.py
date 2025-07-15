import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Helper to fit asymmetric powerlaw

def fit_powerlaw_asymmetric(
    x, y, yerr_lo, yerr_hi,
    xerr_lo=None, xerr_hi=None,
    nwalkers=200, nsteps=5000, nburn=1000,
    show_plots=True,
    freeze_norm=None,   # if not None, fix A=freeze_norm
    freeze_exp =None    # if not None, fix k=freeze_exp
):
    """Helper function that fit a powerlaw y = A * x^k model with asymmetric errors with MCMC

    Args:
        x (array_like): x-variable values, shape (N,).
        y (array_like): y-variable values, shape (N,).
        yerr_lo (array_like): Lower uncertainty on y, shape (N,).
        yerr_hi (array_like): Upper uncertainty on y, shape (N,).
        xerr_lo (array_like, optional): Lower uncertainty on x, shape (N,). If None, x errors are ignored. Defaults to None.
        xerr_hi (array_like, optional): Upper uncertainty on x, shape (N,). If None, x errors are ignored. Defaults to None.
        nwalkers (int, optional): Number of walkers in the MCMC ensemble sampler. Defaults to 200.
        nsteps (int, optional): Total number of MCMC steps to run for each walker. Defaults to 5000.
        nburn (int, optional): Number of initial steps per walker to discard as burn-in. Defaults to 1000.
        show_plots (bool, optional): If True, displays diagnostic plots (trace plots, corner plot, fit overlay). Set to False to suppress plotting. Defaults to True.

    Raises:
        ValueError: When the lower or upper uncertainties of y are zero-values

    Returns:
        tuple: 
            A_m (float): Median value of the normalization A.
            A_minus_err (float): Absolute lower uncertainty on A.
            A_plus_err (float): Absolute upper uncertainty on A.
            k_m (float): Median value of the exponent k.
            k_minus_err (float): Absolute lower uncertainty on k.
            k_plus_err (float): Absolute upper uncertainty on k.
            chain (numpy.ndarray): Flattened MCMC chain of shape (nwalkers*(nsteps-nburn), ndim) containing [A, k] samples.
    """
    # you can’t freeze both
    if freeze_norm is not None and freeze_exp is not None:
        raise ValueError("Can only freeze one of norm or exponent.")

    # require non-zero y-errors
    if (np.array(yerr_hi)==0).any() or (np.array(yerr_lo)==0).any():
        raise ValueError("Please supply non-zero y-errors.")

    # propagate x-errors into y-errors (as before)...
    if xerr_lo is not None and xerr_hi is not None:
        u, v = np.log10(x), np.log10(y)
        A_mat = np.vstack([u, np.ones_like(u)]).T
        k0, logA0 = np.linalg.lstsq(A_mat, v, rcond=None)[0]
        yerr_lo = np.hypot(yerr_lo, abs((10**logA0)*(x**k0)*(k0*xerr_lo/x)))
        yerr_hi = np.hypot(yerr_hi, abs((10**logA0)*(x**k0)*(k0*xerr_hi/x)))

    # define the model functions
    def model_both(theta, xx):
        logA, k = theta
        return (10**logA)*xx**k
    def model_k_only(theta, xx):
        k, = theta
        return freeze_norm * xx**k
    def model_A_only(theta, xx):
        logA, = theta
        return (10**logA) * xx**freeze_exp

    # priors
    def log_prior_both(theta):
        logA, k = theta
        return 0.0 if (-10<logA<10 and -5<k<5) else -np.inf
    def log_prior_k(theta):
        k, = theta
        return 0.0 if -5<k<5 else -np.inf
    def log_prior_A(theta):
        logA, = theta
        return 0.0 if -10<logA<10 else -np.inf

    # likelihood (identical structure for all three)
    def log_likelihood(theta, xx, yy, slo, shi, model_fn):
        ypred = model_fn(theta, xx)
        resid = yy - ypred
        # pick correct sigma
        chi2_lo = (resid/shi)**2 + np.log(2*np.pi*shi**2)
        chi2_hi = (resid/slo)**2 + np.log(2*np.pi*slo**2)
        chi2    = np.where(resid<0, chi2_lo, chi2_hi)
        return -0.5*np.sum(chi2)

    # posterior wrappers
    if freeze_norm is not None:
        ndim = 1
        model_fn = model_k_only
        log_prior_fn = log_prior_k
    elif freeze_exp is not None:
        ndim = 1
        model_fn = model_A_only
        log_prior_fn = log_prior_A
    else:
        ndim = 2
        model_fn = model_both
        log_prior_fn = log_prior_both

    def log_post(theta, xx, yy, slo, shi):
        lp = log_prior_fn(theta)
        if not np.isfinite(lp): return -np.inf
        return lp + log_likelihood(theta, xx, yy, slo, shi, model_fn)

    # initialize walkers
    u, v = np.log10(x), np.log10(y)
    A_mat = np.vstack([u, np.ones_like(u)]).T
    k0, logA0 = np.linalg.lstsq(A_mat, v, rcond=None)[0]

    if ndim==2:
        p0 = np.column_stack([
            np.random.normal(logA0, 0.1, size=nwalkers),
            np.random.normal(   k0, 0.1, size=nwalkers),
        ])
    elif freeze_norm is not None:  # sample k only
        p0 = np.random.normal(k0, 0.1, size=(nwalkers,1))
    else:                         # sample logA only
        p0 = np.random.normal(logA0,0.1, size=(nwalkers,1))

    # run sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_post,
                                    args=(x, y, yerr_lo, yerr_hi))
    sampler.run_mcmc(p0, nsteps, progress=True)
    chain = sampler.get_chain(discard=nburn, flat=True)

    # summarize
    if ndim==2:
        logA_m, k_m = np.median(chain, axis=0)
        logA_lo, logA_hi = np.percentile(chain[:,0],[16,84]) - logA_m
        k_lo,    k_hi    = np.percentile(chain[:,1],[16,84]) - k_m
    elif freeze_norm is not None:
        k_m = np.median(chain[:,0])
        k_lo, k_hi = np.percentile(chain[:,0],[16,84]) - k_m
        logA_m, logA_lo, logA_hi = np.log10(freeze_norm), 0.0, 0.0
    else:
        logA_m = np.median(chain[:,0])
        logA_lo, logA_hi = np.percentile(chain[:,0],[16,84]) - logA_m
        k_m, k_lo, k_hi = freeze_exp, 0.0, 0.0

    A_m = 10**logA_m
    A_minus = A_m - 10**(logA_m-logA_lo)
    A_plus  = 10**(logA_m+logA_hi) - A_m

    # optional plots (unchanged)
    if show_plots:
        # === Posterior plot ===
        if ndim == 2:
            # both logA and k were free, so make a corner plot
            fig_corner = corner.corner(
                chain,
                labels=["logA","k"],
                truths=[logA_m, k_m],
                show_titles=True
            )
            fig_corner.suptitle("Posterior on (logA, k)", y=1.02)
        else:
            # only one parameter was free—just histogram it
            fig_hist, ax = plt.subplots()
            ax.hist(chain[:,0], bins=50, histtype='step')
            val = k_m if freeze_norm is not None else logA_m
            label = 'k' if freeze_norm is not None else 'logA'
            ax.axvline(val, color='C1', linestyle='--')
            ax.set_xlabel(label)
            ax.set_title(f"Posterior of {label}")
        # === Fit overlay ===
        plt.figure()
        if xerr_lo is not None and xerr_hi is not None:
            plt.errorbar(x, y, yerr=[yerr_lo,yerr_hi], xerr=[xerr_lo,xerr_hi],
                         fmt='o', capsize=3)
        else:
            plt.errorbar(x, y, yerr=[yerr_lo,yerr_hi], fmt='o', capsize=3)
        xx = np.logspace(np.log10(x.min()), np.log10(x.max()), 200)
        plt.plot(xx, A_m*xx**k_m, 'C1-', lw=2,
                 label=fr"$y \propto x^{{{k_m:.2f}}}$")
        plt.xscale('log'); plt.yscale('log')
        plt.xlabel('x'); plt.ylabel('y')
        plt.legend()
        plt.title('Data + Power‐law fit')
        plt.show()

    return A_m, abs(A_minus), abs(A_plus), k_m, abs(k_lo), k_hi, chain


def predict_with_errors(
    A_m, A_err_lo, A_err_hi,
    k_m, k_err_lo, k_err_hi,
    x_m, x_err_lo = 0, x_err_hi = 0
):
    """
    Helper function to propagate asymmetric errors for y = A * x^k, including error in x,
    using a Taylor‐expansion + quadrature approach.

    Args:
        A_m (float): Best-fit normalization A.
        A_err_lo (float): Lower uncertainty on A (A_m − A_minus).
        A_err_hi (float): Upper uncertainty on A (A_plus − A_m).
        k_m (float): Best-fit exponent k.
        k_err_lo (float): Lower uncertainty on k.
        k_err_hi (float): Upper uncertainty on k.
        x_m (float): Best-fit x value.
        x_err_lo (float): Lower uncertainty on x.
        x_err_hi (float): Upper uncertainty on x.

    Returns:
        tuple:
            y_m (float): Nominal value of y = A_m * x_m**k_m.
            y_err_lo (float): Lower uncertainty on y.
            y_err_hi (float): Upper uncertainty on y.
    """
    # Nominal value
    y_m = A_m * x_m**k_m

    # Partial contributions
    # from A
    delta_A_plus  = x_m**k_m * A_err_hi
    delta_A_minus = x_m**k_m * A_err_lo
    # from k
    delta_k_plus  = A_m * x_m**k_m * np.log(x_m) * k_err_hi
    delta_k_minus = A_m * x_m**k_m * np.log(x_m) * k_err_lo
    # from x
    delta_x_plus  = A_m * k_m * x_m**(k_m - 1) * x_err_hi
    delta_x_minus = A_m * k_m * x_m**(k_m - 1) * x_err_lo

    # Combine all in quadrature
    y_err_hi = np.sqrt(delta_A_plus**2  + delta_k_plus**2  + delta_x_plus**2)
    y_err_lo = np.sqrt(delta_A_minus**2 + delta_k_minus**2 + delta_x_minus**2)

    return y_m, y_err_lo, y_err_hi

def compute_chi2_powerlaw(
    x, y, ylo, yhi,
    A, A_err_lo, A_err_hi,
    k, k_err_lo, k_err_hi,
    xlo=None, xhi=None,
    dof=None, plot_resid=True
):
    """
    Helper function to compute χ² and reduced χ² for y = A * x^k. 
    This function to show the user how fit their fitting is.

    Args:
        x (array_like): x-variable values, shape (N,).
        y (array_like): y-variable values, shape (N,).
        ylo (array_like): Lower uncertainty on y, shape (N,).
        yhi (array_like): Upper uncertainty on y, shape (N,).
        A (float); Model normalization.
        A_err_lo (float): Lower uncertainty on A
        A_err_hi (float): Upper uncertainty on A
        k (float): Model exponent.
        k_err_lo (float): Lower uncertainty on k
        k_err_hi (float): Upper uncertainty on k
        xlo (array_like): Lower uncertainty on x, shape (N,).
        xhi (array_like): Upper uncertainty on x, shape (N,).
        dof (int or None): Degrees of freedom. If None, defaults to N − 2.
        plot_resid (bool): If True, show a residual plot.

    Returns:
        tuple:
            chi2 (float): Total χ².
            chi2_red (float): Reduced χ² = χ²/dof.
    """
    # A_err_lo, A_err_hi, k_err_lo, k_err_hi = 0,0,0,0
    # Propagate A, k (and x, if given) → y_errors
    if (xlo is not None) and (xhi is not None):
        y_m, y_model_lo, y_model_hi = predict_with_errors(
            A_m=A,      A_err_lo=A_err_lo,  A_err_hi=A_err_hi,
            k_m=k,      k_err_lo=k_err_lo,  k_err_hi=k_err_hi,
            x_m=x,      x_err_lo=xlo,       x_err_hi=xhi
        )
        xerr = (xlo, xhi)
    else:
        # no x‐errors: pass zeros
        zeros = np.zeros_like(x)
        y_m, y_model_lo, y_model_hi = predict_with_errors(
            A_m=A,      A_err_lo=A_err_lo,  A_err_hi=A_err_hi,
            k_m=k,      k_err_lo=k_err_lo,  k_err_hi=k_err_hi,
            x_m=x,      x_err_lo=zeros,     x_err_hi=zeros
        )
        xerr = None

    # Compute residuals
    resid = y - y_m

    # Pick the measurement error branch
    sig_y = np.where(resid < 0, yhi, ylo)
    # Pick the model‐param error branch
    sig_model = np.where(resid < 0, y_model_hi, y_model_lo)

    # Total σ = sqrt(σ_meas^2 + σ_model^2)
    sigma = np.hypot(sig_y, sig_model)

    # χ²
    chi2 = np.sum((resid / sigma)**2)

    # Degrees of freedom
    N = x.size
    if dof is None:
        dof = N - 2   # A, k

    chi2_red = chi2 / dof

    # Plotting
    if plot_resid:
        plt.errorbar(
            x, resid,
            yerr=sigma,
            xerr=xerr,
            fmt='o', capsize=3, color='k'
        )
        plt.axhline(0, color='C1', linestyle='--', label='zero residual')
        plt.xscale('log')
        plt.xlabel('x')
        plt.ylabel('Residual = data - model')
        plt.title(f'Residual plot, red χ² = {chi2_red:.2f}')
        plt.legend()
        plt.show()

    return chi2, chi2_red

def fit_powerlaw_simple(x, y, yerr=None, yerr_lo=None, yerr_hi=None, p0=None):
    """Helper function that fit a powerlaw y = norm * x^exp model with scipy.optimize.curve_fit

    Args:
        x (array_like): x-variable values, shape (N,).
        y (array_like): y-variable values, shape (N,).
        yerr (array_like, optional): Symmetric uncertainties on y, shape (N,). If None and both `yerr_lo` and `yerr_hi` are provided, symmetric errors are computed as 0.5*(yerr_lo + yerr_hi). If still None, residuals are weighted by the model prediction. Defaults to None.
        yerr_lo (array_like, optional): Lower uncertainties on y, shape (N,). Used to compute `yerr` if it is None. Defaults to None.
        yerr_hi (array_like, optional): Upper uncertainties on y, shape (N,). Used to compute `yerr` if it is None. Defaults to None.
        p0 (array_like, optional): Initial guess for the parameters (length N). If None, then the initial values will all be 1 (Source: scipy.optimize.curve_fit)
        
    Returns:
        tuple: 
            norm (float): Best-fit normalization
            norm_err (float): Normalization uncertainties
            exp (float): Best-fit exponent
            exp_err (float): Exponent uncertainties
    """
    def model(logx, norm, exp):
        return norm + exp * logx
    
    if yerr is None and yerr_lo is not None and yerr_hi is not None:
        yerr = 0.5*(yerr_lo + yerr_hi)
        
    logx = np.log10(x)
    logy = np.log10(y)
    logyerr = yerr/(y * np.log(10)) if yerr is not None else 0
    
    if (np.array(logyerr) == 0).any():
        popt, pcov = curve_fit(
            model,
            logx, logy,
            p0=p0
        )
    else:
        popt, pcov = curve_fit(
            model,
            logx, logy,
            sigma=logyerr,
            absolute_sigma=True,
            p0=p0
        )
    norm, exp = popt
    norm_err, exp_err = np.sqrt(np.diag(pcov))
    
    norm = 10**(norm)
    norm_err = norm * np.log(10) * norm_err
    
    return norm, norm_err, exp, exp_err


def compute_chi2_powerlaw_simple(x, y, norm, exp,  yerr=None, yerr_lo=None, yerr_hi=None, dof=None, plot_resid=True):
    """Helper function that computes χ² and reduced χ² for model y = norm * x^exp from curve_fit.
    This function to show the user how fit their fitting is.

    Args:
        x (array_like): x-variable values, shape (N,).
        y (array_like): y-variable values, shape (N,).
        norm (float): Model normalization
        exp (float): Model exponent
        yerr (array_like, optional): Symmetric uncertainties on y, shape (N,). If None and both `yerr_lo` and `yerr_hi` are provided, symmetric errors are computed as 0.5*(yerr_lo + yerr_hi). If still None, residuals are weighted by the model prediction. Defaults to None.
        yerr_lo (array_like, optional): Lower uncertainties on y, shape (N,). Used to compute `yerr` if it is None. Defaults to None.
        yerr_hi (array_like, optional): Upper uncertainties on y, shape (N,). Used to compute `yerr` if it is None. Defaults to None.
        dof (int, optional): Degrees of freedom; if None, defaults to N − 2. Defaults to None.
        plot_resid (bool, optional): If True, displays a residual plot with error bars. Defaults to True.

    Returns:
        tuple:
            chi2 (float): Total χ².
            chi2_red (float): Reduced χ² = χ²/dof.
    """
    if yerr is None and yerr_lo is not None and yerr_hi is not None:
        yerr = 0.5*(yerr_lo + yerr_hi)
        
    # Model prediction and residuals
    y_pred = norm * x**exp
    resid = y - y_pred

    # χ² calculation
    chi2 = np.sum((resid / yerr)**2) if yerr is not None else np.sum((resid / y_pred)**2)
    if dof is None:
        dof = x.size - 2
    chi2_red = chi2 / dof

    # Optional residual plot
    if plot_resid:
        plt.errorbar(
            x, resid,
            yerr=yerr,
            fmt='o', capsize=3, color='k'
        )
        plt.axhline(0, color='C1', linestyle='--', label='zero residual')
        plt.xscale('log')
        plt.xlabel('x')
        plt.ylabel('Residual = data - model')
        plt.title(f'Residual plot, red χ² = {chi2_red:.2f}')
        plt.legend()
        plt.show()
        
    return chi2, chi2_red