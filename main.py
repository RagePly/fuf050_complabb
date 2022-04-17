import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.constants as const
import scipy.optimize as optim

def _electron_speed(E): return math.sqrt(2 * E / const.m_e)

def _beta(v): return v / const.c

def rho_ch_wrapper(X):
    rho0_ch, a, b = X
    return lambda r: rho0_ch / (1 + np.exp((r-a)/b))

def mott_wrapper(Zp, Zt, alpha, E):
    beta = _beta(_electron_speed(E))
    A = Zp**2 * Zt**2 * alpha**2 * (const.hbar * const.c)**2 / (4 * beta**4 * E)
    def mott_inner(theta):
        sin2_recip = 1/np.sin(theta/2)**2
        return A * (sin2_recip**2 - beta**2 * sin2_recip)
    return mott_inner

def form_factor_wrapper(Z,E,X):
    p = math.sqrt(2 * E * const.m_e)
    rho_ch = rho_ch_wrapper(X)
    A = 4 * const.pi * const.hbar / Z / const.e # TODO: dont forget q!!!!
    def form_factor_inner(theta):
        q = 2 * p * np.sin(theta)
        def f(single_theta):
            return r * rho_ch(r) * np.sin(q * r/ const.hbar)
        int_res = np.array[integrate.quad(f, 0, np.inf)
        return int_res 
    return form_factor_inner

def theo_cross_section_wrapper(Zp, Zt, alpha, E, X):
    mott = mott_wrapper(Zp, Zt, alpha, E) # TODO: optimize by pre-calulating this X-independent factor
    form_factor = form_factor_wrapper(Zt, E, X) # TODO: verify that Zt (target) is correct 
    return lambda theta: mott(theta) * np.abs(form_factor(theta))**2

def Z_t_theoretical(X):
        rho_ch = rho_ch_wrapper(X)
        y, _ = integrate.quad(lambda r: r**2 * rho_ch(r), 0, np.inf)
        return 4*const.pi/const.e * y

def cost_function_wrapper(Zp, Zt, epsilon, alpha, E, data):
    def lsq_function(X):
        Z_t_eval = Z_t_theoretical(X)
        theta, exp_cross, exp_error = data
        Z_term = (Zt - Z_t_eval)**2 / epsilon
        F = theo_cross_section_wrapper(Zp, Zt, alpha, E, X)(theta)
        return np.sum((F - exp_cross)**2 / exp_error**2) + Z_term
    return lsq_function



E = 250e6 * const.e
alpha = const.e
Zp = 1
Zt = 20
eps = 10e-6

# import data
theta, slc, err = np.loadtxt('data.txt', skiprows=2).T
cost_fun = cost_function_wrapper(Zp, Zt, eps, alpha, E, [np.deg2rad(theta), slc, err])

res = optim.least_squares(cost_fun, (1, 0, 1), bounds=([0, 0,0 ], [np.inf, np.inf, np.inf]))
print(res)

fig, ((data_ax,_), _) = plt.subplots(2,2)
data_ax.scatter(theta, slc, 
        label="Experiment (Frosch et al. (1968))")
data_ax.set_yscale('log')
data_ax.set_ylabel(r'$\mathrm{d}\sigma / \mathrm{d}\Omega$ (mb/sr)')
data_ax.set_xlabel(r'$\theta$ (deg)')
data_ax.grid()
data_ax.legend()
data_ax.set_title("Experimental data")
plt.get_current_fig_manager().window.showMaximized()

plt.show()

