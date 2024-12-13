#%%
import astropy.constants as const
import astropy.units as u
import numpy as np

# calculate critical density 
H0 = 70.0 * u.km / u.s / u.Mpc
#convert to s^-1
H0 = H0.to(1/u.s)

rho_crit = 3 * H0**2 / (8 * np.pi *  const.G)
print(rho_crit.to(u.g/u.cm**3))

omega_bar = 0.05
hydrogen_frac = 0.76
n_h = (omega_bar * rho_crit * hydrogen_frac / const.m_p).to(u.cm**-3)

z = 1200 

n_h = n_h * (1 + z)**3

print(n_h)

# photon density
T = 2.725 * u.K
T_z = T * (1 + z)
n_photon = 0.243 * (const.k_B * T_z)**3 / (const.hbar * const.c)**3
n_photon = n_photon.to(u.cm**-3)
print("{:.2e}".format(n_photon))

print("{:.2e}".format(  ((1 + z)**3 * 1.6 * u.m**-3).to(u.cm**-3)  ) ) 

print("baryon to photon ratio: {:.2e}".format(n_h/n_photon))

#%% now do the same for helium
z = 1600
helium_frac = 0.24
m_he = 4 * const.m_p
n_he = (omega_bar * rho_crit * helium_frac / m_he).to(u.cm**-3)
n_he = n_he * (1 + z)**3
print(n_he)

# photon density
T = 2.725 * u.K
T_z = T * (1 + z)
n_photon = 0.243 * (const.k_B * T_z)**3 / (const.hbar * const.c)**3
n_photon = n_photon.to(u.cm**-3)
print("{:.2e}".format(n_photon))


print("{:.2e}".format(  ((1 + z)**3 * 1.6 * u.m**-3).to(u.cm**-3)  ) ) 

print("baryon to photon ratio: {:.2e}".format(n_he/n_photon))

#%%
x = 0.5 
n_p = n_h * x /(1 - x)
exp_term = np.exp((-13.6 * u.eV) / (const.k_B * T_z))
print(n_p * exp_term  * ((const.m_p* const.k_B * T_z)/ (2 * np.pi * const.hbar**2))**(-1.5))