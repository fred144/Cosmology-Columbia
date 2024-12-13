#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


H0 = 70  # Hubble constant in km/s/Mpc
c = 299792.458  # Speed of light in km/s
Omega_m = 0.3  # Matter density parameter
Omega_Lambda = 0.7  # Dark energy density parameter

# Hubble parameter as a function of redshift z
def H(z):
    return H0 * np.sqrt(Omega_m * (1 + z)**3 + Omega_Lambda)

# Comoving distance integral
def comoving_distance_integrand(z):
    return c / H(z)

# Compute comoving distance r(z)
def comoving_distance(z):
    r, _ = quad(comoving_distance_integrand, 0, z)
    return r  # in Mpc

# Generate data for z and v_rec
z_values = np.linspace(0, 6, 500)  # Redshift range
r_values = np.array([comoving_distance(z) for z in z_values])
v_rec_values = H0 * r_values  # Recessional velocity in km/s

# Plotting

fig,ax = plt.subplots(1,2,figsize=(12, 5))

ax[1].plot(z_values, r_values, label=r"$r(z)$")
ax[1].set(xlabel="Redshift $z$", ylabel="Comoving Distance $r(z)$ (Mpc)")

exceeds_c = np.where(v_rec_values > c)[0][0]

ax[0].axhline(y=c, color='r', linestyle='--')
ax[0].axvline(x=z_values[exceeds_c], color='r', linestyle='--', label= "$z = {:.2f}$".format(z_values[exceeds_c]) )
ax[0].plot(z_values, v_rec_values, label=r'$\Omega_m = 0.3, \Omega_\Lambda = 0.7$')
ax[0].set(xlabel="Redshift $z$", ylabel="Recessional Velocity $v_{\mathrm{rec}}$ (km/s)", ylim=(0, 8e5))

# plot local relation for recession velocity
ax[0].plot(z_values, c*z_values, label=r"Hubble's law $v_{\mathrm{rec}} = c z$", color='k', linestyle=":")

ax[0].legend()
plt.show()
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


H0 = 70  # Hubble constant in km/s/Mpc
H0_s = (H0 * 1000) / 3.086e22  # H0 in s^-1

# Function to compute age of universe for a given omega_lambda
def correct_age_of_universe(omega_lambda):
    omega_m = 1 - omega_lambda
    # Define the integrand for t_0 as a function of omega_lambda
    def integrand(a, omega_m, omega_lambda):
        return 1 / (a**(3/2) * np.sqrt(omega_m * a**(-3) + omega_lambda))

    # Numerically integrate to find t_0
    result, _ = quad(integrand, 0, 1, args=(omega_m, omega_lambda))
    
    # t_0 is the integral divided by H0 (converted to Gyr)
    return result   # Convert from seconds to Gyr

# Create an array of omega_lambda values from 0 to 1
omega_lambda_vals = np.linspace(1e-3, 1 - 1e-3, 100)

# Compute the age of the universe for each value of omega_lambda
t_0_vals = [correct_age_of_universe(omega_lambda) for omega_lambda in omega_lambda_vals]

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

ax.plot(omega_lambda_vals, np.array(t_0_vals) / ( H0_s * 1e9 * 3.154*1e7)  , label=r"$t_0$ vs. $\Omega_\Lambda$", color='blue')
ax.axhspan(0,15, alpha=0.1, color='red')
ax.set_xlabel(r'$\Omega_\Lambda$', fontsize=14)
ax.set_ylabel(r'$t_0$ (Gyr)', fontsize=14)
ax.axhline(y=13.8, color='r', linestyle='--', label=r'$t_0 = 13.8$ Gyr')
ax.set_title(r'Age of the Universe $t_0$ vs. $\Omega_\Lambda$', fontsize=16)
ax.axvline(x=0.7, color='r', linestyle='--', label=r'$\Omega_\Lambda = 0.7$')
ax.legend()
plt.show()

#%%
def time (omega_m, omega_lambda, h0=70):
    h0 = h0 * (1 / 3.08e19)
    t = 2 * np.arcsinh(np.sqrt(omega_lambda/omega_m)) / (3 * h0 * np.sqrt(omega_lambda)) 
    return t / (3.154e+7 * 1e9)

fig,ax = plt.subplots(1,1,figsize=(8, 6))

omega_lambda_vals = np.linspace(1e-3, 1 - 1e-3, 100)
omega_m_vals = 1 - omega_lambda_vals

t_0 = time(omega_m_vals, omega_lambda_vals)
ax.plot(omega_lambda_vals, t_0, label=r"$t_0$ vs. $\Omega_\Lambda$")
# ax.plot(omega_lambda_vals, np.array(t_0_vals) / ( H0_s * 1e9 * 3.154*1e7)  , label=r"$t_0$ vs. $\Omega_\Lambda$")
ax.axhspan(0,15, alpha=0.1, color='red')
ax.set(xlabel=r"$\Omega_\Lambda$", ylabel=r"$t_0$ (Gyr)")
# ax.plot(omega_lambda_vals, 1 / (( H0_s * 1e9 * 3.154*1e7) * omega_lambda_vals**(1/2) )  , label=r"$t_0$ vs. $\Omega_\Lambda$")
ax.axhline(y=13.8, color='r', linestyle='--', label=r'$t_0 = 13.8$ Gyr')
ax.set_title(r'Age of the Universe $t_0$ vs. $\Omega_\Lambda$', fontsize=16)
ax.axvline(x=0.7, color='r', linestyle='--', label=r'$\Omega_\Lambda = 0.7$')