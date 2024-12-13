# %%
import numpy as np
from astropy import units as u
from astropy import constants as c
import matplotlib.pyplot as plt
from scipy.special import kn
from scipy.special import i1, i0

x = 0.07
v = 220
m = 1
dop = 50
timescale = (
    24.7
    * (m / 0.5) ** 0.5
    * (dop / 4) ** 0.5
    * ((x * (1 - x)) / 0.25) ** 0.5
    * (v / 220) ** -1
)
au = 2.85 * u.au * (m / 0.5) ** 0.5 * (dop / 4) ** 0.5 * ((x * (1 - x)) / 0.25) ** 0.5
au = au.to(u.km).value
timescale_lense = 1e9 / 220  # seconds
timescale_lense = timescale_lense / (60 * 60 * 24)
print("{:2e}".format(au))
print("{:2e}".format(timescale_lense))

# %%
R = np.geomspace(0.01, 1e5, 100) * u.pc
Rd = 8.5e3 * u.pc

surface_density = 1e12 * u.Msun / (2 * np.pi * 63 * u.kpc**2)
surface_density = surface_density.to(u.Msun / u.pc**2)

G = c.G.to(u.pc**3 / (u.Msun * u.s**2))
m_sphere_enc = (
    2 * np.pi * surface_density * (Rd**2 - (R * Rd + Rd**2) * np.exp(-R / Rd))
)
v_sphere = np.sqrt((G * m_sphere_enc) / R)


y = R / (2 * Rd)
v_disk = np.sqrt(
    4 * np.pi * G * surface_density * Rd * y**2 * (i0(y) * kn(1, y) - i1(y) * kn(0, y))
)

fractional_difference = (v_disk - v_sphere) / v_disk

fig, ax = plt.subplots(2, 1, dpi=200, figsize=(6, 4), sharex=True)
ax[0].plot(R /1e3, v_disk.to(u.km / u.s),label='Disk')
ax[0].plot(R /1e3, v_sphere.to(u.km / u.s), label = 'Sphere')
ax[1].plot(R /1e3, fractional_difference, label = 'Sphere')
ax[0].set(xlim=(0, 100), ylim=(0, 1000), ylabel="v (km/s)")
ax[1].set(ylabel="(v_disk - v_sphere) / v_disk", xlabel="r (kpc)")
ax[0].legend()
