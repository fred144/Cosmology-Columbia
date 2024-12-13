#%%
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import WMAP7
from astropy import units as u
import astropy.constants as c

theta = np.linspace(0, 2*np.pi, 100)
k=1
G = 1

rho_0 = WMAP7.critical_density(0).to(u.kg/u.m**3)
a = (4*np.pi*G*rho_0)/(3*k) * (1-np.cos(theta))
t = (4*np.pi*G*rho_0)/(3*k*1.5)*(theta-np.sin(theta))

fig, ax = plt.subplots(1,2, figsize=(10,5), dpi=300)
ax[0].plot(theta/np.pi,a, label=r'$a(\theta)$')
ax[0].plot(theta/np.pi,t, label=r'$t(\theta)$', )
ax[0].set(xlabel=r'$\theta [\pi]$', ylabel='values')
ax[1].plot(t,a,label="big crucnh :o")
ax[1].set(xlabel=r'$t$ [arbitrary units]', ylabel='a(t)')
ax[0].legend(title=r"k = 1, G=1, c=1, $\rho_0$ = {:.2}".format(rho_0), frameon=False)
ax[1].legend()
plt.show()