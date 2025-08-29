import openbandparams as obp
import numpy as np
from matplotlib import pyplot as plt

x=np.linspace(0,1,100)

Eg = 1.35-0.72*x+0.12*x**2

Eg_obp = []
for xi in x:
    mat = obp.GaInPAs(As = xi, a=obp.InP.a())
    Eg_obp.append(mat.Eg())

plt.plot(x, Eg, label='Analytical')
plt.plot(x, Eg_obp, label='OpenBandParams')
plt.xlabel('As Composition')
plt.ylabel('Bandgap Energy (eV)')
plt.legend()
plt.show()