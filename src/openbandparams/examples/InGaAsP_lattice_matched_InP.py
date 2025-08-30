import openbandparams as obp
import numpy as np
from matplotlib import pyplot as plt
import sys

x=np.linspace(0,1,100)

Eg = 1.35-0.72*x+0.12*x**2

Eg_obp = []
for xi in x:
    mat = obp.GaInPAs(As = xi, a=obp.InP.a())
    Eg_obp.append(mat.Eg())

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, Eg, label='Analytical')
ax.plot(x, Eg_obp, label='OpenBandParams')
ax.set_xlabel('As Composition')
ax.set_ylabel('Bandgap Energy (eV)')
ax.legend()

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        output_filename = sys.argv[1]
        plt.savefig(output_filename)
    else:
        plt.show()