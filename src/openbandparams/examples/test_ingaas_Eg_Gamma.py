import openbandparams as obp
import numpy as np
from matplotlib import pyplot as plt

### Adjust bowings
# obp.iii_v_zinc_blende_ternaries[5].Eg_Gamma_bowing = 0.65*10 #GaInP
# obp.iii_v_zinc_blende_ternaries[8].Eg_Gamma_bowing = 0.42*10 #GaInAs
# obp.iii_v_zinc_blende_ternaries[9].Eg_Gamma_bowing = 0.19*10 #GaPAs
# obp.iii_v_zinc_blende_ternaries[10].Eg_Gamma_bowing = 0.1*10 #InPAs

#Test for InGaAs
x_ingaas = np.linspace(0, 1, 10)

#https://www.batop.de/information/Eg_InGaAs.html
Eg_ingaas = 1.425-1.501*x_ingaas+0.436*x_ingaas**2 #amount of In

Eg_ingaas_obp = np.zeros_like(x_ingaas)

for i, xi in enumerate(x_ingaas	):
	mat = obp.GaInAs(In = xi)
	Eg_ingaas_obp[i] = mat.Eg()


# ## Test for InGaP
# print('Eg for In0.49Ga0.51P:', obp.GaInP(In = 0.49, Ga = 0.51).Eg())
# print('Eg for In0.49Ga0.51P from literature: 1.9 eV')


# ## Test for InP
# print(f'Bandgap from InP: {obp.InP.Eg(T=300):.3f}')
# print('Bandgap from InP literature: 1.344 eV')

# ## Test for GaP
# print(f'Bandgap from GaP: {obp.GaP.Eg():.3f}')
# print('Bandgap from GaP literature: 2.26 eV')

# ## Test for GaAs
# print(f'Bandgap from GaAs: {obp.GaAs.Eg():.3f}')
# print('Bandgap from GaAs literature: 1.42 eV')

# ## Test for InAs
# print(f'Bandgap from InAs: {obp.InAs.Eg():.3f}')
# print('Bandgap from InAs literature: 0.354 eV')


## Test for InGaP
x_ingap = np.linspace(0, 1, 10) #Amount of In

Eg_ingap = 2.77*(1-x_ingap) + 1.344*x_ingap - 0.65*x_ingap*(1-x_ingap)
Eg_ingap_obp = np.zeros_like(x_ingap)

for i, xi in enumerate(x_ingap):
	mat = obp.GaInP(In = xi)
	Eg_ingap_obp[i] = mat.Eg_Gamma()

## Test for GaPAs
x_gapas = np.linspace(0, 1, 10) #Amount of P (GaAs_1-xP_x)

Eg_gapas = 1.42*(1-x_gapas) + 2.77*x_gapas - 0.21*x_gapas*(1-x_gapas)
Eg_gapas_obp = np.zeros_like(x_gapas)

for i, xi in enumerate(x_gapas):
	mat = obp.GaPAs(P = xi)
	Eg_gapas_obp[i] = mat.Eg_Gamma()

## Test for InPAs


### Test for InGaAsP
y_ingaasp = np.linspace(0, 1, 10)

#Fiedler, F., Schlachetzki, A., 1987. Optical parameters of InP-based waveguides. Solid-State Electronics 30, 73–83. https://doi.org/10.1016/0038-1101(87)90032-3

Eg_ingaasp = 1.35-0.72*y_ingaasp+0.12*y_ingaasp**2
Eg_ingaasp_vurg = 1.353*(1-y_ingaasp) + 0.737*y_ingaasp + 0.13*y_ingaasp**2

Eg_ingaasp_obp = np.zeros_like(y_ingaasp)

for i, yi in enumerate(y_ingaasp):
	mat = obp.GaInPAs(As = yi, a = obp.InP.a())
	Eg_ingaasp_obp[i] = mat.Eg()

mat = obp.GaInPAs(As = 0.53, a = obp.InP.a())
print(mat._type)
print(mat.binaries)
print(mat.ternaries)
print(mat._xyz)
print(mat.Eg())
y1 = 0.53
print(1.35-0.72*y1+0.12*y1**2)


# plt.xlabel('x (Ga composition)')
# plt.ylabel('Eg (eV)')
# plt.title('Bandgap of GaInAs vs Ga composition')


plt.plot(y_ingaasp, Eg_ingaasp, label='Eg InGaAsP')
plt.plot(y_ingaasp, Eg_ingaasp_obp, label='Eg_obp InGaAsP')
plt.plot(y_ingaasp, Eg_ingaasp_vurg, label='Eg_vurg InGaAsP')

# plt.plot(x_ingaas, Eg_ingaas, label='Eg InGaAs')
# plt.plot(x_ingaas, Eg_ingaas_obp, label='Eg_obp InGaAs')

# plt.plot(x_ingap, Eg_ingap, label='Eg InGaP')
# plt.plot(x_ingap, Eg_ingap_obp, label='Eg_obp InGaP')

# plt.plot(x_gapas, Eg_gapas, label='Eg GaPAs')
# plt.plot(x_gapas, Eg_gapas_obp, label='Eg_obp GaPAs')

plt.legend()
plt.grid()
plt.show()