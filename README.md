The main goal of this project is to provide easy access to semiconductor
band parameters for calculations and simulations. Basic functionality
requires only the standard python distribution.

Example scripts are provided for basic usage and for generating common
plots such as bandgap vs. lattice constant and bandgap vs. alloy
composition.

Materials included in this version:
- III-V Zinc Blendes
    - Binaries
        - AlN, GaN, InN,
          AlP, GaP, InP,
          AlAs, GaAs, InAs,
          AlSb, GaSb, InSb
    - Ternaries
        - AlGaN, AlInN, GaInN,
          AlGaP, AlInP, GaInP,
          AlGaAs, AlInAs, GaInAs,
          AlGaSb, AlInSb, GaInSb,
          AlNP, GaNP, InNP,
          AlNAs, GaNAs, InNAs,
          AlPAs, GaPAs, InPAs,
          AlPSb, GaPSb, InPSb,
          AlAsSb, GaAsSb, InAsSb
    - Quaternaries
        - AlNPAs, AlPAsSb,
          GaNPAs, GaPAsSb,
          InNPAs, InPAsSb,
          AlGaInN, AlGaInP, AlGaInAs, AlGaInSb,
          AlGaNP, AlInNP, GaInNP,
          AlGaNAs, AlInNAs, GaInNAs,
          AlGaPAs, AlInPAs, GaInPAs,
          AlGaPSb, AlInPSb, GaInPSb,
          AlGaAsSb, AlInAsSb, GaInAsSb

Parameters included in this version:
- lattice constants
- thermal expansion coefficients
- bandgap energies (direct and indirect)
- Varshni parameters
- split-off energies
- effective masses
- Luttinger parameters
- Kane parameters (Ep and F)
- Valance band offsets
- band deformation potentials
- elastic constants
- alloy bowing parameters
- effects of biaxial strain
- optical refractive index based on doi: 10.1088/0965-0393/4/4/002

The [source code](http://github.com/scott-maddox/openbandparams) and [documentation](http://scott-maddox.github.io/openbandparams) are graciously hosted by GitHub.

## Known issues

1. The quaternary alloys interpolation scheme must be improved. If we use as a benchmark the GaInPAs quaternary, the current estimate for the bandgap energy for an InP lattice matched alloy with an As concentration of 0.53 is 1.06eV, when it is well known to be 1.00eV. The current implementation follows a weighted sum of terms, rather than a polynomial expansion. A good place to start is this: https://pubs.aip.org/aip/jap/article/136/21/215105/3324155/Interpolation-of-compound-semiconductor-alloy

## Points to improve

### Include III-V mobility values

The package is built in such a way that the expansion to other material properties can be made seamlessly. Therefore, it would be quite advantageous to extend this package to include mobility values from https://pubs.aip.org/aip/jap/article/87/6/2890/489121/Empirical-low-field-mobility-model-for-III-V.