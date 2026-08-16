# Fundamental Constants in IS
const ELECTRON_ENERGY_AT_REST_eV = 0.511e6 # eV electron energy at rest
const EMIN = 1e6 # minimum eV energy on the energy spectrum Padovani2021
const EMAX = 1e11 # maximum eV energy on the energy spectrum Padovani2021
const E_0 = 710e6 # eV Padovani 2021 SKA cf eq6
const C_m = 2.99792458e8 # speed of light in m.s^-1

# in CGS
const E_CHARGE = 4.8032e-10 # electron charge in cm^3/2 g^1/2 s^-1 (esu-CGS)
const M_e = 9.109e-28 # electron mass in g
const M_p = 1.6726231e-24 # proton mass in g
const C = 2.99792458e10 # speed of light in cm.s^-1
const K_B = 1.380649e-16 # Boltzmann constant in cm^2 g s^-2 K^-1
const J_0 = 2.1e18 # prefactor of j_e in e^-1 s^-1 cm^-2 sr^-1, Padovani,Galli 2018 SKA
const RM_PREFACTOR = 0.81 #microG^-1 pc^-1 cm^-3

# Unit conversion
const EV_TO_ERG = 1.60218e-12 # eV->erg conversion factor
const JY_TO_MJY = 1e-6 # Jansky to MegaJansky conversion factor
const Wm2Hz_to_Jy = 1e26 # W.m^2.Hz to Jansky conversion factor
const MILI_TO_MICRO = 1e3 # mili -> mu conversion factor
const PARSEC_TO_CM = 3.086e18 # pc -> cm conversion factor

# prefactors
const PRE_NU_C = 3 * E_CHARGE / (4 * pi * M_e * C) # prefactor = 4.19 MHz Padovani,Galli 2018 SKA cf eq 3
const PRE_P = (sqrt(3) * E_CHARGE^3) / (2 * M_e * C^2) # prefactor eq.2 Padovani+ 2021
