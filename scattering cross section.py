import math
import numpy as np
import cmath

epsilon_inf = 3.7
#plasma frequency in radians per second
plasma_freq = 1.37e16
gamma = 8e13
#refractive index
refractive_index = 1
radius = 20e-9
prefactor = 4 * math.pi * radius**4

def main():
    scattering = []
    absorption = []
    extinction = []

    wavelengths = np.linspace(300, 450, 100)
    e1 = e1_calc(refractive_index)
    for i in wavelengths:
        wavelength = wavelengths[i]
        e2 = e2_calc(wavelength)
        alpha = prefactor * (e2 -e1)/(e2  + 2*e1)
        alpha_mod = complex_modulo(alpha)
        k_0 = wavelength/(2*math.pi)
        scattering_prefactor = (k_0**4)/(6*math.pi)
        scattering.append(scattering_prefactor * alpha_mod)
        absorption.append(k_0 * alpha.imag)
        extinction.append(scattering_prefactor * alpha_mod + k_0 * alpha.imag)

    print(scattering)



def e2_calc(frequency):
    z = complex(epsilon_inf - (((plasma_freq**2)*frequency)/(((frequency**4)) + ((gamma**2)* (frequency**2)))),(gamma*frequency*(plasma_freq**2))/((plasma_freq**4) + (gamma**2)* (frequency**2)))
    return z

def e1_calc(refractive_index):
    e1 = complex(refractive_index**2,0)
    return e1

def complex_modulo(z):
    z_real = z.real
    z_imag = z.imag
    return math.sqrt(z_real**2+z_imag**2)



if __name__ == '__main__':
    main()