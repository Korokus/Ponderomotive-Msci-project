import math
import numpy as np
import cmath
import matplotlib.pyplot as plt
epsilon_inf = 3.7
#plasma frequency in radians per second
plasma_freq = 1.37e16
gamma = 8e13
#refractive index
radius = 20e-9
prefactor = 4 * math.pi * radius**3
sigma_geom = math.pi * radius**2
c = 3e8
def main():
    refractive_indices = np.linspace(1,1.5,6)
    wavelengths = np.linspace(300e-9, 450e-9, 1000)
    values1 = extinction_calc(refractive_indices[0],wavelengths)
    values2 = extinction_calc(refractive_indices[1],wavelengths)
    values3 = extinction_calc(refractive_indices[2],wavelengths)
    values4 = extinction_calc(refractive_indices[3],wavelengths)
    values5 = extinction_calc(refractive_indices[4],wavelengths)
    values6 = extinction_calc(refractive_indices[5],wavelengths)
    graph(values1, values2, values3, values4, values5, values6, wavelengths)
    
def extinction_calc(refractive_index,wavelengths):
    extinction = []
    e1 = e1_calc(refractive_index)
    for i in range(1000):
        wavelength = wavelengths[i]
        w = 2*math.pi * c / wavelength
        e2 = e2_calc(w)
        alpha = prefactor * (e2 -e1)/(e2  + 2*e1)
        alpha_mod = complex_modulo(alpha)
        k_0 = (2*math.pi)/wavelength
        scattering_prefactor = (k_0**4)/(6*math.pi)
        scattering = scattering_prefactor * alpha_mod**2
        absorption = (k_0 * alpha.imag)
        extinction.append((scattering + absorption)/sigma_geom)
    return extinction

def e2_calc(w):
    z = epsilon_inf - (plasma_freq**2)/(w**2 + complex(0,gamma*w))
    return z

def e1_calc(refractive_index):
    e1 = refractive_index**2
    return e1

def complex_modulo(z):
    return math.sqrt((z * np.conjugate(z)))

def graph(values1, values2, values3, values4, values5, values6, wavelengths):
    plt.figure(0)
    plt.plot(wavelengths, values1, label = "1.0")
    plt.plot(wavelengths, values2, label = "1.1")
    plt.plot(wavelengths, values3, label = "1.2")
    plt.plot(wavelengths, values4, label = "1.3")
    plt.plot(wavelengths, values5, label = "1.4")
    plt.plot(wavelengths, values6, label = "1.5")
    plt.title("Extinction of silver nanoparticle")
    plt.ylabel("sig_ext/sig_geom")
    plt.xlabel("wavelength/nm")
    plt.legend(loc="upper right")
    plt.savefig("Extinction of silver nanoparticle", dpi = 300, bbox_inches = 'tight')


if __name__ == '__main__':
    main()