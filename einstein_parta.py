import numpy as np
import matplotlib.pyplot as plt
import warnings

# calculate and return mass and radius values
def calculate_mass_radius(x_values):

    mass_values = []
    radius_values = []
    m_P_values = []
    fractional_binding_values = []
    
    for x in x_values:
        m = [0, 0]; m_p = [0, 0]; v = [0, 0]
        r = [0, 1e-3]; h = 1e-3; k = 100
        d = [x, x]; p = [k * x**2, k * x**2]

        index = 1
        while d[-1] > 0:
            r_current = r[index] + h
            r.append(r_current)

            m_current = (m[index-1] + 8 * h * np.pi * r[index]**2 * d[index])
            m.append(m_current)

            v_derivative_current = (2 * (m[index] + 4 * np.pi * r[index]**3 * 
                                    d[index]) / (r[index] * (r[index]-2 * m[index])))
            v_current = (v[index-1] + v_derivative_current * 2 * h)
            v.append(v_current)

            p_current = (p[index-1] - h * (d[index] + 
                                           p[index]) * v_derivative_current)
            p.append(p_current)

            d_current = (np.sqrt(p_current) / np.sqrt(k))
            d.append(d_current)

            m_p_current = (m_p[index-1] + 8 * h * np.pi * 
                           (1 - (2 * m[index] / r[index]))**(-1/2) * 
                           r[index]**2 * d[index])
            m_p.append(m_p_current)

            index += 1

        fractional_binding_values.append((m_p[-2] - m[-2]) / m[-2])
        m_P_values.append(m_p[-2])
        mass_values.append(m[-2])
        radius_values.append(r[-2] / 1.5)

    return mass_values, radius_values, fractional_binding_values


def mytests():
    warnings.filterwarnings("ignore")

    print("Tests for Einstein Part A)\n")

    start_point = 10**(-5)
    end_point = 5 * 10**(-2)
    num_values = 100

    x_values = np.linspace(start_point, end_point, num_values)

    mass_values, radius_values, _ = calculate_mass_radius(x_values)

    print("Plotting Mass vs Radius for Neutron Stars\n")
    plt.figure(figsize=(12, 8))
    plt.plot(radius_values, mass_values)
    plt.xlabel('R')
    plt.ylabel('M')
    plt.title("Mass vs Radius for Neutron Stars")
    plt.grid()
    plt.tight_layout()
    plt.show()


def main():
    mytests()


if __name__ == '__main__':
    main()