import numpy as np
import matplotlib.pyplot as plt
import warnings

def calculate_mass_for_K(x_values, k):
    mass_values = []

    for x in x_values:
        m = [0, 0]; v = [0, 0]; r = [0, 0.001]
        h = 0.001; d = [x, x]
        p = [k * x**2, k * x**2]

        index = 1

        while d[-1] > 0:
            r_current = r[index] + h
            r.append(r_current)

            m_current = (m[index - 1] + 8 * h * np.pi * r[index]**2 * d[index])
            m.append(m_current)

            v_derivative_current = (2 * (m[index] + 4 * np.pi * r[index]**3 * 
                                d[index]) / (r[index] * (r[index] - 2 * m[index])))
            v_current = (v[index - 1] + v_derivative_current * 2 * h)
            v.append(v_current)

            p_current = (p[index - 1] - h * (d[index] + 
                                             p[index]) * v_derivative_current)
            p.append(p_current)

            d_current = (np.sqrt(p_current) / np.sqrt(k))
            d.append(d_current)

            index += 1

        mass_values.append(m[-2])

    return max(mass_values)

def mytests():
    warnings.filterwarnings("ignore")

    print("Tests for Einstein Part D)\n")
    start_value = 100
    end_value = 300
    num_values = 20

    k_values = np.linspace(start_value, end_value, num_values)
    y_values = []
    maximum_passed = False

    for k in k_values:
        start_point = 1.25 * 10**(-3)/k
        end_point = 2.5 / k
        num_points = num_values

        x_values = np.linspace(start_point, end_point, num_points)

        y = calculate_mass_for_K(x_values, k)
        y_values.append(y)
        if y > 2.14 and not maximum_passed:
            print(f"\nMass: {y:.5f} for k={k:.2f} is the first one above maximum\n") 
            maximum_passed = True
            continue
        print(f"Mass: {y:.5f} for k={k:.2f}")
        

    experimental_maximum = 2.14
    print("\nPlotting M_max vs K for Neutron Stars\n")
    plt.figure(figsize=(12, 8))
    plt.plot(k_values, y_values, label='Max Mass for K')
    plt.plot(k_values, experimental_maximum * np.ones(num_values),
             label='M = 2.14', linestyle='dashed')
    plt.xlabel('K')
    plt.ylabel('M_max(K)')
    plt.title("M_max vs K for Neutron Stars")
    plt.legend()
    plt.tight_layout()
    plt.grid()
    plt.show()


def main():
    mytests()


if __name__ == '__main__':
    main()