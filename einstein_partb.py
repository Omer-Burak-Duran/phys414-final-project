import numpy as np
import matplotlib.pyplot as plt
import warnings
from einstein_parta import calculate_mass_radius


def mytests():
    warnings.filterwarnings("ignore")

    print("Tests for Einstein Part B)\n")

    start_point = 10**(-5)
    end_point = 5 * 10**(-2)
    num_values = 100

    x_values = np.linspace(start_point, end_point, num_values)

    _, radius_values, fractional_binding_values = calculate_mass_radius(x_values)

    print("Plotting Fractional Binding Energy Values for Neutron Stars\n")
    plt.figure(figsize=(12, 8))
    plt.plot(radius_values, fractional_binding_values)
    plt.xlabel('R')
    plt.ylabel('Fractional Binding Energy')
    plt.title("Fractional Binding Energy Values for Neutron Stars")
    plt.grid()
    plt.tight_layout()
    plt.show()


def main():
    mytests()


if __name__ == '__main__':
    main()