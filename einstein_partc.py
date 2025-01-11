import numpy as np
import matplotlib.pyplot as plt
import warnings
from einstein_parta import calculate_mass_radius


def mytests():
    warnings.filterwarnings("ignore")

    print("Tests for Einstein Part C)\n")

    start_point = 10**(-5)
    end_point = 5 * 10**(-2)
    num_values = 100

    x_values = np.linspace(start_point, end_point, num_values)

    cutoff_point = 27
    cutoff_value = 1.3238

    mass_values, _, _ = calculate_mass_radius(x_values)

    stable_x_values = x_values[:cutoff_point]
    stable_y_values = mass_values[:cutoff_point]

    unstable_x_values = x_values[cutoff_point:]
    unstable_y_values = mass_values[cutoff_point:]

    highest_y_values = cutoff_value * np.ones(num_values) 


    print("Plotting Stable and Unstable points for Neutron Stars\n")
    plt.figure(figsize=(12, 8))
    plt.plot(stable_x_values, stable_y_values, label='Stable')
    plt.plot(unstable_x_values, unstable_y_values, 
             label='Unstable', linestyle='dashed')
    plt.plot(x_values, highest_y_values, 
             label=f'Highest Stable M = {cutoff_value:.2f}')
    plt.xlabel('rho_c')
    plt.ylabel('Mass')
    plt.title("Stable and Unstable points for Neutron Stars")
    plt.legend()
    plt.grid()
    plt.tight_layout
    plt.show()


def main():
    mytests()


if __name__ == '__main__':
    main()