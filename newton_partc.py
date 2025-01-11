import numpy as np
import matplotlib.pyplot as plt
from newton_partb import read_csv_data, calculate_radius
from numpy.polynomial.polynomial import polyfit


# calculate and plot the power-law dependence
def plot_power_law_dependence(mass_values, radius_values):
    num_values = len(mass_values)

    # cutoff values to try
    cutoff_values = [-0.5, -1.0]

    y_values1 = cutoff_values[0] * np.ones(num_values)
    y_values2 = cutoff_values[1] * np.ones(num_values)

    x_values = np.log(radius_values)
    y_values = list(map(lambda x: np.log(x), mass_values))

    print("Plotting log M vs log R to see power-law dependence\n")
    plt.figure(figsize=(12, 6))
    plt.plot(x_values, y_values1, label=f"y = {cutoff_values[0]}")
    plt.plot(x_values, y_values2, label=f"y = {cutoff_values[1]}")
    plt.plot(x_values, y_values, label="log M values")
    plt.xlabel("log R")
    plt.ylabel("log M")
    plt.title("log M vs log R to see power-law dependence")
    plt.legend()
    plt.tight_layout
    plt.grid()
    plt.show()

# return the filtered data to contain only samples with low-mass
def get_filtered_data(mass_values, log_g_values, radius_values, cutoff=-1.1):
    filtered_mass_values = []
    filtered_log_g_values = []
    filtered_radius_values = []

    num_values = len(radius_values)

    for i in range(num_values):
        if np.log(mass_values[i]) < cutoff:
            filtered_mass_values.append(mass_values[i])
            filtered_log_g_values.append(log_g_values[i])
            filtered_radius_values.append(radius_values[i])

    return filtered_mass_values, filtered_log_g_values, filtered_radius_values


# calculate lane-emden with q = 3 and n = 3/2
def calculate_lane_emden(h = 1e-4):
    count = 1
    previous_f = 1
    current_f = 1
    next_f = 1
    result = [1, 1]

    while next_f > 0:
        next_f = ((h * previous_f + 
                   2 * count * h * current_f - 
                   count * h * previous_f - 
                   h**3 * count * current_f**(3 / 2)) 
                   / (h * (count + 1)))
        
        result.append(next_f)
        previous_f = current_f
        current_f = next_f

        count += 1

    ksi = count * h
    return result, ksi


# calculate and return the value of rho_c
def calculate_rho_c(mass_values, radius_values, h = 1e-4):
    result, ksi = calculate_lane_emden(h)

    d_ksi = (result[-1] - result[-3]) / (2 * h)

    rho_c_values = []
    num_values = len(mass_values)

    for i in range(num_values):
        rho_c_values.append(mass_values[i] * ksi / 
                            (4 * np.pi * radius_values[i]**(3) * (-d_ksi)))
        
    return rho_c_values


# calculate the fit to find the value of q and K_star
def calculate_fit(mass_values, radius_values):
    y_values = np.log(mass_values)
    x_values = np.log(radius_values)

    result = polyfit(x_values, y_values, deg=1)

    A = np.exp(result[0]) 
    q = 5 / 2 - 5 / (2 * result[1] - 4)

    return A, q


# calculate the value of proportionality constant using the fit
def calculate_proportionality_constant(A, h=1e-4):
    solar_mass = float(1.988 * 10 ** 30)
    average_radius = 6.371 * 10 ** 6
    gravitational_constant = 6.67384e-11
    modified_g = (gravitational_constant * solar_mass / average_radius**3)
    constant_term = (2.5)**3 * (np.pi * 4)**(-2)

    result, ksi = calculate_lane_emden(h)

    d_ksi = (result[-1] - result[-3]) / (2 * h)
    B = constant_term * d_ksi * (-(ksi)**5)
    unit_scale = (average_radius * solar_mass**(1 / 3))

    K = modified_g * (A / B)**(1 / 3)
    scaled_K = (gravitational_constant * (A / B)**(1 / 3) * unit_scale)

    return K, scaled_K

# calculate the values of C, D, and K analytically
def analytical_C_D_K():
    atomic_mass = float(1.66053907 * 10**(-27))
    electron_mass = float(9.1093847 * 10**(-31))
    speed_of_light = float(2.99792458 * 10**8)
    reduced_planck_constant = float(1.054571817 * 10**(-34))

    C = ((electron_mass**4 * speed_of_light**5) / 
         (24 * np.pi**2 * reduced_planck_constant**3))
    
    D = ((atomic_mass * electron_mass**3 * speed_of_light**3 * 2) / 
          (3 * np.pi**2 * reduced_planck_constant**3))
    
    K = 1.6 * C / (D**(5 / 3))

    return C, D, K


def mytests():

    print("Tests for Newton Part C)\n")

    # read the csv data
    mass_values, log_g_values = read_csv_data()
    
    # calculate the radius values
    radius_values = calculate_radius(mass_values, log_g_values)

    # plot the power-law dependence
    plot_power_law_dependence(mass_values, radius_values)

    # get the filtered values
    filtered_mass_values, filtered_log_g_values, filtered_radius_values = \
        get_filtered_data(mass_values, log_g_values, radius_values)
    
    # calculate and plot the rho_c values
    rho_c_values = calculate_rho_c(filtered_mass_values, filtered_radius_values)

    print("Plotting rho_c vs Mass for filtered White Dwarf Data\n")
    plt.figure(figsize=(12, 8))
    plt.plot(filtered_mass_values, rho_c_values)
    plt.xlabel('Mass')
    plt.ylabel('rho_c')
    plt.title('rho_c vs Mass for filtered White Dwarf Data')

    plt.tight_layout()
    plt.grid()
    plt.show()

    # calculate the fit and fing the value of q
    A, q = calculate_fit(filtered_mass_values, filtered_radius_values)

    print(f"Calculated value for q: {q}\n" + 
          f"Nearest integer value is: {np.rint(q)}\n")

    # find the value of K_star
    K, scaled_K = calculate_proportionality_constant(A)
    print(f"K in scaled units: {K}\nK in SI units: {scaled_K}\n")

    # calculate the values of C, D, and K analytically
    exact_C, exact_D, exact_K = analytical_C_D_K()
    print("Analytical values for C, D, and K are:\n" + 
          f"C: {exact_C}\nD: {exact_D}\nK: {exact_K}\n")

    # compare the analytical and numerical values of K
    K_error = 100 * np.abs(scaled_K - exact_K) / exact_K
    print(f"Error for K is {K_error:.2f}%\n")


def main():
    mytests()


if __name__ == '__main__':
    main()