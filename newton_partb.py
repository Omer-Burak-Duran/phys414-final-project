import numpy as np
import matplotlib.pyplot as plt
import csv


# function to read the data
# extracts and returns the values of log_g and mass
# sorted according to mass values
def read_csv_data():
    wd_data = open('white_dwarf_data.csv')
    reader = csv.reader(wd_data, delimiter=',')
    skip_header = True
    log_g_values = []
    mass_values = []
    for row in reader:
        if skip_header:
            skip_header = False
            continue
        log_g_values.append(float(row[1]))
        mass_values.append(float(row[2]))
    zipped_values = list(zip(mass_values, log_g_values))
    sorted_values = list(zip(*sorted(zipped_values, key=lambda x: x[0])))
    return sorted_values[0], sorted_values[1]

# function to calculate and return the scaled radius values
def calculate_radius(mass_values, log_g_values):
    solar_mass = float(1.988 * 10**(30))
    gravitational_constant = 6.67384e-11
    radius_unit_converter = 1/(6.371 * 10**(6))
    radii_average_earth_radii = np.zeros(len(log_g_values))
    for i in range(len(log_g_values)):
        g_SI = 10**(log_g_values[i] - 2)
        radii_average_earth_radii[i] = \
            (np.sqrt(gravitational_constant * solar_mass * 
                     mass_values[i] / g_SI) * radius_unit_converter)
    return radii_average_earth_radii


def mytests():
    print("Tests for Newton Part B)\n")
    # read the csv data
    mass_values, log_g_values = read_csv_data()
    
    # calculate the radius values
    radius_values = calculate_radius(mass_values, log_g_values)

    # plot mass vs radius
    print("Plotting Mass vs Radius for White Dwarf Data")
    plt.figure(figsize=(12, 8))
    plt.plot(radius_values, mass_values)
    plt.xlabel('R')
    plt.ylabel('M')
    plt.title('Mass vs Radius for White Dwarf Data')

    plt.tight_layout()
    plt.grid()
    plt.show()

def main():
    mytests()
    
if __name__ == '__main__':
    main()