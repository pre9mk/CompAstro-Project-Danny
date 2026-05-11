"""
Testing module with defaulted variable values.

This script mainly serves as a test for potentials.py.
It verifies the accuracy of the code for our three galaxy components by
comparing simulated values with their expected values for the same default parameters.

Does not return values to the overarching program
Does not save any files, but generates a 1D slice of the potential wells
along the x-axis for visual assessment.

Author: Danny Western
Institution: University of Virginia, Department of Astronomy
Course: ASTR 5470 - Computational Astrophysics (Dr. Shane Davis)

"""

import numpy as np
import matplotlib.pyplot as plt
from potentials import MiyamotoNagai, hernquist, nfw_halo

G = 4.498e-12

#Miyamoto-Nagai test
def test_miyamoto_nagai():
    """
    Test function for the Miyamoto-Nagai Disk

    Defaulted Parameters:
    Mass: M = 5e10 solar masses
    Scale length: a = 2 kpc
    Scale height: b = 0.28 kpc
    """
    
    print("Miyamoto-Nagai Test Start")

    #Milky Way Value Sources: 
    #https://science.nasa.gov/missions/hubble/what-does-the-milky-way-weigh-hubble-and-gaia-investigate/ - mass
    #https://arxiv.org/abs/astro-ph/9710197 - scale length
    #https://www.mdpi.com/2075-4434/11/3/77 - scale height

    #Assigns variables to values from the dictionary produced by potentials.py and the defaulted values
    potential = MiyamotoNagai(M=5e10, a=2, b=0.28)
    phi = potential["evaluate"](8, 0, 0)
    accel = potential["acceleration"](8, 0, 0)


    #Hand-calculated values for given parameters
    expected_phi = -0.026945
    expected_ax = -0.003094

    #https://docs.python.org/3/reference/simple_stmts.html - "assert" documentation
    #https://numpy.org/doc/2.3/reference/generated/numpy.isclose.html - numpy.isclose() documentation

    #Checks if the outputted values matche the predetermined expected values
    print("Expected PE: ", expected_phi)
    print("Expected Acceleration: ", expected_ax)
    
    assert np.isclose(phi, expected_phi, rtol=1e4), f"MN PE failed: {phi}"
    assert np.isclose(accel[0], expected_ax, rtol=1e4), f"MN acceleration failed: {ax}"
    print("Simulated PE: ", phi)
    print("Simulated Acceleration: ", accel)
    print("Passed! Output matches/similar to analytic solution")
    print("")

#NFW halo test
def test_nfw():
    """
    Test function for the Navarro-Frenk-White dark matter halo

    Defaulted Parameters:
    M = 1e12 solar masses
    r_s = 16
    """
    
    print("NFW Test Start")

    #https://ui.adsabs.harvard.edu/abs/2012PASJ...64...75S/abstract - Milky Way values
    potential = nfw_halo(M=1e12, r_s=16)
    phi = potential["evaluate"](16, 0, 0)

    expected_phi = -0.194861
    print("Expected PE: ", expected_phi)
    
    assert np.isclose(phi, expected_phi, rtol=1e4), f"NFW PE failed: {phi}"
    
    print("Simulated PE: ", phi)
    print("Passed! Output matches/similar to analytic solution")
    print("")


#Hernquist test
def test_hernquist():
    """
    Test function for the Navarro-Frenk-White dark matter halo

    Defaulted Parameters:
    M = 1e10 solar masses
    c = 1
    """
    
    print("Start Hernquist Bulge Test")

    #https://www.aanda.org/articles/aa/abs/2016/03/aa27500-15/aa27500-15.html - Milky Way values
    potential = hernquist(M=1e10, c=1)
    phi = potential["evaluate"](1, 0, 0)
    accel = potential["acceleration"](1, 0, 0)

    expected_phi = -0.02249
    expected_ax = -0.011245
    print("Expected PE: ", expected_phi)
    print("Expected Acceleration: ", expected_ax)
    
    assert np.isclose(phi, expected_phi, rtol=1e4), f"Hernquist PE failed: {phi}"
    assert np.isclose(accel[0], expected_ax, rtol=1e4), f"Hernquist acceleration failed: {ax}"

    print("Simulated PE: ", phi)
    print("Simulated Acceleration: ", accel)
    print("Passed! Output matches/similar to analytic solution")
    print("")


def plot_potentials():
    """Generates a 1D slice of the potential wells along the x-axis."""
    print("Generating plots")

    #Create an array of x-coordinates from 0.1 kpc to 50 kpc
    x_vals = np.linspace(0.1, 50, 500)

    #Initialize components
    disk_pot = MiyamotoNagai(M=6.82e10, a=3, b=0.28)["evaluate"]
    bulge_pot = hernquist(M=1e10, c=1)["evaluate"]
    halo_pot = nfw_halo(M=1e12, r_s=16)["evaluate"]

    # Calculate the potential energy arrays
    disk_phi = [disk_pot(x, 0, 0) for x in x_vals]
    bulge_phi = [bulge_pot(x, 0, 0) for x in x_vals]
    halo_phi = [halo_pot(x, 0, 0) for x in x_vals]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, disk_phi, label="Miyamoto-Nagai Disk", color="blue")
    plt.plot(x_vals, bulge_phi, label="Hernquist Bulge", color="orange")
    plt.plot(x_vals, halo_phi, label="NFW Halo", color="green")
    
    plt.title("1D Galactic Potential Wells (X-Axis Slice)")
    plt.xlabel("Distance from Galactic Center (kpc)")
    plt.ylabel("Potential Energy (Phi)")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.show()
    

#Run tests
if __name__ == "__main__":
    print("Running Validation Tests")
    test_miyamoto_nagai()
    test_hernquist()
    test_nfw()
    print("All tests passed successfully!")

#Generate the plot of potentials
plot_potentials()
