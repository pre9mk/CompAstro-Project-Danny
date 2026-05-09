"""
Module that calculates the potentials and accelerations of the following components:

Miyamoto-Nagai Disk: 
def Miyamoto-Nagai(M, a, b)

Hernquist Bulge: 
def hernquist(M, c)

Navarro-Frenk-White Dark Matter Halo:
def nfw_halo(M, r_s)

The component potentials and accelerations are then combined into a total and returned as a Python dictionary.

Tools:
def evaluate(x, y, z) evaluates and returns the potentials as single floats
def acceleration(x, y, z)  evaluates and returns the 3D acceleration vectors as 1D numpy arrays

Author: Danny Western
Institution: University of Virginia, Department of Astronomy
Course: ASTR 5470 - Computational Astrophysics (Dr. Shane Davis)

"""

import numpy as np

G = 4.498e-12

def MiyamotoNagai(M, a, b):
    """
    Creates a Miyamoto-Nagai disk potential

    Parameters:
    Mass of disk: M (solar masses)
    Scale length of disk: a (kpc)
    Scale height of disk: b (kpc)

    Returns (as a Python dictionary):
    Potential
    Acceleration (numpy array)
    
    """

    #Evaluates and returns the potential
    def evaluate(x, y, z):
        R_sq = x**2 + y**2
        z_term = np.sqrt(z**2 + b**2)
        denominator = np.sqrt(R_sq + (a + z_term)**2)

        potential = -(G * M) / denominator

        return potential

    #Evaluates acceleration and returns a numpy array
    def acceleration(x, y, z):
        R_sq = x**2 + y**2
        z_term = np.sqrt(z**2 + b**2)
        a_z_term = a + z_term

        #Calculate base denominator shared by all the components
        base_denominator = (R_sq + a_z_term**2)**1.5

        #x and y accelerations
        ax = -(G * M * x) / base_denominator
        ay = -(G * M * y) / base_denominator

        #z acceleration w/ extra factor
        az_numerator = -(G * M * z * (a + z_term))
        az_denominator = z_term * base_denominator
        az = az_numerator / az_denominator

        return np.array([ax, ay, az])

    return {"evaluate": evaluate, "acceleration": acceleration}

def hernquist(M, c):
    """
    Creates a Hernquist bulge potential

    Parameters:
    Mass of central bulge: M (solar masses)
    Scale radius of the bulge: c (kpc)

    Returns (as a Python dictionary):
    Potential
    Acceleration (numpy array)
    
    """
    
    def evaluate(x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2)

        potential = -(G * M) / (r + c)
        return potential

    def acceleration(x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2)

        if r == 0:
            return np.array([0, 0, 0])

        magnitude = -(G * M) / ((r + c)**2)

        ax = magnitude * (x/r)
        ay = magnitude * (y/r)
        az = magnitude * (z/r)

        return np.array([ax, ay, az])

    return {"evaluate": evaluate, "acceleration": acceleration}

def nfw_halo(M, r_s):
    """
    Creates a Navarro-Frenk-White dark matter halo potential

    Parameters:
    Mass of halo: M (solar masses)
    Scale radius of halo: r_s (kpc)

    Returns (as a Python dictionary):
    Potential
    Acceleration (numpy array)
    
    """

    def evaluate(x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2)

        if r == 0:
            return -(G * M) / r_s

        potential = -(G * M/r) * np.log(1 + (r/r_s))
        return potential

    def acceleration(x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2)
        
        if r == 0:
            return np.array([0.0, 0.0, 0.0])
            
        ratio = r / r_s
        term1 = ratio / (1 + ratio)
        
        
        term2 = np.log(1 + ratio)
        
        magnitude = (G * M / r**2) * (term1 - term2)
        
        ax = magnitude * (x / r)
        ay = magnitude * (y / r)
        az = magnitude * (z / r)
        
        return np.array([ax, ay, az])

    return {"evaluate": evaluate, "acceleration": acceleration}


def total_potential(components):
    """
    Takes the potentials and returns a single unified evaluate and acceleration function

    Parameters:
    Components: All 3 individual components of the galaxy

    Returns (as a Python dictionary):
    Total Potential
    Total Acceleration (numpy array)
    
    """

    def evaluate(x, y, z):
        total_pot = 0

        for comp in components:
            total_pot += comp["evaluate"](x, y, z)

        return total_pot

    def acceleration(x, y, z):
        total_acc = np.array([0.0, 0.0, 0.0])

        for comp in components:
            total_acc += comp["acceleration"](x, y, z)

        return total_acc

    return {"evaluate": evaluate, "acceleration": acceleration}
