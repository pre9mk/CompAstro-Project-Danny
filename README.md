# Orbit Integrator
A repository containing the code for my ASTR 5470 Computational Astrophysics final project

## Overview
This code is a multi-component numerical solver that simulates stellar trajectories within a 3D galactic potential.
We use a 4th-Order Runge-Kutta (RK4) integrator to calculate the orbital path of a star influenced by three components:

1. A Miyamoto-Nagai Disk
2. A Hernquist Bulge
3. A Navarro-Frenk-White (NFW) Dark Matter Halo

## Dependencies
This code utilizes Python 3 and the following libraries:

1. NumPy - https://numpy.org/doc/
2. Matplotlib - https://matplotlib.org/stable/index.html

## Files
We include 4 primary Python (.py) files:

1. main.py: Main execution file that takes user input for variables, then calls functions from integrators.py and potentials.py to calculate a simulated trajectory 
2. test.py: Test file with three tests - one for each component. Defaulted, simulated variables are used in the calculation and are compared with already calculated expected values
3. integrators.py: File containing the RK4 integrator and a function to take derivatives of input variables
4. potentials.py: File containing the equations for each component. The potentials and accelerations are calculated for each component and returned as numpy arrays.

The two executable files can be run by entering the following in the command line:
python main.py

OR

python test.py



