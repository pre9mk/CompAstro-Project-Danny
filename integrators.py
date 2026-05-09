"""
Module containing the required methods to solve the equations of motion for a star
within a gravitational field.

This implements a 4th-Order Runge-Kutta (RK4) integrator to advance a stellar state vector
[x, y, z, vx, vy, vz] through time.

Author: Danny Western
Institution: University of Virginia, Department of Astronomy
Course: ASTR 5470 - Computational Astrophysics (Dr. Shane Davis)

"""

import numpy as np

def get_derivatives(t, U, acceleration_func):
    """Calculates the derivative of the state vector U"""

    x, y, z, vx, vy, vz = U

    accel = acceleration_func(x, y, z)

    ax, ay, az = accel[0], accel[1], accel[2]

    return np.array([vx, vy, vz, ax, ay, az])

def rk4(t, U, h, acceleration_func):
    """Advances U by timestep h using RK4"""

    k1 = get_derivatives(t, U, acceleration_func)
    k2 = get_derivatives(t + (h/2), U + (h/2) * k1, acceleration_func)
    k3 = get_derivatives(t + (h/2), U + (h/2) * k2, acceleration_func)
    k4 = get_derivatives(t + h, U + h * k3, acceleration_func)

    U_next = U + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

    return U_next
    