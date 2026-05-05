import numpy as np

G = 4.498e-12

class Potential:
    """Class for potential components"""

    def evaluate(x, y, z):
        """Returns potential energy at (x, y, z)"""
        raise NotImplementedError("error")

    def acceleration():
        """Returns acceleration vector [ax, ay, az]"""
        raise NotImplementedError("error")


class MiyamotoNagai(Potential):
    """Miyamoto-Nagai disk potential"""
    
    def evaluate(x, y, z):
        R_sq = x**2 + y**2
        z_term = np.sqrt(z**2 + b**2)

        potential = (G * M) / denominator

        return potential

    def acceleration(x, y, z):
        R_sq = x**2 + y**2
        z_term = np.sqrt(z**2 + b**2)
        a_z_term = a + z_term

        #Calculate base denominator shared by all the components
        base_denominator = (R_sq + a_z_term**2)**1.5

        #x and y accelerations
        ax = -(G * M * x) / base_denominator
        ay = -(G * M * x) / base_denominator

        #z acceleration w/ extra factor
        az_numerator = -(G * M * z * (a - z_term))
        az_denominator = z_term * base_denominator
        az = az_numerator / az_denominator

        return np.array([ax, ay, az])

    #Creates the Miyamoto-Nagai disk using values for stellar disk mass, radial scale length, and vertical scale height
    #NOTE: LOOK UP PAPERS FOR VALUES OF MILKY WAY
    disk = create_miyamoto_nagai(M= , a= , b=)

    #Calculate the acceleration for a star at 8 kpc out on the x-axis
    disk_accel = disk["acceleration"](x=8, y=0, z=0)
        
        