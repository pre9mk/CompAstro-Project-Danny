import numpy as np
from potentials import MiyamotoNagai, hernquist, nfw_halo
from integrators import rk4
#Constants
#Length: kiloparsecs kpc
#Mass: solar masses M_sun
#Time: megayears Myr

#Gravitational constant in kpc^3 / (M_sun * Myr^2)
G = 4.498e-12

#sanity check - miyamoto nagai
print("Test Start")

#value sources: 
#https://science.nasa.gov/missions/hubble/what-does-the-milky-way-weigh-hubble-and-gaia-investigate/ - mass
#https://arxiv.org/abs/astro-ph/9710197
#https://www.mdpi.com/2075-4434/11/3/77
potential = MiyamotoNagai(M=5e10, a=2.1, b=0.28)

eval_tool = potential["evaluate"]
accel_tool = potential["acceleration"]

#test a particle sitting on the x-axis at 8 kpc
phi = eval_tool(8, 0, 0)
accel = accel_tool(8, 0, 0)

print("Potential energy: ", phi)
print("acceleration: ", accel[0])

