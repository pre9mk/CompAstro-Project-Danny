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


#NFW halo test
print("nfw test start")

halo_pot = nfw_halo(M=1e12, r_s=16)
halo_eval = halo_pot["evaluate"]

phi_nfw = halo_eval(16, 0, 0)

print("phi_nfw = ", phi_nfw)


#hernquist test
bulge_pot = hernquist(M=1e10, c=1)
bulge_eval = bulge_pot["evaluate"]
bulge_accel = bulge_pot["acceleration"]

phi_hern = bulge_eval(1, 0, 0)
accel_hern = bulge_accel(1, 0, 0)

print("PE: ", phi_hern)
print("x acceleration: ", accel_hern[0])


"""miyamoto nagai test error log (after fixing fatal errors):
run1 - bad, PE should be negative

Test Start
Potential energy:  0.026945360907484145
acceleration:  -0.003094304799293085

run2 - better, signs are fixed
Test Start
Potential energy:  -0.026945360907484145
acceleration:  -0.003094304799293085

calculating by hand gives: 
distance formula = 8.35 kpc --> G * M = 0.2249 --> phi = -(0.2249/8.3465) = -0.026945 --> ax = (-0.22495 * 8)/581.458 = -0.003094

miyamoto nagai works!

"""



"""nfw test error log
first run fatally errored - did not properly return "evaluate" and "accelerate" (fixed this for hernquist too)

run2 - math error (divided by 0)
RuntimeWarning: divide by zero encountered in log

  potential = -(G * M/r) * np.log(1 - (r/r_s))

phi_nfw =  inf

run3 - code works (and verified by doing the calculation myself)
nfw test start
phi_nfw =  -0.19486100113491459

"""



"""hernquist bulge error log
run1 - actually ran, no fatal errors. math error is clearly present
PE:  0.02249
x acceleration:  -0.02249

run2 - good, matches calculations
PE:  -0.02249
x acceleration:  -0.011245
"""
