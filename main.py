import numpy as np
import matplotlib.pyplot as plt
from potentials import MiyamotoNagai, hernquist, nfw_halo, total_potential
from integrators import rk4
#Constants
#Length: kiloparsecs kpc
#Mass: solar masses M_sun
#Time: megayears Myr

#Gravitational constant in kpc^3 / (M_sun * Myr^2)
G = 4.498e-12

#User input for values of each component
print("Miyamoto-Nagai Disk")
m_disk = float(input("Enter disk mass (M_sun, e.g., 5e10): "))
a_disk = float(input("Enter disk scale length 'a' (kpc): "))
b_disk = float(input("Enter disk scale height 'b' (kpc): "))
disk_pot = MiyamotoNagai(M=m_disk, a=a_disk, b=b_disk)

print("Hernquist Bulge")
m_bulge = float(input("Enter bulge mass (M_sun): "))
c_bulge = float(input("Enter scale radius 'c' (kpc): "))
bulge_pot = hernquist(M=m_bulge, c=c_bulge)

print("NFW Dark Matter Halo")
m_halo = float(input("Enter halo mass (M_sun): "))
rs_halo = float(input("Enter halo scale radius r_s (kpc): "))
halo_pot = nfw_halo(M=m_halo, r_s=rs_halo)

#combine components
galaxy = total_potential([disk_pot, bulge_pot, halo_pot])
accel_tool = galaxy["acceleration"]
eval_tool = galaxy["evaluate"]

#set initial conditions
x0 = float(input("Initial x position (kpc): "))
v0_kms = float(input("Initial tangential velocity (km/s): "))
v0_units = v0_kms * 0.0010227 #converting to kpc/Myr

#initial state vector: [x, y, z, vx, vy, vz]
#z is 0.1 to avoid being perfectly trapped on the x-axis
U_current = np.array([x0, 0, 0.1, 0, v0_units, 0])

#integration setup
num_steps = 50000
h = 1

times = np.zeros(num_steps)
states = np.zeros((num_steps, 6))
energies = np.zeros(num_steps)

#run integration
print("Starting Integration")

for i in range(num_steps):
    times[i] = i * h
    states[i] = U_current

    x, y, z, vx, vy, vz = U_current
    v_squared = vx**2 + vy**2 + vz**2
    phi_total = eval_tool(x, y, z)

    energies[i] = (0.5 * v_squared) + phi_total

    #Use rk4 to move the star forward in time
    U_current = rk4(t=times[i], U=U_current, h=h, acceleration_func=accel_tool)

#make plots
print("Integration complete. Generating plots...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

#plot 1: x vs y
x_coords = states[:, 0]
y_coords = states[:, 1]

ax1.plot(x_coords, y_coords, color="blue", linewidth=0.5)
ax1.set_title("Stellar Orbit (x-y plane)")
ax1.set_xlabel("x Position (kpc)")
ax1.set_ylabel("y Position (kpc)")
ax1.grid(True)
ax1.axis("equal")

#plot 2: energy error
E_initial = energies[0]
energy_error = np.abs((energies - E_initial) / E_initial)

ax2.plot(times, energy_error, color="red", linewidth=1)
ax2.set_title("Energy Conservation Error")
ax2.set_xlabel("Time (Myr)")
ax2.set_ylabel("|dE| / E0")
ax2.set_yscale("log")
ax2.grid(True)

plt.tight_layout()
plt.show()

#sanity check - miyamoto nagai
#print("Test Start")

#value sources: 
#https://science.nasa.gov/missions/hubble/what-does-the-milky-way-weigh-hubble-and-gaia-investigate/ - mass
#https://arxiv.org/abs/astro-ph/9710197
#https://www.mdpi.com/2075-4434/11/3/77
#potential = MiyamotoNagai(M=5e10, a=2.1, b=0.28)

#eval_tool = potential["evaluate"]
#accel_tool = potential["acceleration"]

#test a particle sitting on the x-axis at 8 kpc
#phi = eval_tool(8, 0, 0)
#accel = accel_tool(8, 0, 0)

#print("Potential energy: ", phi)
#print("acceleration: ", accel[0])


#NFW halo test
#print("nfw test start")

#halo_pot = nfw_halo(M=1e12, r_s=16)
#halo_eval = halo_pot["evaluate"]

#phi_nfw = halo_eval(16, 0, 0)

#print("phi_nfw = ", phi_nfw)


#hernquist test
#bulge_pot = hernquist(M=1e10, c=1)
#bulge_eval = bulge_pot["evaluate"]
#bulge_accel = bulge_pot["acceleration"]

#phi_hern = bulge_eval(1, 0, 0)
#accel_hern = bulge_accel(1, 0, 0)

#print("PE: ", phi_hern)
#print("x acceleration: ", accel_hern[0])


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
