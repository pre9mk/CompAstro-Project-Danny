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

#Put data in a .txt file and export it
print("Exporting Data...")
filename = "stellar_orbit_data.txt"

with open(filename, "w") as f:
    #Header row
    f.write("Time(Myr)    x(kpc)    y(kpc)    z(kpc)    vx(kpc/Myr)    vy(kpc/Myr)    vz(kpc/Myr)    Energy_Error\n")

    #Loop through and return the first 1000 lines of data
    for i in range(0, num_steps, 50):
        t_val = times[i]
        x_val, y_val, z_val, vx_val, vy_val, vz_val = states[i]
        e_err = energy_error[i]

        f.write(f"{t_val:<12.1f} {x_val:<11.6f} {y_val:<11.6f} {z_val:<11.6f} "
                f"{vx_val:<12.6f} {vy_val:<12.6f} {vz_val:<12.6f} {e_err:<12.6e}\n")

print(f"Data successfully saved to {filename}")