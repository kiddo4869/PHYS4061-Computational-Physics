import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd

def plot(filename):
    data = np.loadtxt(f"{filename}.txt", skiprows=1)

    time = data[:,1]
    total_energy = data[:,2]
    potential_energy = data[:,3]
    kinetic_energy = data[:,4]
    k1 = data[:,5]
    k2 = data[:,6]
    k3 = data[:,7]

    # Plot the potential energy
    plt.clf()
    plt.plot(time, potential_energy, label="Potential Energy")
    plt.legend()

    plt.xlabel("Time (year)")
    plt.ylabel("Energy (M_sun AU^2 / year^2)")
    plt.title("Potential Energy vs Time")
    plt.savefig(f"{filename}_PE.png")

    # Plot the kinetic energy
    plt.clf()
    plt.plot(time, kinetic_energy, label="Total Kinetic Energy")
    plt.plot(time, k1, label="Kinetic Energy 1")
    plt.plot(time, k2, label="Kinetic Energy 2")
    plt.plot(time, k3, label="Kinetic Energy 3")
    plt.legend()

    plt.xlabel("Time (year)")
    plt.ylabel("Energy (M_sun AU^2 / year^2)")
    plt.title("Kinetic Energy vs Time")
    plt.savefig(f"{filename}_KE.png")

    # Plot the total energy
    plt.clf()
    plt.plot(time, total_energy, label="Total Energy")
    plt.legend()

    plt.xlabel("Time (year)")
    plt.ylabel("Energy (M_sun AU^2 / year^2)")
    plt.title("Total Energy vs Time")
    plt.savefig(f"{filename}_E.png")


def main():
    plot("euler_energies")
    plot("leapfrog_energies")

if __name__ == '__main__':
    main()