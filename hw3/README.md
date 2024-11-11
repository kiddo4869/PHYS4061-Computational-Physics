# Three Body Simulator

This program simulates the gravitational interactions of a three-body system

## Features

- Simulates using both the Euler method and the leapfrog method
- Outputs trajectory over time
- Outputs the total potential energy, total kinetic energy, and individual kinetic energy over time 

## Requirements

- C++ compiler
- Python with the following libraries:
  - "numpy"
  - "matplotlib"

## Usage

1. Compile the program:
   g++ ThreeBodySimulation.cpp -o ThreeBodySimulation

2. Run the program:
   ./ThreeBodySimulation

3. Input numbers of the row and column of the matrix, and the corresponding elements:
   e.g.
   Please input the time period (in year): 50
   Please input the time step (in year): 0.001

4. Plot the KE vs time and PE vs time in two different integration methods
   python plot.py

## Outputs
   - Trajectory files (in XYZ format):
      - euler.xyz
      - leapfrog.xyz

   - Energy data over time:
      - euler_enegies.txt
      - leapfrog_enegies.txt

   - Energy plots:
      - euler_enegies_E.png (Total Energy)
      - euler_enegies_KE.png (Kinetic Energy)
      - euler_enegies_PE.png (Potential Energy)
      - leapfrog_enegies_E.png (Total Energy)
      - leapfrog_enegies_KE.png (Kinetic Energy)
      - leapfrog_enegies_PE.png (Potential Energy)