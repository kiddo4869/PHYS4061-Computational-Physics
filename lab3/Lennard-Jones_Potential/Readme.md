# Lennard-Jones Potential Calculator

This program calculates the Lennard-Jones potential based on the neighbor list of atoms for different crystal structures within a specified cutoff distance.

## Features

- Supports Simple Cubic (SC), Body-Centered Cubic (BCC), Face-Centered Cubic (FCC), and Diamond structures
- Computes and outputs reciprocal lattice vectors for each structure
- Outputs a list of nearest neighbors for each structure to a text file (.txt)
- Calculates the total energy (cohesive energy)

## Requirements

- C++ compiler
- Standard C++ libraries

## Usage

1. Compile the program:
   g++ Lennard-Jones_Potential.cpp -o Lennard-Jones_Potential

2. Run the program:
   ./Lennard-Jones_Potential

3. Input the lattice constant and the number of periods in the x, y, and z directions when prompted:
   e.g.
   Please input the number of period in x y z direction: 4 4 4

4. The program will generate the following files:
   - Neighbor_List_Diamond_OverCounted.txt

## Outputs
The output files contain a table with the following format:
   e.g.
	   Label_1,Label_2,Distance
      0,1,3.75474
      0,2,3.75474
      0,3,3.75474

Each column separated by a comma (,) in the text file
Each row represents a pair of atoms and their corresponding distance value. The columns represent:
   - Label_1, Label_2: The labels of atoms
   - Distance: The distance between the two atoms

## Terminal Output
When you run the program, you will see the following output in the terminal:

   Number of pairs included: 19968
   epsilon = 0.0100659 eV
   sigma = 3.401 A
   Total Lattice Energy using Lennard-Jones potential = -20.831 eV
   Cohesive Energy per atom = -0.08 eV/at
   Lattice Energy per atom = -0.0813712 eV/at
   Difference = -0.00137118 eV/at
   Percentage Error = 1.71398 %

## Parameters for Ar atoms
   - **Unit Cell**: [1]
      - Structure: FCC
      - lattice constant: 5.31 A       (from page 20)
      - Cohesive Energy:  -4.63 eV/at  (from page 50)

   - **Lennard-Jones Potential Calculation**: [2] (from page 4)
      - epsilon: 116.81 eV/k
      - sigma:   3.401 A

## Reference
   [1] Charles Kittel. Introduction to Solid State Physics. Eighth Edition. Retrieved from ./reference Introduction_to_Solid_State_Physics.pdf
   [2] Seung-Kyo Oh. Modified Lennard-Jones Potentials with a Reduced Temperature-Correction Parameter for Calculating Thermodynamic and Transport Properties: Noble Gases and Their Mixtures (He, Ne, Ar, Kr, and Xe). Hindawi Publishing Corporation Journal of Thermodynamics, Volume 2013. Retrieved from ./reference/Introduction_to_Solid_State_Physics.pdf



