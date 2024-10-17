# Lattice Structure Generator

This program generates atomic lattice structures in ".xyz" format.

The structures is in below:
- Simple Cubic (SC)
- Body-Centered Cubic (BCC)
- Face-Centered Cubic (FCC)
- Diamond

## Features

- Generate lattice structures based on user input for lattice constant and dimensions.
- Output atomic coordinates in `.xyz` format.

## Requirements

- C++ compiler

## Usage

1. Compile the program:
   g++ Generate_Structure.cpp -o Generate_Structure

2. Run the program:
   ./Generate_Structure

3. Input the lattice constant and the number of periods in the x, y, and z directions when prompted:
   e.g.
   Please input the lattice constant: 1
   Please input the number of period in x y z direction: 3 3 3   

4. The program will generate the following files:
   - SC.xyz
   - BCC.xyz
   - FCC.xyz
   - Diamond.xyz