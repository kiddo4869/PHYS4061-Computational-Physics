# Neighbor List Generator

This program generates a neighbor list of atoms for different crystal structures within a specified cutoff distance.

## Features

- Correct the code from Lab1
- Supports structures of Simple Cubic (SC), Body-Centered Cubic (BCC), Face-Centered Cubic (FCC), and Diamond
- Computes and outputs the reciprocal lattice vectors for each structure
- Outputs a list of nearest neighbors for each structure to a text file (.txt)

## Requirements

- C++ compiler
- Standard C++ libraries

## Usage

1. Compile the program:
   g++ Interatomic_Potential.cpp -o Interatomic_Potential

2. Run the program:
   ./Interatomic_Potential

3. Input the lattice constant and the number of periods in the x, y, and z directions when prompted:
   e.g.
   Please input the lattice constant: 1
   Please input the number of period in x y z direction: 4 4 4

4. The program will generate the following files:
   - Neighbor_List_SC.txt
   - Neighbor_List_BCC.txt
   - Neighbor_List_FCC.txt
   - Neighbor_List_Diamond.txt

## Outputs
The output files will contain a table with the following format:

e.g.
	Label_1,Label_2,Distance
	1,0,1
        2,1,1
        3,0,1
        3,2,1
        4,0,1
        5,1,1
        5,4,1
        6,2,1
        6,5,1

Each column separated by a comma (,) in the text file

Each row represents a pair of atoms and their corresponding distance value. The columns represent:
   - Label_1, Label_2: The labels of atoms
   - Distance: The distance between the two atoms



