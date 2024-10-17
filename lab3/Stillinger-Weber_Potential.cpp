#define _USE_MATH_DEFINES
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<cstdlib>
#include<tuple>
#include<cmath>

using namespace std;

//////////// GLOBAL VARIABLES ////////////

double a;                            // lattice constant
int n_x, n_y, n_z;                   // number of periods in x, y, and z directions
int n_atoms;
vector<double> x, y, z;              // vectors storing coordinates of atoms in x, y, and z directions

// Newly added for lab2
vector<double> a1(3), a2(3), a3(3);  // unit cell lattice vectors
vector<double> b1(3), b2(3), b3(3);  // reciprocal lattice vectors

vector<double> coor_1(3);            // (x1, y1, z1) displacement, cartesian coordinates before applying PBC
vector<double> coor_2(3);            // (x2, y2, z2) displacement, cartesian coordinates after applying PBC
vector<double> coor_frac(3);         // (n1, n2, n3) displacement, fractional coordinates

double dist;                         // distance of two atoms
double cut_off;                      // distance cut off
double cut_off_delta;                // added to cut off

double epsilon;
double sigma;
vector<int> label_1_list, label_2_list;
vector<double> dist_list;

//////////// FUNCTIONS ////////////

// Functions to get atoms positions for different structures (from lab1)
void GetAtomicPositions(vector<vector<double>> b)
{
    // resize the coordinate arrays
    x.resize(n_atoms);
    y.resize(n_atoms);
    z.resize(n_atoms);

    // initial positions in x, y, and z directions
    double x_0=0.0, y_0=0.0, z_0=0.0;

    // initialize atom counter
    int n = 0;

    // loop over each period in x, y, and z directions
    for (int i=0; i<n_x; i++)
    {
        for (int j=0; j<n_y; j++)
        {
            for (int k=0; k<n_z; k++)
            {
                // base position
                double base_x = x_0 + i * a;
                double base_y = y_0 + j * a;
                double base_z = z_0 + k * a;

                // calculate the atom positions
                for (const auto& b_i : b) {
                    x[n] = base_x + b_i[0] * a;
                    y[n] = base_y + b_i[1] * a;
                    z[n] = base_z + b_i[2] * a;
                    n++;
                }       
            }
        }
    }
}

void GetAtomicPositions_SC()
{
    // write the number of atoms
    n_atoms=n_x*n_y*n_z;

    // create the basis according to the structure
    vector<vector<double>> b = {
        {0.0, 0.0, 0.0}
    };

    // get atoms' positions
    GetAtomicPositions(b);
}

void GetAtomicPositions_BCC()
{
    // write the number of atoms
    n_atoms=n_x*n_y*n_z*2;

    // create the basis according to the structure
    vector<vector<double>> b = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    };

    // get atoms' positions
    GetAtomicPositions(b);
}

void GetAtomicPositions_FCC()
{
    // write the number of atoms
    n_atoms=n_x*n_y*n_z*4;

    // create the basis according to the structure
    vector<vector<double>> b = {
        {0.0, 0.0, 0.0},
        {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5},
        {0.5, 0.5, 0.0}
    };

    // get atoms' positions
    GetAtomicPositions(b);
}

void GetAtomicPositions_Diamond()
{
    // write the number of atoms
    n_atoms=n_x*n_y*n_z*8;

    // create the basis according to the structure
    vector<vector<double>> b = {
        {0.0, 0.0, 0.0},
        {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5},
        {0.5, 0.5, 0.0},
        {0.25, 0.25, 0.25},
        {0.25, 0.75, 0.75},
        {0.75, 0.25, 0.75},
        {0.75, 0.75, 0.25}
    };

    // get atoms' positions
    GetAtomicPositions(b);
}

// Functions to calculate cross and dot products (from lab0)
vector<double> cross(vector<double> v1, vector<double> v2)
{
    vector<double> v3(3);
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return v3;
}

double dot(vector<double> v1, vector<double> v2)
{
    double sum=0;
    for (int i=0; i<3; i++){
        sum+= v1[i] * v2[i];
    }
    return sum;
}

// Newly added for lab2
// Function to calculate reciprocal lattice vectors b1, b2, b3 from unit cell lattice vectors a1, a2, a3
void CalculateReciprocalLatticeVectors()
{
    vector<double> h1 = cross(a2, a3);   // h1 = a2 × a3
    double d1 = dot(a1, h1);             // d1 = a1 · (a2 × a3)
    for (int i=0; i<h1.size(); i++)
    { 
        b1[i] = (h1[i] / d1);            // b1 = a2 × a3 / [a1 · (a2 × a3)]
    } 

    vector<double> h2 = cross(a3, a1);   // h2 = a3 × a1
    double d2 = dot(a2, h2);             // d2 = a2 · (a3 × a1)
    for (int i=0; i<h2.size(); i++)
    {
        b2[i] = (h2[i] / d2);            // b2 = a3 × a1 / [a2 · (a3 × a1)]
    }
    
    vector<double> h3 = cross(a1, a2);   // h3 = a1 × a2
    double d3 = dot(a3, h3);             // d3 = a3 · (a1 × a2)
    for (int i=0; i<h3.size(); i++)
    {
        b3[i] = (h3[i] / d3);            // b3 = a1 × a2 / [a3 · (a1 × a2)]
    }
}

// Function to apply periodic boundary condition (PBC) on displacement coordinates
void ApplyPeriodicBoundaryCondition()
{
    // calculate fractional coordinates from cartesian coordinates
    coor_frac[0] = dot(b1, coor_1);
    coor_frac[1] = dot(b2, coor_1);
    coor_frac[2] = dot(b3, coor_1);
    
    // Apply PBC on fractional coordinates
    for (int i=0; i<coor_frac.size(); i++) {
        coor_frac[i] -= round(coor_frac[i]);
    }

    // convert fractional coordinates back to cartesian coordinates
    coor_2[0] = dot(a1, coor_frac);
    coor_2[1] = dot(a2, coor_frac);
    coor_2[2] = dot(a3, coor_frac);
}

// Function to write neighbor list for different structures
void WriteNeighborlist(string structure_type)
{
    // Clean the vector container
    label_1_list.clear();
    label_2_list.clear();
    dist_list.clear();

    string filename = "Neighbor_List_";
    filename += structure_type;
    filename += "_OverCounted";
    filename += ".txt";
    ofstream MyFile(filename);
    
    // write the column names
    MyFile << "Label_1,Label_2,Distance\n";

    // loop for every atom pairs
    for (int idx_i=0; idx_i<n_atoms; idx_i++)
    {
        for (int idx_j=0; idx_j<idx_i; idx_j++)
        {   
            // calculate a displacement coordinate of pair
            coor_1[0] = x[idx_j] - x[idx_i];
            coor_1[1] = y[idx_j] - y[idx_i];
            coor_1[2] = z[idx_j] - z[idx_1];
                
            // apply periodic boundary condition (PBC)
            ApplyPeriodicBoundaryCondition();

            // calculate the distance of the pair after PBC
            dist = sqrt(pow(coor_2[0], 2) + pow(coor_2[1], 2) + pow(coor_2[2], 2));
                
            // include the pair if the distance is less than cut-off
            if (dist < (cut_off + cut_off_delta))
            {
                MyFile << idx_i << "," << idx_j << "," << dist <<"\n";
                label_1_list.push_back(idx_i);
                label_2_list.push_back(idx_j);
                dist_list.push_back(dist);
            } 
        }
    }
    cout << "Number of pairs included: " << dist_list.size() << "\n";
    MyFile.close();
}

void WriteNeighborlist_OverCounted(string structure_type)
{
    // Clean the vector container
    label_1_list.clear();
    label_2_list.clear();
    dist_list.clear();

    string filename = "Neighbor_List_";
    filename += structure_type;
    filename += "_OverCounted";
    filename += ".txt";
    ofstream MyFile(filename);
    
    // write the column names
    MyFile << "Label_1,Label_2,Distance\n";

    // loop for every atom pairs
    for (int idx_i=0; idx_i<n_atoms; idx_i++)
    {
        for (int idx_j=0; idx_j<n_atoms; idx_j++)
        {   
            if (idx_i != idx_j)
            {
                // calculate a displacement coordinate of pair
                coor_1[0] = x[idx_j] - x[idx_i];
                coor_1[1] = y[idx_j] - y[idx_i];
                coor_1[2] = z[idx_j] - z[idx_1];
                
                // apply periodic boundary condition (PBC)
                ApplyPeriodicBoundaryCondition();

                // calculate the distance of the pair after PBC
                dist = sqrt(pow(coor_2[0], 2) + pow(coor_2[1], 2) + pow(coor_2[2], 2));
                
                // include the pair if the distance is less than cut-off
                if (dist < (cut_off + cut_off_delta))
                {
                    MyFile << idx_i << "," << idx_j << "," << dist <<"\n";
                    label_1_list.push_back(idx_i);
                    label_2_list.push_back(idx_j);
                    dist_list.push_back(dist);
                }
            } 
        }
    }
    cout << "Number of pairs included: " << dist_list.size() << "\n";
    MyFile.close();
}

// Function to pack the sub-functions for different structures
void GenerateNeighborlist(string structure_type)
{
    // get the atomic positions in the unit cell
    if (structure_type == "SC"){
        GetAtomicPositions_SC();}
    else if (structure_type == "BCC"){
        GetAtomicPositions_BCC();}
    else if (structure_type == "FCC"){
        GetAtomicPositions_FCC();}
    else if (structure_type == "Diamond"){
        GetAtomicPositions_Diamond();}
    else{
        throw invalid_argument("Invalid structure type");}

    // calculate the reciprocal lattice vectors
    cout << "calculating the reciprocal lattice vectors...\n";
    CalculateReciprocalLatticeVectors();

    // write neighbor list
    cout << "writing neighbor list... with cut off = " << cut_off << " nm\n";
    WriteNeighborlist(structure_type);
}

void GenerateNeighborlist_OverCounted(string structure_type)
{
    // get the atomic positions in the unit cell
    if (structure_type == "SC"){
        GetAtomicPositions_SC();}
    else if (structure_type == "BCC"){
        GetAtomicPositions_BCC();}
    else if (structure_type == "FCC"){
        GetAtomicPositions_FCC();}
    else if (structure_type == "Diamond"){
        GetAtomicPositions_Diamond();}
    else{
        throw invalid_argument("Invalid structure type");}

    // calculate the reciprocal lattice vectors
    cout << "calculating the reciprocal lattice vectors...\n";
    CalculateReciprocalLatticeVectors();

    // write neighbor list
    cout << "writing neighbor list... with cut off = " << cut_off << " nm\n";
    WriteNeighborlist_OverCounted(structure_type);
}

// Newly added for lab3
double A, B, a;
double V_2(double r)
{
    double V = (A * B * pow(r, -4) - A) * exp(1/(r - a));
    return V;
}

double get_len(vector<double> vec)
{
    return pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2);
}

double theta(vector<double> v_ij, vector<double> v_ik, double r_ij, double r_ik)
{
    return acos(dot(v_ij, v_ik) / (r_ij * r_ik));
}

double V_3(double r_ij, double r_ik, double theta_jik)
{
    double gamma_ij;
    double gamma_ik;
    double lambda;
    double V = lambda * exp(gamma_ij / (r_ij - a) + gamma_ik / (r_ik - a)) * pow(cos(theta_jik) + double(1)/double(3), 2);
    return V;
}

double Stillinger_Weber_Potential(double r_ij, double r_ik, double theta_jik)
{
    double V = V_2(r_ij) + V_2(r_ij, r_ik, theta_jik);
    return V;
}

double Total_Energy()
{
    // Initialize the energy E
    double E=0.0;

    // Calculate the Lennard-Jones potential according to the neighbor list
    int n = dist_list.size();
    for (int i=0; i<n; i++)
    {
        E += Stillinger_Weber_Potential(dist_list[i]);
    }

    return E;
}


int main()
{   
    // read number of period n_x, n_y, and n_z
    cout << "Please input the number of period in x y z direction: ";
    cin >> n_x >> n_y >> n_z;

    cout << "\n----------C atoms----------\n";
    // set structure type and lattice constant for Ni atoms
    string structure_type = "FCC";
    a = 0.35238;               // (unit: nm)

    // set unit cell lattice vectors
    a1 = {1.0*a*n_x, 0.0*a*n_y, 0.0*a*n_z};
    a2 = {0.0*a*n_x, 1.0*a*n_y, 0.0*a*n_z};
    a3 = {0.0*a*n_x, 0.0*a*n_y, 1.0*a*n_z};

    // set potential parameters
    epsilon = 0.66092;         // (unit: eV/at)
    sigma = 0.223949;          // (unit: nm)
    cut_off = 0.85 * a;        // (unit: nm)
    cut_off_delta = 0.0;
    
    // calculate the LennardJones potential
    GenerateNeighborlist(structure_type);
    double E = Total_Energy();

    cout << "Total Lattice Energy using Lennard-Jones potential = " << E << " eV\n";
    cout << "Lattice Energy per atom = " << double(E / n_atoms) << " eV/at\n";
    cout << "Lattice Energy per atom = " << double(E / n_atoms) * 96.5 << " kJ/mol\n";
        
    return 0;
}


