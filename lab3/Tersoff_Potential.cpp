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

vector<double> coor_ij(3);           // (x1, y1, z1) displacement, cartesian coordinates before applying PBC
vector<double> coor_ik(3);           // (x2, y2, z2) displacement, cartesian coordinates after applying PBC
vector<double> coor_frac(3);         // (n1, n2, n3) displacement, fractional coordinates

double dist, r_ij, r_ik;             // distance of two atoms
double cut_off;                      // distance cut off
double cut_off_delta;                // added to cut off

double epsilon;
double sigma;
vector<int> label_1_list, label_2_list;
vector<double> dist_list;

// Terstoff potential parameters
double A_i, B_i, lmd_i, mu_i, R_i, S_i, beta_i, n_i, c_i, d_i, h_i;
double A_ij, B_ij, lmd_ij, mu_ij, R_ij, S_ij;
double b_ij, X_ij, xi_ij, omega_ik;
double g_ijk, theta_ijk;
double E, E_coh;

//////////// FUNCTIONS ////////////

// Functions to get atoms positions for different structures (from lab1)
void GetAtomicPositions(vector<vector<double>> b) {
    x.resize(n_atoms);                // resize the coordinate arrays
    y.resize(n_atoms);
    z.resize(n_atoms);
    double x_0=0.0, y_0=0.0, z_0=0.0; // initial positions in x, y, and z directions
    int n = 0;                        // initialize atom counter

    // loop over each period in x, y, and z directions
    for (int i=0; i<n_x; i++) {
        for (int j=0; j<n_y; j++) {
            for (int k=0; k<n_z; k++) {
                double base_x = x_0 + i * a;     // base position
                double base_y = y_0 + j * a;
                double base_z = z_0 + k * a;

                for (const auto& b_i : b) {      // calculate the atom positions
                    x[n] = base_x + b_i[0] * a;
                    y[n] = base_y + b_i[1] * a;
                    z[n] = base_z + b_i[2] * a;
                    n++;
                }       
            }
        }
    }
}

void GetAtomicPositionsByStructure(string structure_type) {
    vector<vector<double>> b;
    if (structure_type == "SC") {
        n_atoms=n_x*n_y*n_z;           // write the number of atoms
        b = {                          // create the basis according to the structure
            {0.0, 0.0, 0.0}
        };
    }
    else if (structure_type == "BCC") {
        n_atoms=n_x*n_y*n_z*2;           // write the number of atoms
        b = {                            // create the basis according to the structure
            {0.0, 0.0, 0.0},
            {0.5, 0.5, 0.5}
        };
    }
    else if (structure_type == "FCC") {
        n_atoms=n_x*n_y*n_z*4;           // write the number of atoms
        b = {                            // create the basis according to the structure
            {0.0, 0.0, 0.0},
            {0.0, 0.5, 0.5},
            {0.5, 0.0, 0.5},
            {0.5, 0.5, 0.0}
        };
    }
    else if (structure_type == "Diamond") {
        n_atoms=n_x*n_y*n_z*8;           // write the number of atoms
        b = {                            // create the basis according to the structure
            {0.0, 0.0, 0.0},
            {0.0, 0.5, 0.5},
            {0.5, 0.0, 0.5},
            {0.5, 0.5, 0.0},
            {0.25, 0.25, 0.25},
            {0.25, 0.75, 0.75},
            {0.75, 0.25, 0.75},
            {0.75, 0.75, 0.25}
        };
    }
    else {
        throw invalid_argument("Invalid structure type");
    }

    GetAtomicPositions(b);        // get atoms' positions
}

// Functions to calculate cross and dot products (from lab0)
vector<double> cross(vector<double> v1, vector<double> v2) {
    vector<double> v3(3);
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return v3;
}

double dot(vector<double> v1, vector<double> v2) {
    double sum=0.0;
    for (int i=0; i<3; i++) {
        sum+= v1[i] * v2[i];
    }
    return sum;
}

// Newly added for lab2
// Function to calculate reciprocal lattice vectors b1, b2, b3 from unit cell lattice vectors a1, a2, a3
void CalculateReciprocalLatticeVectors() {
    vector<double> h1 = cross(a2, a3);   // h1 = a2 x a3
    double d1 = dot(a1, h1);             // d1 = a1 • (a2 x a3)
    for (int i=0; i<h1.size(); i++) { 
        b1[i] = (h1[i] / d1);            // b1 = a2 x a3 / [a1 • (a2 x a3)]
    } 

    vector<double> h2 = cross(a3, a1);   // h2 = a3 x a1
    double d2 = dot(a2, h2);             // d2 = a2 • (a3 x a1)
    for (int i=0; i<h2.size(); i++) {
        b2[i] = (h2[i] / d2);            // b2 = a3 x a1 / [a2 • (a3 x a1)]
    }
    
    vector<double> h3 = cross(a1, a2);   // h3 = a1 x a2
    double d3 = dot(a3, h3);             // d3 = a3 • (a1 x a2)
    for (int i=0; i<h3.size(); i++) {
        b3[i] = (h3[i] / d3);            // b3 = a1 x a2 / [a3 • (a1 x a2)]
    }
}

// Function to apply periodic boundary condition (PBC) on displacement coordinates
vector<double> get_disp(int idx_1, int idx_2) {
    // from idx_1 to idx_2
    return {x[idx_2] - x[idx_1], y[idx_2] - y[idx_1], z[idx_2] - z[idx_1]};
}

vector<double> ApplyPeriodicBoundaryCondition(vector<double> coor_1) {
    // calculate fractional coordinates from cartesian coordinates
    coor_frac[0] = dot(b1, coor_1);
    coor_frac[1] = dot(b2, coor_1);
    coor_frac[2] = dot(b3, coor_1);
    
    // Apply PBC on fractional coordinates
    for (int i=0; i<coor_frac.size(); i++) {
        if (coor_frac[i] != 0.5) {
            coor_frac[i] -= round(coor_frac[i]);
        }
    }

    // convert fractional coordinates back to cartesian coordinates
    return {dot(a1, coor_frac), dot(a2, coor_frac), dot(a3, coor_frac)};
}

// Newly added for lab3: Eenrgy calculation
double get_len(vector<double> vec) {
    return sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
}

double theta(vector<double> v_ij, vector<double> v_ik, double r_ij, double r_ik) {
    return acos(dot(v_ij, v_ik) / (r_ij * r_ik));
}

double f_C(double r_ij) {
    if (r_ij <= R_ij) {
        return 1;
    }
    else if (r_ij >= S_ij) {
        return 0;
    }
    else { // R_ij <= r_ij <= S_ij
        return 0.5 + 0.5 * cos(M_PI * (r_ij - R_ij) / (S_ij - R_ij));
    }
}

double f_R(double r_ij) {
    cout << "A_ij = " << A_ij << ", lmd_ij = " << lmd_ij << ", r_ij = " << r_ij << "\n";
    return A_ij * exp(-lmd_ij * r_ij);
}

double f_A(double r_ij) {
    //cout << "B_ij = " << B_ij << ", mu_ij = " << mu_ij << ", r_ij = " << r_ij << "\n";
    return - B_ij * exp(-mu_ij * r_ij);
}

double g(double theta_ijk) {
    return 1 + pow(c_i, 2) / pow(d_i, 2) - pow(c_i, 2) / (pow(d_i, 2) + pow(h_i - cos(theta_ijk), 2));
}

double Tersoff_Potential(double r_ij) {
    double V_ij = f_C(r_ij) * (f_R(r_ij) + b_ij * f_A(r_ij));
    // cout << "f_C = " << f_C(r_ij) << ", f_R = " << f_R(r_ij) << ", f_A = " << f_A(r_ij) << ", b_ij = " << b_ij << ", V_ij = " << V_ij << "\n";
    return V_ij;
}

// Function to write neighbor list for different structures
void WriteNeighborlist_OverCounted(string structure_type) {
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

    // Initialize the energy E
    E = 0.0;

    // loop for every atom pairs
    for (int idx_i=0; idx_i<n_atoms; idx_i++) {
        for (int idx_j=0; idx_j<n_atoms; idx_j++) {   
            if (idx_i != idx_j) {
                
                coor_ij = get_disp(idx_i, idx_j);                   // calculate a displacement coordinate of pair
                coor_ij = ApplyPeriodicBoundaryCondition(coor_ij);  // apply periodic boundary condition (PBC)
                r_ij = get_len(coor_ij);                            // calculate the distance of the pair after PBC
                
                if (r_ij < (cut_off + cut_off_delta)) {             // include the pair if the distance is less than cut-off
                    xi_ij = 0.0;
                    for (int idx_k=0; idx_k<n_atoms; idx_k++) {
                        if (idx_k != idx_i && idx_k != idx_j) {

                            coor_ik = get_disp(idx_i, idx_k);
                            coor_ik = ApplyPeriodicBoundaryCondition(coor_ik);
                            r_ik = get_len(coor_ik);

                            if (r_ik < (cut_off + cut_off_delta)) {
                                theta_ijk = theta(coor_ij, coor_ik, r_ij, r_ik);
                                g_ijk = g(theta_ijk);
                                xi_ij += f_C(r_ik) * omega_ik * g_ijk;
                            }
                        }
                    }

                    b_ij = X_ij * pow(1 + pow(beta_i * xi_ij, n_i), -0.5 / n_i);
                    E += Tersoff_Potential(r_ij);

                    MyFile << idx_i << "," << idx_j << "," << r_ij <<"\n";
                    label_1_list.push_back(idx_i);
                    label_2_list.push_back(idx_j);
                    dist_list.push_back(r_ij);
                }
            } 
        }
    }

    // Divide the total energy by 2 to avoid double counting
    E /= 2;

    cout << "Number of pairs included: " << dist_list.size() << "\n";
    MyFile.close();
}

// Function to pack the sub-functions for different structures
void GenerateNeighborlist(string structure_type) {
    // get the atomic positions in the unit cell
    GetAtomicPositionsByStructure(structure_type);

    // calculate the reciprocal lattice vectors
    cout << "calculating the reciprocal lattice vectors...\n";
    CalculateReciprocalLatticeVectors();

    // write neighbor list
    cout << "writing neighbor list... with cut off = " << cut_off << " nm\n";
    WriteNeighborlist_OverCounted(structure_type);
}

// Simplied version of GetParameters function
void GetParameters_Single(string element) {
    if (element == "C") {
        A_i = 1.3936 * pow(10, 3);  // (unit: eV)
        B_i = 3.467 * pow(10, 2);   // (unit: eV)
        lmd_i = 3.4879;             // (unit: A)
        mu_i = 2.2119;              // (unit: A)
        beta_i = 1.5724 * pow(10, -7);
        n_i = 7.2751 * pow(10, -1);
        c_i = 3.8049 * pow(10, 4);
        d_i = 4.384;
        h_i = -5.7058 * pow(10, -1);
    }
    else if (element == "Si") {
        A_i = 1.8308 * pow(10, 3);  // (unit: eV)
        B_i = 4.7118 * pow(10, 2);  // (unit: eV)
        lmd_i = 2.4799;             // (unit: A)
        mu_i = 1.7322;              // (unit: A)
        beta_i = 1.1000 * pow(10, -6);
        n_i = 7.8734 * pow(10, -1);
        c_i = 1.0039 * pow(10, 5);
        d_i = 1.6217 * pow(10, 1);
        h_i = -5.9825 * pow(10, -1);
    }
    else {
        throw invalid_argument("Invalid element");
    }

    // For single element
    A_ij = sqrt(A_i * A_i);
    B_ij = sqrt(B_i * B_i);
    R_ij = sqrt(R_i * R_i);
    S_ij = sqrt(S_i * S_i);
    lmd_ij = 0.5 * (lmd_i + lmd_i);
    mu_ij = 0.5 * (mu_i + mu_i);
    omega_ik = 1.0;
    X_ij = 1.0;
}

int main(){   
    // read number of period n_x, n_y, and n_z
    cout << "Please input the number of period in x y z direction: ";
    cin >> n_x >> n_y >> n_z;

    /*
    cout << "\n----------C atoms----------\n";
    a = 3.56683;                    // (unit: A)
    
    // set unit cell lattice vectors
    a1 = {1.0*a*n_x, 0.0*a*n_y, 0.0*a*n_z};
    a2 = {0.0*a*n_x, 1.0*a*n_y, 0.0*a*n_z};
    a3 = {0.0*a*n_x, 0.0*a*n_y, 1.0*a*n_z};
    cut_off = n_x * a;              // (unit: A)
    cut_off_delta = 0.1 * a;        // (unit: A)
    R_i = cut_off;                  // (unit: A)
    S_i = cut_off + cut_off_delta;  // (unit: A)
    
    // calculate the Terstoff potential parameters
    GetParameters_Single("C");
    GenerateNeighborlist("Diamond");

    cout << "Lattice constant = " << a << " A\n";
    cout << "First nearest neighbor distance = " << a * sqrt(3) / 4 << " A\n";
    cout << "R_ij = " << R_ij << " A\n";
    cout << "S_ij = " << S_ij << " A\n";

    double E_coh = -7.37;
    cout << "Total Lattice Energy using Tersoff potential = " << E << " eV\n";
    cout << "Cohesive Energy per atom = " << E_coh << " eV/at\n";
    cout << "Lattice Energy per atom = " << double(E / n_atoms) << " eV/at\n";
    cout << "Difference = " << double(E / n_atoms) - E_coh << " eV/at\n";
    cout << "Percentage Error = " << (double(E / n_atoms) - E_coh) / E_coh * 100 << " %\n";
    */

    cout << "\n----------Si atoms----------\n";
    a = 5.43;                       // (unit: A)
    
    // set unit cell lattice vectors
    a1 = {1.0*a*n_x, 0.0*a*n_y, 0.0*a*n_z};
    a2 = {0.0*a*n_x, 1.0*a*n_y, 0.0*a*n_z};
    a3 = {0.0*a*n_x, 0.0*a*n_y, 1.0*a*n_z};
    cut_off = 2.7;              // (unit: A)
    cut_off_delta = 0.3;        // (unit: A)
    R_i = cut_off;                  // (unit: A)
    S_i = cut_off + cut_off_delta;  // (unit: A)
    

    // calculate the Terstoff potential parameters
    GetParameters_Single("Si");
    GenerateNeighborlist("Diamond");

    cout << "Lattice constant = " << a << " A\n";
    cout << "First nearest neighbor distance = " << a * sqrt(3) / 4 << " A\n";
    cout << "R_ij = " << R_ij << " A\n";
    cout << "S_ij = " << S_ij << " A\n";

    E_coh = -4.63;
    cout << "Total Lattice Energy using Tersoff potential = " << E << " eV\n";
    cout << "Cohesive Energy per atom = " << E_coh << " eV/at\n";
    cout << "Lattice Energy per atom = " << double(E / n_atoms) << " eV/at\n";
    cout << "Difference = " << double(E / n_atoms) - E_coh << " eV/at\n";
    cout << "Percentage Error = " << (double(E / n_atoms) - E_coh) / E_coh * 100 << " %\n";

    return 0;
}


