#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<cstdlib>

using namespace std;

//////////// GLOBAL VARIABLES ////////////

double a;                            // lattice constant
int n_x, n_y, n_z;                   // number of periods in x, y, and z directions

double x_0=0.0, y_0=0.0, z_0=0.0;    // initial positions in x, y, and z directions
int n_atoms;                         // number of atoms
vector<double> x, y, z;              // vectors storing coordinates of atoms in x, y, and z directions

//////////// FUNCTIONS ////////////

// Functions to get atoms positions for different structures
void WriteAtomicPositions(vector<vector<double>> b, string structure_type)
{
    // resize the coordinate arrays
    x.resize(n_atoms);
    y.resize(n_atoms);
    z.resize(n_atoms);

    ofstream MyFile(structure_type + ".xyz");
    cout << "Generating " << structure_type << ".xyz file" << "\n";

    // write the number of atoms
    MyFile << n_atoms << "\n\n";

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

                    MyFile << structure_type << " " << x[n] << " " << y[n] << " " << z[n] << "\n";

                    n++;
                }       
            }
        }
    }
}

void WriteAtomicPositions_SC()
{
    // write the number of atoms
    n_atoms=n_x*n_y*n_z;

    // create the basis according to the structure
    vector<vector<double>> b = {
        {0.0, 0.0, 0.0}
    };

    // write atoms' positions
    WriteAtomicPositions(b, "SC");
}

void WriteAtomicPositions_BCC()
{
    // write the number of atoms
    n_atoms=n_x*n_y*n_z*2;

    // create the basis according to the structure
    vector<vector<double>> b = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    };

    // write atoms' positions
    WriteAtomicPositions(b, "BCC");
}

void WriteAtomicPositions_FCC()
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

    // write atoms' positions
    WriteAtomicPositions(b, "FCC");
}

void WriteAtomicPositions_Diamond()
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

    // write atoms' positions
    WriteAtomicPositions(b, "Diamond");
}

int main()
{
    // read lattice constant a and number of period n_x, n_y, and n_z
    cout << "Please input the lattice constant: ";
    cin >> a;
    cout << "Please input the number of period in x y z direction: ";
    cin >> n_x >> n_y >> n_z;
    cout << n_x << " " << n_y << " " << n_z << "\n";

    // generate structures
    WriteAtomicPositions_SC();
    WriteAtomicPositions_BCC();
    WriteAtomicPositions_FCC();
    WriteAtomicPositions_Diamond();

    cout << "The SC, BCC, FCC, and Diamond crystalline structure has been output!\n";
        
    return 0;
}


