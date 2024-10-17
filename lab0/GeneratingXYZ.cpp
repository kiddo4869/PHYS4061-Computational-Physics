// Write a function that writes a XYZ file from a global array of positions
// Submit both the source code and the file

// Requirement:
// .xyz file is properly formatted, able to be displayed in VMD
// The atomic species and coordinates can be hardcoded or randomly generated. 

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

using namespace std;

int n_frame=3;
int n_atom=5;
double x, y, z;
double MAX=(double)RAND_MAX;

int main()
{
    ofstream MyFile("test.xyz");
	
    for (int i=0; i<n_frame; i++){
	MyFile << n_atom << "\n\n";
	for (int j=0; j<n_atom; j++){
	    x=rand()/MAX;
	    y=rand()/MAX;
	    z=rand()/MAX;
	    MyFile << "Si " << x  << " " << y << " " << z << "\n";
	}
    }

    MyFile.close();

    return 0;
}
