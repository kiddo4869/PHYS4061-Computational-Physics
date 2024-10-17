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

double x_0=0.0, y_0=0.0, z_0=0.0;    // initial positions in x, y, and z directions
int n;                               // counter of atoms
int n_atoms;
double x, y;
double MAX=(double)RAND_MAX;

vector<double> a1(3), a2(3), a3(3);  // unit cell lattice vectors
vector<double> b1(3), b2(3), b3(3);  // reciprocal lattice vectors

double s0 = 0;
double ds = 0;

//////////// FUNCTIONS ////////////

// Functions to integrate
double fofx(double x) {
    return pow(x, 8);
}

void TrapezoidalRule(int N)
{
    double a=-1, b=1;
    double dx = abs(b - a) / N;
    double x1, x2;
    y = 0;

    for (int i=0; i<N; i++)
    {
        x1 = a + dx * i;
        x2 = a + dx * (i + 1);
        y += (fofx(x1) + fofx(x1)) * dx / 2;
    }
}

void RandomSampling(int N)
{
    s0 = 0;
    ds = 0;
    for (int i=0; i<N; i++)
    {
        x = rand()/MAX;
        y = fofx(x);
        s0 += y;
        ds += y * y;
    }

    s0 /= N;
    s0 *= 2;
    //ds /= N;
    //ds = sqrt(abs(ds - pow(s0, 2)) / N);
}


int main()
{
    int N;
    cout << "Please input the no. of slices (or no. of sample points): ";
    cin >> N;

    double I_exact = 2.0/9.0;
    
    cout << "Trapezoidal Rule\n";
    TrapezoidalRule(N);
    cout << "N = " << N << "\n";
    cout << "I = " << y << "\n";
    cout << "squared deviation = " << pow(y - I_exact, 2) << "\n";

    cout << "---------------------\n";

    cout << "Random Sampling\n";
    RandomSampling(N);
    cout << "N = " << N << "\n";
    cout << "I = " << s0 << "\n";
    cout << "squared deviation = " << pow(s0 - I_exact, 2) << "\n";
            
    return 0;
}


