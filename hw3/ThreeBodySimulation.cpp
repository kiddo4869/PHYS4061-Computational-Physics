#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<cmath>

using namespace std;

//////////// GLOBAL VARIABLES ////////////
int n_bodies = 3;                       // Number of bodies
double G = 39.478;                      // Gravitational constant (AU^3 / (Msun * yr^2))
double m1, m2, m3;                      // Masses of three bodies
vector<double> x1(3), x2(3), x3(3);     // Positions of three bodies
vector<double> x1_(3), x2_(3), x3_(3);  // New positions of three bodies
vector<double> v1(3), v2(3), v3(3);     // Velocities of three bodies
vector<double> v1_(3), v2_(3), v3_(3);  // Velocities of three bodies
vector<double> a1(3), a2(3), a3(3);     // Accelerations of three bodies

double t, h, T;                      // Time, time step and total time
int n_steps;                         // Number of time steps
double K, K1, K2, K3, V;             // Total kinetic energy, individual kinetic and potential energies

void initialize(){    
    // Masses (in solar masses)
    m1 = 1.0;      // Sun (1 solar mass)
    m2 = 0.001;    // Earth-like (approximate)
    m3 = 0.001;    // Another planet

    // Positions (in AU)
    x1 = {0.0, 0.0, 0.0};           // Sun at center
    x2 = {1.0, 0.0, 0.0};           // Earth at 1 AU
    x3 = {1.5, 0.0, 0.0};           // Another planet at 1.5 AU

    // Velocities (in AU/year)
    v1 = {0.0, 0.0, 0.0};
    v2 = {0.0, 6.28318, 0.0};       // Approximately Earth's orbital velocity
    v3 = {0.0, 5.13, 0.0};          // Velocity for 1.5 AU orbit

    t = 0.0;

    cout << "\nInitizlizing the three-body system...\n";
    cout << "Gravitational constant: " << G << " AU^3 / (Msun * yr^2)\n";
    cout << "Masses (in solar masses): " << m1 << " " << m2 << " " << m3 << "\n";
    cout << "\nInitial positions of three bodies (in AU): \n";
    cout << "Mass 1: " << x1[0] << " " << x1[1] << " " << x1[2] << "\n";
    cout << "Mass 2: " << x2[0] << " " << x2[1] << " " << x2[2] << "\n";
    cout << "Mass 3: " << x3[0] << " " << x3[1] << " " << x3[2] << "\n";
    cout << "\nInitial velocities of three bodies (in AU/year): \n";
    cout << "Mass 1: " << v1[0] << " " << v1[1] << " " << v1[2] << "\n";
    cout << "Mass 2: " << v2[0] << " " << v2[1] << " " << v2[2] << "\n";
    cout << "Mass 3: " << v3[0] << " " << v3[1] << " " << v3[2] << "\n";
    cout << "\nTotal time: " << T << " years\n";
    cout << "Time step: " << h << " years\n";
}

// Functions to calculate the acceleration of body i due to body j
vector<double> a(vector<double> x_i, vector<double> x_j, double m_j) {
    vector<double> a_ij(3);
    vector<double> x_ij(3);
    double x_ij_abs = 0.0;

    // body 1
    for (int i=0; i<3; i++) {
        x_ij[i] = x_j[i] - x_i[i];
        x_ij_abs += pow(x_ij[i], 2);
    }
    x_ij_abs = sqrt(x_ij_abs);

    for (int i=0; i<3; i++) {
        a_ij[i] = G * (m_j * x_ij[i] / pow(x_ij_abs, 3));
    }
    return a_ij;
}

vector<double> add(vector<double> v1, vector<double> v2) {
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

// Functions to simulate the three-body system
void update_euler(){

    // Calculate the acceleration of three bodies
    a1 = add(a(x1, x2, m2), a(x1, x3, m3));
    a2 = add(a(x2, x1, m1), a(x2, x3, m3));
    a3 = add(a(x3, x1, m1), a(x3, x2, m2));

    // Find the next steps for the velocities and positions of three bodies
    for (int i=0; i<3; i++) {
        v1_[i] = v1[i] + h * a1[i];
        v2_[i] = v2[i] + h * a2[i];
        v3_[i] = v3[i] + h * a3[i];
    }

    for (int i=0; i<3; i++) {
        x1_[i] = x1[i] + h * v1[i];
        x2_[i] = x2[i] + h * v2[i];
        x3_[i] = x3[i] + h * v3[i];
    }

    // Update the positions and velocities of three bodies
    x1 = x1_;
    x2 = x2_;
    x3 = x3_;
    v1 = v1_;
    v2 = v2_;
    v3 = v3_;
}

void self_start_leapfrog(){
    a1 = add(a(x1, x2, m2), a(x1, x3, m3));
    a2 = add(a(x2, x1, m1), a(x2, x3, m3));
    a3 = add(a(x3, x1, m1), a(x3, x2, m2));

    // self start using Euler method
    // from v_0 to v_1/2
    for (int i=0; i<3; i++) {
        v1[i] = v1[i] + 0.5 * h * a1[i];
        v2[i] = v2[i] + 0.5 * h * a2[i];
        v3[i] = v3[i] + 0.5 * h * a3[i];
    }
}

void update_leapfrog(){

    // from x_0 to x_1 (we need v_1/2)
    for (int i=0; i<3; i++) {
        x1_[i] = x1[i] + h * v1[i];
        x2_[i] = x2[i] + h * v2[i];
        x3_[i] = x3[i] + h * v3[i];
    }

    // calculate a_1
    a1 = add(a(x1_, x2_, m2), a(x1_, x3_, m3));
    a2 = add(a(x2_, x1_, m1), a(x2_, x3_, m3));
    a3 = add(a(x3_, x1_, m1), a(x3_, x2_, m2));

    // from v_1/2 to v_3/2 (we need a_1)
    for (int i=0; i<3; i++) {
        v1_[i] = v1[i] + h * a1[i];
        v2_[i] = v2[i] + h * a2[i];
        v3_[i] = v3[i] + h * a3[i];
    }

    // Update the positions and velocities of three bodies
    // x_0 <-- x_1
    x1 = x1_;
    x2 = x2_;
    x3 = x3_;

    // v_1/2 <-- v_3/2
    v1 = v1_;
    v2 = v2_;
    v3 = v3_;
}

void calculate_energy(){
    // Calculate the kinetic energy of three bodies
    K1 = 0.5 * m1 * (pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    K2 = 0.5 * m2 * (pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
    K3 = 0.5 * m3 * (pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));
    K = K1 + K2 + K3;

    // Calculate the potential energy of three bodies
    V = - G * m1 * m2 / sqrt(pow(x2[0] - x1[0], 2) + pow(x2[1] - x1[1], 2) + pow(x2[2] - x1[2], 2))
        - G * m1 * m3 / sqrt(pow(x3[0] - x1[0], 2) + pow(x3[1] - x1[1], 2) + pow(x3[2] - x1[2], 2))
        - G * m2 * m3 / sqrt(pow(x3[0] - x2[0], 2) + pow(x3[1] - x2[1], 2) + pow(x3[2] - x2[2], 2));
}

void WriteTrajectory(string filename, string format="xyz"){
    ofstream Trajectory(filename+"."+format);
    ofstream Energies(filename+"_energies.txt");

    cout << "----------------------------------------\n";
    cout << "Writing trajectory to " << filename+"."+format << "\n";
    cout << "Method: " << filename << "\n";

    // initialize the system
    initialize();
    
    // write the number of bodies
    Energies << "step" << " " << "time" << " " << "Total" << " " << "Potential" << " " << "Kinetic" << " " << "Kinetic_1" << " " << "Kinetic_2" << " " << "Kinetic_3" << "\n";

    n_steps = T / h;
    cout << "Number of steps: " << n_steps << "\n";

    // self start for leapfrog
    if (filename == "leapfrog") {
        self_start_leapfrog();
        t += h;
    }

    for (int n=0; n<n_steps; n++) {
        // write the number of bodies
        Trajectory << n_bodies << "\n\n";

        // write the positions of three bodies
        Trajectory << "Body" << " " << x1[0] << " " << x1[1] << " " << x1[2] << "\n";
        Trajectory << "Body" << " " << x2[0] << " " << x2[1] << " " << x2[2] << "\n";
        Trajectory << "Body" << " " << x3[0] << " " << x3[1] << " " << x3[2] << "\n";

        // calculate the energies of three bodies
        calculate_energy();
        Energies << n << " " << t << " " << K + V << " " << V << " " << K << " " << K1 << " " << K2 << " " << K3 << "\n";

        // update the positions of three bodies
        if (filename == "euler") {
            update_euler();
        }
        else if (filename == "leapfrog") {
            update_leapfrog();
        }
        else {
            throw invalid_argument("Invalid filename");
        }

        t += h;
    }

    Trajectory.close();
    Energies.close();
}

int main()
{
    cout << "Please input the time period (in year): ";
    cin >> T;

    cout << "Please input the time step (in year): ";
    cin >> h;

    WriteTrajectory("euler", "xyz");
    WriteTrajectory("leapfrog", "xyz");

    return 0;
}


