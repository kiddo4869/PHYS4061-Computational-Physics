//Write two functions that accept two 3D vectors and calculate their dot product and cross product, respectively, by using the technique of passing references to store the results.
//Attention: Dot product gives a scalar, and cross product gives a vector.
//Tips:
//For Cross product, you can pass 3 arrays, a, b, and c, by references and store a x b in  c 
//For Dot product, it is relatively simple,

#include<iostream>
#include<vector>

using namespace std;

vector<double> v1 = {1, 0, 0};
vector<double> v2 = {0, 1, 0};
vector<double> cp(3);
double dp;

vector<double> CrossProduct(vector<double> v1, vector<double> v2)
{
    vector<double> v3(3);
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return v3;
}

double DotProduct(vector<double> v1, vector<double> v2)
{
    double sum=0;
    for (int i=0; i<3; i++){
        sum+= v1[i] * v2[i];
    }
    return sum;
}


int main()
{
    dp = DotProduct(v1, v2);
    cp = CrossProduct(v1, v2);

    cout << "Dot Product: " << dp << "\n";
    cout << "Cross Product: " << cp[0] << " " << cp[1] << " " << cp[2] << "\n";
    return 0;
}


