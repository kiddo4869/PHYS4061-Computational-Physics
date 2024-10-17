#include<iostream>
#include<vector>
using namespace std;

// Class of using the Faddeev-Leverrier method to obtain the inversion of a matrix
class Faddeev
{
    public:
        int n;                              // Matrix size
        vector<vector<double>> A;           // Input matrix
        vector<vector<double>> D;           // Inverse matrix

        vector<double> C;
        vector<vector<vector<double>>> S;

        Faddeev(vector<vector<double>> input_matrix, int size)
        {
            n = size;
            A.resize(n);
            A = input_matrix;
        }

        void PrintSquareMatrix(vector<vector<double>> SM);
        void PrintInputMatrix();
        void PrintInverseMatrix();
        double TraceOfMatrix(vector<vector<double>> a);
        vector<vector<double>> ProductOfMatrices(vector<vector<double>> a, vector<vector<double>> b);
        vector<vector<double>> FaddeevLeverrierMethod();
};

// Functions to print the square matrix on terminal
void Faddeev::PrintSquareMatrix(vector<vector<double>> SM)
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            cout << SM[i][j] << " ";
        }
            cout << "\n";
    }
};

void Faddeev::PrintInputMatrix()
{
    cout << "\nThe input matrix is\n";
    PrintSquareMatrix(A);
};


void Faddeev::PrintInverseMatrix()
{
    cout << "\nAnd its inverse matrix is\n";
    PrintSquareMatrix(D);
};

// Function to return the trace of matrix
double Faddeev::TraceOfMatrix(vector<vector<double>> a)
{
    double sum=0.0;
    for (int i=0; i<n; i++)
    {
        sum += a[i][i];
    }
    return sum;
};

// Function to return the product of matrices
vector<vector<double>> Faddeev::ProductOfMatrices(vector<vector<double>> a, vector<vector<double>> b)
{
    vector<vector<double>> c(n, vector<double>(n, 0));
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            for (int k=0; k<n; k++)
            {
            c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}

// Function to return the inversion of matrix
vector<vector<double>> Faddeev::FaddeevLeverrierMethod()
{

    C.resize(n);
    S.resize(n, vector<vector<double>>(n, vector<double>(n, 0)));
    D.resize(n, vector<double>(n, 0));

    // Initialize S_0 = I
    for (int i=0; i<n; i++)
    {
        S[0][i][i] = 1;
    }
    
    // Calculate S_k = AS_{k-1} + C_{n-k}I
    for (int k=1; k<n; k++)
    {
        S[k] = ProductOfMatrices(A, S[k-1]);
        C[n-k] = -TraceOfMatrix(S[k]) / k;
        for (int i=0; i<n; i++)
        {
            S[k][i][i] += C[n-k];
        }
    }

    // Obtain the inverse A^{-1} = -\frac{1}{C_0}S_{n-1}
    C[0] = -TraceOfMatrix(ProductOfMatrices(A, S[n-1])) / n;
    for (int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
                D[i][j] = -S[n-1][i][j] / C[0];
        }
    }

    return D;
}

int main()
{
    // Read the matrix size
    int N;
    cout << "Please input the size of square matrix: ";
    cin >> N;
    
    // Read the matrix elements
    vector<vector<double>> Matrix(N, vector<double>(N, 0));
    for (int i=0; i<N; i++)
    {
        cout << "Please input the row " << i + 1 << " of matrix : ";
        for (int j=0; j<N; j++)
        {
            cin >> Matrix[i][j];
        }
    }

    // Create the matrix calculator from Faddeev class
    Faddeev matrix_calculator_1(Matrix, N);
    matrix_calculator_1.PrintInputMatrix();

    // Calculate the inverse of matrix
    vector<vector<double>> I_Matrix;
    I_Matrix = matrix_calculator_1.FaddeevLeverrierMethod();
    matrix_calculator_1.PrintInverseMatrix();

    // Check the inverse of the inverse matrix is the matrix
    Faddeev matrix_calculator_2(I_Matrix, N);
    matrix_calculator_2.PrintInputMatrix();

    vector<vector<double>> II_Matrix;
    II_Matrix = matrix_calculator_2.FaddeevLeverrierMethod();
    matrix_calculator_2.PrintInverseMatrix();

    return 0;
}


