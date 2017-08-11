#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

struct Element
{
    int l, c;
    double v;
}typedef Element;

void Usage (char *argv[])
{
    cout << "=================================================================" << endl;
    cout << "Usage:> " << argv[0] << " <matrix_file> <rhs_file> <method>" << endl;
    cout << "<method> = Number of the method to solve the linear system" << endl;
    cout << "\t1 = Full Pivot LU Decomposition" << endl;
    cout << "\t2 = QR Decomposition with Column pivoting" << endl;
    cout << "\t3 = QR Decomposition using Householder" << endl;
    cout << "=================================================================" << endl;
}

MatrixXf readMatrix (const char mName[])
{
    ifstream in(mName);
    int n, m, k;
    string str;
    getline(in,str);
    in >> n >> m >> k;

    MatrixXf mat(n,m);

    int p, q;
    double v;
    for (int i = 0; i < k; i++)
    {   
        in >> p >> q >> v;
        p--; q--;
        mat(p,q) = v;
    } 
    in.close();
    return mat;
}

VectorXf readRHS (const char rhsName[])
{
    ifstream in(rhsName);
    int n, k;
    string str;
    getline(in,str);
    in >> n >> k;

    VectorXf rhs(n);

    double v;
    for (int i = 0; i < n; i++)
    {   
        in >> v;
        rhs(i) = v;
    } 
    in.close();
    return rhs;
}

int main (int argc, char *argv[])
{
    if (argc-1 < 3)
    {
        Usage(argv);
        exit(EXIT_FAILURE);
    }
    else
    {
        MatrixXf A = readMatrix(argv[1]);
        cout << "Matrix A" << endl;
        cout << A << endl;

        VectorXf b = readRHS(argv[2]);
        cout << "Vector b" << endl;
        cout << b << endl;

        VectorXf x;
        int method = atoi(argv[3]);
        switch (method)
        {
            case 1: {
                        x = A.fullPivLu().solve(b);
                        break;
                    }
                    
            case 2: {
                        x = A.colPivHouseholderQr().solve(b);
                        break;
                    }
            case 3: {
                        x = A.householderQr().solve(b);
                        break;
                    }
        }
        cout << "Solution x" << endl;
        cout << x  << endl;

        double error = (A*x - b).norm() / b.norm();
        cout << "Relative error of the system = " << error << endl;
    }
  
}
