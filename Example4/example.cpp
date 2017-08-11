#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

void Usage (char *argv[])
{
    cout << "=================================================================" << endl;
    cout << "Usage:> " << argv[0] << " <matrix_file> <rhs_file>" << endl;
    cout << "=================================================================" << endl;
}


int main (int argc, char *argv[])
{
    if (argc-1 < 2)
    {
        Usage(argv);
        exit(EXIT_FAILURE);
    }
    else
    {
        vector<T> coeff;
        int n, m, k;
        int l, c;
        double val;
        string str;
        ifstream in(argv[1]);
        getline(in,str);
        in >> n >> m >> k;
        for (int i = 0; i < k; i++)
        {
            in >> l >> c >> val;
            l--; c--;
            coeff.push_back(T(l,c,val));
        }
        in.close();

        SpMat A(n,m);
        A.setFromTriplets(coeff.begin(), coeff.end());
        cout << A << endl;

        in.open(argv[2]);
        getline(in,str);
        in >> n >> k;
        VectorXf b(n);
        for (int i = 0; i < n; i++)
            in >> b(i);
        cout << b << endl;
        in.close();

        SparseLU<SparseMatrix<double> >   solver;
        
        // /solver.analyzePattern(A);
    }
  
}
