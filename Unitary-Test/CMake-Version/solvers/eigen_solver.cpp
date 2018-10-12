#include "eigen_solver.h"

VectorXd read_rhs (const char rhsName[])
{
    int n, k;
    string str;
    
    ifstream in(rhsName);
    in >> n;
    VectorXd b(n);
    for (int i = 0; i < n; i++)
        in >> b(i);
    in.close();

    return b;
}

SpMat read_matrix (const char mName[])
{
    vector<T> coeff;
    int n, m;
    int l, c;
    double val;
    string str;

    ifstream in(mName);
    in >> n >> m;
    while (in >> l >> c >> val)
    {
        coeff.push_back(T(l,c,val));
    }
    in.close();

    // Build the sparse matrix by using the coefficient vector
    SpMat A(n,m);
    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();

    return A;
}