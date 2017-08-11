#include <iostream>
#include <cassert>
#include <vector>
#include <QImage>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Triplet<double> T;          // Triplet (row,col,value)

// Insert the coefficient in the proper way of the problem
// If the index 'i' or 'j' of the cell is on the boundary we know the solution, so we reduce the value on the RHS
// If the index 'i' and 'j' are not in the boundary then is unknown coefficient, so we add a Triplet to make the SparseMatrix
void insertCoefficient(int id, int i, int j, double w, std::vector<T>& coeffs,
                       Eigen::VectorXd& b, const Eigen::VectorXd& boundary)
{
  int n = int(boundary.size());
  int id1 = i+j*n;
  if(i==-1 || i==n) b(id) -= w * boundary(j);       // constrained coefficient
  else  if(j==-1 || j==n) b(id) -= w * boundary(i); // constrained coefficient
  else  coeffs.push_back(T(id,id1,w));              // unknown coefficient
}

void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n)
{
  b.setZero();
  Eigen::ArrayXd boundary = Eigen::ArrayXd::LinSpaced(n, 0,M_PI).sin().pow(2);
  for(int j = 0; j < n; ++j)
  {
    for(int i = 0; i < n; ++i)
    {
      int id = i+j*n;
      insertCoefficient(id, i-1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i+1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j-1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j+1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j,    4, coefficients, b, boundary);
    }
  }
}

void saveAsBitmap(const VectorXd &x, int n, const char *filename)
{
  // transform the solution, which is in [0,1] to values in [0,255].
  Array<unsigned char,Eigen::Dynamic,Eigen::Dynamic> bits = (x*255).cast<unsigned char>();
  QImage img(bits.data(), n,n,QImage::Format_Indexed8);
  img.setColorCount(256);
  for(int i=0;i<256;i++) img.setColor(i,qRgb(i,i,i));
  img.save(filename);
}

void Usage (const char pName[])
{
    cout << "====================================" << endl;
    cout << "Usage:> " << pName << " <name_image.jpg>" << endl;
    cout << "Example: " << pName << " test.jpg" << endl;
    cout << "====================================" << endl;
}

int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        Usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    else
    {
        int n = 300;  // size of the image
        int m = n*n;  // number of unknows (=number of pixels)
        // Assembly:
        vector<T> coefficients;            // list of non-zeros coefficients
        VectorXd b(m);                     // the right hand side-vector resulting from the constraints
        buildProblem(coefficients, b, n);
        
        // Build the sparse matrix using the non-zero coefficients
        SpMat A(m,m);
        A.setFromTriplets(coefficients.begin(), coefficients.end());

        // Solving:
        SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
        VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side
        //cout << b << endl;
        // Export the result to a file:
        saveAsBitmap(x, n, argv[1]);
    }
    return 0;
}