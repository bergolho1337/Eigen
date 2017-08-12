#include <iostream>
 
#include <Eigen/Sparse>
 
int main(int, char *[])
{
  // Create matrix
  typedef Eigen::SparseMatrix<double> SparseMatrixType;
  SparseMatrixType A(2, 2);
 
  // Create the right-hand-side vector
  Eigen::VectorXd b(2);
 
  // Fill matrix
  A.coeffRef(0, 0) += 1;
  A.coeffRef(1, 1) += 3;
 
  // Fill vector
  b[0] = 1;
  b[1] = 2;
 
  // Solve the (symmetric) system
  Eigen::SimplicialLDLT<SparseMatrixType> sparseSolver(A);
  Eigen::VectorXd x = sparseSolver.solve(b);
  if(sparseSolver.info() != Eigen::Success)
  {
    throw std::runtime_error("Decomposition failed!");
  }
 
  std::cout << "Result: " << x << std::endl;
 
  return EXIT_SUCCESS;
}