#include "../include/heat.h"

double initialCondition (double x, double L)
{
    return sin(M_PI*x/L);
}

HeatSolver1D::HeatSolver1D (int argc, char *argv[])
{
    alfa = atof(argv[1]);
    L = atof(argv[2]);
    dx = atof(argv[3]);
    tmax = atof(argv[4]);
    dt = atof(argv[5]);
    ic = initialCondition;
    nx = nearbyint(L/dx);
    nt = nearbyint(tmax/dt);
}

void HeatSolver1D::solve ()
{
    SpMat A(nx+1,nx+1);
    buildMatrix(A);
    //cout << A << endl;
    
    VectorXd x(nx+1);
    setInitialCondition(x);
    //cout << x << endl;

    VectorXd b(nx+1);
    b = x;

    SparseLU<SpMat> sparseSolver(A);
    #ifdef VTK
    vector< VectorXd > solution;
    #endif

    double d = 1 / dt;
    // Time loop
    for (int k = 0; k <= nt; k++)
    {
        if (k % PLOT_STEP == 0) writeSolution(k,x);
        #ifdef VTK
        solution.push_back(x);
        #endif

        // Calculate RHS and set the boundary conditions
        b(0) = bc_west;
        b(nx) = bc_east;
        b *= d;

        // Solve the current timestep
        x = sparseSolver.solve(b);
        
        // Advance to next timestep
        b = x;
    }

    #ifdef VTK
    writeVTK(solution);
    #endif
}

void HeatSolver1D::setInitialCondition (VectorXd &b)
{
    for (int i = 0; i <= nx; i++)
    {
        double x = i*dx;
        b(i) = ic(x,L);
    }
}

void HeatSolver1D::computeCoeff (double &a, double &b, double &c, double &d)
{
    a = (1.0/dt)+(2*alfa/(dx*dx));
    b = -alfa/(dx*dx);
    c = b;
    d = 1.0/dt;
}

void HeatSolver1D::insertCoefficient (int i, int j, double val, vector<T> &coeff)
{
    // Make sure that the column index is not out of eange
    if (j >= 0 && j <= nx)
        coeff.push_back(T(i,j,val));
}

void HeatSolver1D::buildMatrix (SpMat &A)
{
    int n = nx + 1;
    vector<T> coeff;
    double a, b, c, d;
    computeCoeff(a,b,c,d);
    insertCoefficient(0,0,1,coeff);                 // Boundary coefficient
    for (int i = 1; i < n-1; i++)
    {
        insertCoefficient(i,i,a,coeff);
        insertCoefficient(i,i-1,b,coeff);
        insertCoefficient(i,i+1,c,coeff);
    }
    insertCoefficient(n-1,n-1,1,coeff);             // Boundary coefficient
    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();
}

void HeatSolver1D::printHeatSolver1D ()
{
    cout << "Alfa = " << alfa << endl;
    cout << "L = " << L << endl;
    cout << "dx = " << dx << endl;
    cout << "tmax = " << tmax << endl;
    cout << "dt = " << dt << endl;
}

ostream& operator<<(ostream& ost, const HeatSolver1D& h)
{
    ost << "nt = " << h.nt << endl;
    ost << "nx = " << h.nx << endl;
    ost << "Alfa = " << h.alfa << endl;
    ost << "L = " << h.L << endl;
    ost << "dx = " << h.dx << endl;
    ost << "tmax = " << h.tmax << endl;
    ost << "dt = " << h.dt << endl;
    ost << "bc_west = " << bc_west << endl;
    ost << "bc_east = " << bc_east << endl;
    return ost;
}

void HeatSolver1D::writeSolution (int iter, VectorXd x)
{
    double px;
    ostringstream ss;
    ss << "Data/data" << iter << ".dat";
    string filename(ss.str());
    ofstream out(filename.c_str());
    for (int i = 0; i <= nx; i++)
    {
        px = i*dx;
        out << px << " " << x(i) << endl;
    }
         
    out.close();
}

void HeatSolver1D::writeVTK (const vector<VectorXd> sol)
{
    ofstream out("VTK/solution.vtk");
    out << "# vtk DataFile Version 3.0" << endl;
    out << "Heat Equation 1D" << endl;
    out << "ASCII" << endl;
    out << "DATASET STRUCTURED_POINTS" << endl;
    out << "DIMENSIONS 2 " << (nx+1) << " " << (nt+1) << endl;
    out << "ORIGIN 0 0 0" << endl;
    out << "SPACING 1 1 1" << endl;
    out << "POINT_DATA " << (nt+1)*(nx+1)*2 << endl;
    out << "SCALARS solution float" << endl;
    out << "LOOKUP_TABLE default" << endl;
    for (int k = 0; k <= nt; k++)
        for (int i = 0; i <= nx; i++)
            out << sol[k][i] << " " << sol[k][i] << endl;
    out.close();
}

void Usage (char *argv[])
{
    cout << "=================================================================" << endl;
    cout << "Usage:> " << argv[0] << " <ALPHA> <L> <dx> <tmax> <dt>" << endl;
    cout << "<ALPHA> = Diffusion coefficient" << endl;
    cout << "<L> = Size of the bar" << endl;
    cout << "<dx> = Size of the space discretization" << endl;
    cout << "<tmax> = Maximum simulation time" << endl;
    cout << "<dt> = Size of the time discretization" << endl;
    cout << "=================================================================" << endl;
}