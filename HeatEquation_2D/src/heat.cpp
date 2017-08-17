#include "../include/heat.h"

double initialCondition (double x, double y)
{
    if (x <= 0) return 4;
    else        return 0;
}

HeatSolver2D::HeatSolver2D (int argc, char *argv[])
{
    alfa = atof(argv[1]);
    L = atof(argv[2]);
    dx = atof(argv[3]);
    H = atof(argv[4]);
    dy = atof(argv[5]);
    tmax = atof(argv[6]);
    dt = atof(argv[7]);
    ic = initialCondition;
    nx = nearbyint(L/dx);
    ny = nearbyint(H/dy);
    nt = nearbyint(tmax/dt);
}

void HeatSolver2D::solve ()
{
    SpMat A((nx+1)*(ny+1),(nx+1)*(ny+1));
    buildMatrix(A);
    //cout << A << endl;
    
    VectorXd x((nx+1)*(ny+1));
    setInitialCondition(x);
    //cout << x << endl;

    VectorXd b((nx+1)*(ny+1));
    b = x;

    SparseLU<SpMat> sparseSolver(A);

    double c = 1 / dt;
    // Time loop
    for (int k = 0; k <= nt; k++)
    {
        #ifdef VTK
        if (k % PLOT_STEP == 0) writeVTK(k,x);
        #endif

        // Calculate RHS and set the boundary conditions
        buildRHS(b,c);

        // Solve the current timestep
        x = sparseSolver.solve(b);
        
        // Advance to next timestep
        b = x;
    }
}

void HeatSolver2D::buildRHS (VectorXd &b, double c)
{
    b *= c;

    // Set north BC
    for (int i = 1; i < nx; i++) b(i) = bc_north;
    // Set south BC
    for (int i = 1; i < nx; i++) b(nx*(ny+1)+i) = bc_south;
    // Set east BC
    for (int i = 1; i < nx; i++) b(i*(ny+1)+nx) = bc_east;
    // Set west BC
    for (int i = 1; i < nx; i++) b(i*(ny+1)) = bc_west;
    // Set northwest BC
    b(0) = bc_north + bc_west;
    // Set northeast BC
    b(nx) = bc_north + bc_east;
    // Set southwest BC
    b(nx*(ny+1)) = bc_south + bc_west;
    // Set southeast BC
    b(nx*(ny+1)+nx) = bc_south + bc_east;

}

void HeatSolver2D::setInitialCondition (VectorXd &b)
{
    for (int i = 0; i <= nx; i++)
    {
        double y = i*dy;
        for (int j = 0; j <= ny; j++)
        {
            int k = i * (ny+1) + j;
            double x = j*dx;
            b(k) = ic(x,y);
        }
        
    }
}

void HeatSolver2D::computeCoeff (double &a, double &b, double &c)
{
    a = (1.0/dt)+(4*alfa/(dx*dx));
    b = -alfa/(dx*dx);
    c = 1.0/dt;
}

void HeatSolver2D::insertBC_Coefficient (vector<T> &coeff)
{
    insertCoefficient(0,0,1,coeff);                             // Upper left
    insertCoefficient(nx,nx,1,coeff);                           // Upper right
    insertCoefficient(nx*(nx+1),nx*(nx+1),1,coeff);             // Bottom left
    insertCoefficient(nx*(nx+1)+nx,nx*(nx+1)+nx,1,coeff);       // Bottom right
    for (int i = 1; i < nx; i++)
    {
        int north = i;
        int south = (ny + 1)*nx + i;
        int east = i*(nx + 1) + nx;
        int west = i*(nx+1);
        insertCoefficient(north,north,1,coeff);
        insertCoefficient(south,south,1,coeff);
        insertCoefficient(east,east,1,coeff);
        insertCoefficient(west,west,1,coeff);
    }
}

void HeatSolver2D::insertCoefficient (int i, int j, double val, vector<T> &coeff)
{
    int n = (nx + 1)*(ny + 1);
    // Make sure that the index is not out of range
    if (j >= 0 && j <= n && i >= 0 && i <= n)
        coeff.push_back(T(i,j,val));
        
}

void HeatSolver2D::buildMatrix (SpMat &A)
{
    vector<T> coeff;
    double a, b, c;
    computeCoeff(a,b,c);
    insertBC_Coefficient(coeff);                    // Boundary condition coefficients
    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            int center = i*(ny+1)+j;
            int north = (i-1)*(ny+1)+j;
            int south = (i+1)*(ny+1)+j;
            int east = i*(ny+1)+(j+1);
            int west = i*(ny+1)+(j-1);
            //cout << "Node " << center << endl;
            //cout << "North = " << north << " South = " << south << " West = " << west << " East = " << east << endl;
            insertCoefficient(center,center,a,coeff);
            insertCoefficient(center,north,b,coeff);
            insertCoefficient(center,south,b,coeff);
            insertCoefficient(center,east,b,coeff);
            insertCoefficient(center,west,b,coeff);
        }
        
    }
    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();
}

void HeatSolver2D::printHeatSolver2D ()
{
    cout << "Alfa = " << alfa << endl;
    cout << "L = " << L << endl;
    cout << "dx = " << dx << endl;
    cout << "tmax = " << tmax << endl;
    cout << "dt = " << dt << endl;
}

ostream& operator<<(ostream& ost, const HeatSolver2D& h)
{
    ost << "nt = " << h.nt << endl;
    ost << "nx = " << h.nx << endl;
    ost << "ny = " << h.ny << endl;
    ost << "Alfa = " << h.alfa << endl;
    ost << "L = " << h.L << endl;
    ost << "dx = " << h.dx << endl;
    ost << "H = " << h.H << endl;
    ost << "dy = " << h.dy << endl;
    ost << "tmax = " << h.tmax << endl;
    ost << "dt = " << h.dt << endl;
    ost << "bc_west = " << bc_west << endl;
    ost << "bc_east = " << bc_east << endl;
    ost << "bc_north = " << bc_north << endl;
    ost << "bc_south = " << bc_south << endl;
    return ost;
}

void HeatSolver2D::writeSolution (int iter, VectorXd x)
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

void HeatSolver2D::writeVTK (int k, const VectorXd x)
{
    ostringstream ss;
    ss << "VTK/solution" << k << ".vtk";
    string filename(ss.str());
    ofstream out(filename.c_str());
    out << "# vtk DataFile Version 3.0" << endl;
    out << "Heat Equation 2D" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << (ny+1)*(nx+1) << " float" << endl;
    for (int i = 0; i <= ny; i++)
    {
        double y = i*dy;
        for (int j = 0; j <= nx; j++)
        {
            double x = j*dx;
            out << x << " " << y << " 0" << endl;
        }
    }
    out << "VERTICES " << (ny+1)*(nx+1) << " " << (ny+1)*(nx+1)*2 << endl;
    for (int i = 0; i < (ny+1)*(nx+1); i++)
        out << "1 " << i << endl;
    out << "POINT_DATA " << (ny+1)*(nx+1) << endl;
    out << "SCALARS solution float 1" << endl;
    out << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i <= ny; i++)
        for (int j = 0; j <= nx; j++)
            out << x(j*ny+i) << endl;
    out.close();
}

void Usage (char *argv[])
{
    cout << "=================================================================" << endl;
    cout << "Usage:> " << argv[0] << " <ALPHA> <L> <dx> <H> <dy> <tmax> <dt>" << endl;
    cout << "<ALPHA> = Diffusion coefficient" << endl;
    cout << "<L> = Length of the plate" << endl;
    cout << "<dx> = Size of the space discretization on the x direction" << endl;
    cout << "<H> = Height of the plate" << endl;
    cout << "<dy> = Size of the space discretization on the y direction" << endl;
    cout << "<tmax> = Maximum simulation time" << endl;
    cout << "<dt> = Size of the time discretization" << endl;
    cout << "=================================================================" << endl;
}