#include <vector>

class navierstokes
{
private:
      int nx; // amount of cells in x direction
      int ny; // amount of cells in y direction
      int iMin;
      int iMax;
      int jMin;
      int jMax;
      int Lx;
      int Ly;

      float dx;
      float dy;
      float dxi;
      float dyi;
public:
      navierstokes(/* args */);
      ~navierstokes();

      void createMesh();
};

navierstokes::navierstokes(/* args */)
{
      nx = 50;
      ny = 50;
      iMin = 0;
      iMax = nx - 1;
      jMin = 0;
      jMax = ny - 1;
      Lx = iMax + 1;
      Ly = jMax + 1;
}

navierstokes::~navierstokes()
{
}

std::vector<float> linspace(int startX, int endX, int steps)
{
      float stepSize = (endX - startX) / (steps - 1);
      std::vector<float> coords;

      for (int i = 0; i < steps; i++)
      {
            coords.push_back(startX + (stepSize * i));
      }

      return coords;
}

void navierstokes::createMesh()
{
      float stepSizeX = (Lx - 0) / (nx);
      float stepSizeY = (Ly - 0) / (ny);

      std::vector<float> x = linspace(0, Lx, nx + 1);
      std::vector<float> y = linspace(0, Ly, ny + 1);
      std::vector<float> xm = linspace(0 + stepSizeX, Lx - stepSizeX, nx);
      std::vector<float> ym = linspace(0 + stepSizeY, Ly - stepSizeY, ny);

      dx = x.at(iMin + 1) - x.at(iMin);
      dy = y.at(jMin + 1) - y.at(jMin);
      dxi = 1 / dx;
      dyi = 1 / dy;
}

/*
Code overview:

Set input parameters: viscosity, density, number of grid points, time information, and boundary conditions
• Create the index extents and the computational grid (see Section 2)
• Initialize any arrays you use to allocate the memory
• Create the Laplacian operator (see Section 6)
• Loop over time (use a for or while loop)
      – Update time t = t + ∆t
      – Apply boundary conditions to the velocity field (see Section 8)
      – Perform the predictor step to find u∗ and v∗ (see Sections 4 and 5)
      – Form the right-hand-side of the Poisson equation (see Section 6)
      – Solve for the pressure using pv = L\R and convert the pressure vector pv into a matrix p(i, j)
      (see Section 6)
      – Perform the corrector step to find un+1 and vn+1 (see Section 7)
      – Plot the velocity field and the pressure field
• End Simulation

We gebruiken classes dus de variables zet je in de class.

/* ∇ · u = 0
∇· is the divergence of the fluid, which is 0 for an incompressible fluid
u = [u,v], the velocity vector (2d)

(∂u/∂t) + u · ∇u = − (1/ρ) ∇p + ν∇²u

t = time
ρ = density = air density
p = pressure = air pressure op vlieghoogte (lookuptable?)
v =  kinematic viscosity (misschien onnodig)
*/