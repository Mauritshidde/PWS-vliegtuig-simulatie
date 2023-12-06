#include <vector>
#include <raylib.h>
#include <math.h>

class NavierStokes
{
private:
      int nx; // amount of cells in x direction
      int ny; // amount of cells in y direction
      int iMin; // minimaal 1
      int iMax;
      int jMin;
      int jMax;
      int Lx;
      int Ly;

      float uBottom, uTop, vLeft, vRight;
      float dx;
      float dy;
      float dxi;
      float dyi;
      float nu; 
      float Re;

      std::vector<float> boundaryMinU, boundaryMaxU, boundaryMinV, boundaryMaxV;
      std::vector<std::vector<float>> v, u, p;

public:
      NavierStokes(/* args */);
      ~NavierStokes();

      void boundaryConditions();
      void calc();
      void createMesh();
};

NavierStokes::NavierStokes()
{
      Re = 100; // Reynolds number
      nu = 1 / Re;
      nx = 50;
      ny = 50;
      iMin = 0;
      iMax = nx - 1;
      jMin = 1; // jMin has to be 1;
      jMax = ny - 1;
      Lx = 100;
      Ly = 100;
}

NavierStokes::~NavierStokes()
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

void NavierStokes::createMesh()
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

      for (int i = iMin-1; i < iMax+2; i++)
      {
            std::vector<float> helper;
            for (int j=jMin-1; j < jMax+2; j++) {
                  helper.push_back(0);
            }
            u.push_back(helper);
      }

      for (int i = iMin-1; i < iMax+2; i++)
      {
            std::vector<float> helper;
            for (int j=jMin-1; j < jMax+2; j++) {
                  helper.push_back(0);
            }
            v.push_back(helper);
      }

      for (int i = iMin-1; i < iMax+2; i++)
      {
            std::vector<float> helper;
            for (int j=jMin-1; j < jMax+2; j++) {
                  helper.push_back(0);
            }
            p.push_back(helper);
      }

      // for (int i=0; i < xm.size(); i++) {
      //       std::vector<float> helper;
      //       for (int j=0; j < ym.size(); j++) {
      //             helper.push_back(0);
      //       }
      //       p.push_back(helper);
      // }

      // for (int i=0; i < x.size(); i++) {
      //       std::vector<float> helper;
      //       for (int j=0; j < y.size(); j++) {
      //             helper.push_back(0);
      //       }
      //       v.push_back(helper);
      // }

      // for (int i=0; i < x.size(); i++) {
      //       std::vector<float> helper;
      //       for (int j=0; j < y.size(); j++) {
      //             helper.push_back(0);
      //       }
      //       u.push_back(helper);
      // }

      for (int i=iMin; i < iMax+1; i++) {
            boundaryMinU.push_back(0);
            boundaryMaxU.push_back(0);
      }

      for (int i=jMin; i < jMax+1; i++) {
            boundaryMinU.push_back(0);
            boundaryMaxU.push_back(0);
      }
}


void NavierStokes::boundaryConditions() {
      for (int i=iMin; i < iMax+1; i++) {
            u.at(i).at(jMin-1) = u.at(i).at(jMin) - 2 * (u.at(i).at(jMin) - uBottom);
            u.at(i).at(jMax+1) = u.at(i).at(jMax) - 2 * (u.at(i).at(jMax) - uBottom);
            // boundaryMinU.at(i) = u.at(i).at(jMin) - 2 * (u.at(i).at(jMin) - uBottom);
            // boundaryMaxU.at(i) = u.at(i).at(jMax) - 2 * (u.at(i).at(jMax) - uTop);
      }

      for (int j=jMin; j < jMax+1; j++) {
            u.at(iMin-1).at(j) = u.at(iMin).at(j) - 2 * (u.at(iMin).at(j) - uBottom);
            u.at(iMax+1).at(j) = u.at(iMax).at(j) - 2 * (u.at(iMax).at(j) - uBottom);
            // boundaryMinV.at(i) = u.at(iMin).at(i) - 2 * (u.at(iMin).at(i) - uBottom);
            // boundaryMaxV.at(i) = u.at(iMax).at(i) - 2 * (u.at(iMax).at(i) - uTop);
      }
}

void NavierStokes::calc() {
      float time = 0;
      float maxTime = 1000;
      float dT = 0.01;
      
      while(time < maxTime) {
            time += dT;
            boundaryConditions();

            // for (int i=jMin; i < jMax; i++) {
            //       for (int j=iMin+1; j < iMax; j++) {

            //       }
            // }

            for (int j=jMin; j < jMax+1; j++) {
                  for (int i=iMin+1; i < iMax+1; j++) {
                        float vHere = 0.25 * (v.at(i-1).at(j) + v.at(i-1).at(j+1) + v.at(i).at(j) + v.at(i).at(j+1));
                        float a = (nu * (u.at(i-1).at(j) -2 * u.at(i).at(j) + u.at(i+1).at(j)) * pow(dxi, 2));
                        float b = nu * (u.at(i).at(j-1) -2 * u.at(i).at(j) + u.at(i).at(j+1) * pow(dyi, 2));
                        float c = -u.at(i).at(j) * (u.at(i+1).at(j) - u.at(i-1).at(j)) * 0.5 * dxi;
                        float d = -vHere * (u.at(i).at(j+1) - u.at(i).at(j-1)) * 0.5 * dyi;
                        float us = u.at(i).at(j) + dT * (a + b + c + d);

                  }
            }

            for (int j=jMin+1; j < jMax+1; j++) {
                  for (int i=iMin; i < iMax+1; i++) {
                        float uHere = 0.25 * (u.at(i).at(j-1) + u.at(i).at(j) + u.at(i+1).at(j-1) + u.at(i+1).at(j));
                        float a = (nu * (v.at(i-1).at(j) -2 * v.at(i).at(j) + v.at(i+1).at(j)) * pow(dxi, 2));
                        float b = nu * (v.at(i).at(j-1) -2 * v.at(i).at(j) + v.at(i).at(j+1) * pow(dyi, 2));
                        float c = -uHere * (v.at(i+1).at(j) - v.at(i-1).at(j)) * 0.5 * dyi;
                        float d = -v.at(i).at(j) * (v.at(1).at(j+1) - v.at(i).at(j-1)) * 0.5 * dxi;
                        float vs = v.at(i).at(j) + dT * (a + b + c + d);
                  }
            }

            // float changeInU = 
      }
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