#include "navierstokes.h"

NavierStokes::NavierStokes()
{
      maxTime = 1000;
      dT = 0.01;
      Re = 100; // Reynolds number
      nu = 1 / Re;

      nx = 5;
      ny = 5;
      iMin = 1;
      iMax = nx - 1;
      jMin = 1; // jMin has to be 1;
      jMax = ny - 1;

      Lx = 100;
      Ly = 100;

      createMesh();

      L = zeros(nx * ny, nx * ny);
      L = calcL(L);
      detOfL = determinantOfMatrix(L, L.size());
      for (int i=0; i < L.size(); i++) {
            for (int j=0; j < L.size(); j++) {
                  std::cout << L.at(i).at(j) << " ";
            }
            std::cout << std::endl;
      }
      // invL = getInverse(L);
      invL = invertMatrix(L, L.size());
      std::cout << "past l \n";
} 



NavierStokes::~NavierStokes()
{
}

std::vector<std::vector<float>> NavierStokes::calcL(std::vector<std::vector<float>> L)
{
      for (int j = 1; j < ny + 1; j++)
      {
            for (int i = 0; i < nx; i++)
            {
                  L.at(i + (j - 1) * nx).at(i + (j - 1) * nx) = 2 * pow(dxi, 2) + 2 * pow(dyi, 2);

                  for (int k = i - 1; k < i + 2; k += 2)
                  {
                        if (k > 0 && k <= nx)
                        {
                              L.at(i + (j - 1) * nx).at(i + (j - 1) * nx) = -pow(dxi, 2);
                        }
                        else
                        {
                              L.at(i + (j - 1) * nx).at(i + (j - 1) * nx) = L.at(i + (j - 1) * nx).at(i + (j - 1) * nx) - pow(dxi, 2);
                        }
                        std::cout << " k";
                  }
                  std::cout << " i";

                  for (int k = j - 1; k < j + 2; k += 2)
                  {
                        if (k > 0 && k <= ny)
                        {
                              L.at(i + (j - 1) * nx).at(i + (k - 1) * nx) = -pow(dyi, 2);
                        }
                        else
                        {
                              L.at(i + (j - 1) * nx).at(i + (j - 1) * nx) = L.at(i + (j - 1) * nx).at(i + (j - 1) * nx) - pow(dyi, 2);
                        }
                  }
            }
            for (int i = 0; i < L.at(0).size(); i++)
            {
                  L.at(0).at(i) = 0;
            }
            L.at(0).at(0) = 1;
      }
      return L;
}

void NavierStokes::plot()
{
      std::vector<int> meshXCoordinates, meshYCoordinates;
      std::vector<int> xVelComponent, yVelComponent;

      for (int i = iMin; i < iMax; i++)
      {
            for (int j = jMin; j < jMax; i++)
            {
                  meshXCoordinates.push_back(j);
                  meshYCoordinates.push_back(i);
                  xVelComponent.push_back(u.at(j).at(i));
                  yVelComponent.push_back(v.at(j).at(i));
            }       
      }
      mat::quiver(meshXCoordinates, meshYCoordinates, xVelComponent, yVelComponent);
      mat::show();
}


void NavierStokes::generateVectors()
{
      for (int i = iMin - 1; i < iMax + 2; i++)
      {
            std::vector<float> helper;
            for (int j = jMin - 1; j < jMax + 2; j++)
            {
                  helper.push_back(0);
            }
            u.push_back(helper);
      }

      for (int i = iMin - 1; i < iMax + 2; i++)
      {
            std::vector<float> helper;
            for (int j = jMin - 1; j < jMax + 2; j++)
            {
                  helper.push_back(0);
            }
            v.push_back(helper);
      }

      for (int i = iMin - 1; i < iMax + 2; i++)
      {
            std::vector<float> helper;
            for (int j = jMin - 1; j < jMax + 2; j++)
            {
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

      for (int i = iMin; i < iMax + 1; i++)
      {
            boundaryMinU.push_back(0);
            boundaryMaxU.push_back(0);
      }

      for (int i = jMin; i < jMax + 1; i++)
      {
            boundaryMinU.push_back(0);
            boundaryMaxU.push_back(0);
      }

      int n = 0;
      for (int j = jMin; j < jMax + 1; j++)
      {
            for (int i = iMin; i < iMax + 1; i++)
            {
                  n++;
                  R.push_back(0);
            }
      }

      n = 0;
      for (int j = jMin; j < jMax + 1; j++)
      {
            for (int i = iMin; i < iMax + 1; i++)
            {
                  n++;
                  pv.push_back(0);
            }
      }

      // for (int j=jMin; j < jMax+1; j++) {
      //       std::vector<float> helper;
      //       for (int i=iMin+1; i < iMax+1; i++) {
      //             helper.push_back(0);
      //       }

      // }

      // for (int j=jMin+1; j < jMax+1; j++) {
      //       for (int i=iMin; i < iMax+1; i++) {

      //       }
      // }
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

      generateVectors();
}

void NavierStokes::boundaryConditions()
{
      for (int i = iMin; i < iMax + 1; i++)
      {
            u.at(i).at(jMin - 1) = u.at(i).at(jMin) - 2 * (u.at(i).at(jMin) - uBottom);
            u.at(i).at(jMax + 1) = u.at(i).at(jMax) - 2 * (u.at(i).at(jMax) - uBottom);
            // boundaryMinU.at(i) = u.at(i).at(jMin) - 2 * (u.at(i).at(jMin) - uBottom);
            // boundaryMaxU.at(i) = u.at(i).at(jMax) - 2 * (u.at(i).at(jMax) - uTop);
      }

      for (int j = jMin; j < jMax + 1; j++)
      {
            u.at(iMin - 1).at(j) = u.at(iMin).at(j) - 2 * (u.at(iMin).at(j) - uBottom);
            u.at(iMax + 1).at(j) = u.at(iMax).at(j) - 2 * (u.at(iMax).at(j) - uBottom);
            // boundaryMinV.at(i) = u.at(iMin).at(i) - 2 * (u.at(iMin).at(i) - uBottom);
            // boundaryMaxV.at(i) = u.at(iMax).at(i) - 2 * (u.at(iMax).at(i) - uTop);
      }
}



void NavierStokes::calc()
{
      float time = 0;

      while (time < maxTime)
      {
            time += dT;
            boundaryConditions();

            // for (int i=jMin; i < jMax; i++) {
            //       for (int j=iMin+1; j < iMax; j++) {

            //       }
            // }

            for (int j = jMin; j < jMax + 1; j++)
            {
                  for (int i = iMin + 1; i < iMax + 1; i++)
                  {
                        std::cout << "ja 1" << std::endl;
                        float vHere = 0.25 * (v.at(i - 1).at(j) + v.at(i - 1).at(j + 1) + v.at(i).at(j) + v.at(i).at(j + 1));
                        float a = (nu * (u.at(i - 1).at(j) - 2 * u.at(i).at(j) + u.at(i + 1).at(j)) * pow(dxi, 2));
                        float b = nu * (u.at(i).at(j - 1) - 2 * u.at(i).at(j) + u.at(i).at(j + 1) * pow(dyi, 2));
                        float c = -u.at(i).at(j) * (u.at(i + 1).at(j) - u.at(i - 1).at(j)) * 0.5 * dxi;
                        float d = -vHere * (u.at(i).at(j + 1) - u.at(i).at(j - 1)) * 0.5 * dyi;
                        us.at(i).at(j) = u.at(i).at(j) + dT * (a + b + c + d); // nieuwe s over tijd
                  }
            }
            std::cout << "ja 1 2" << std::endl;


            for (int j = jMin + 1; j < jMax + 1; j++)
            {
                  for (int i = iMin; i < iMax + 1; i++)
                  {
                        std::cout << "ja 2" << std::endl;
                        float uHere = 0.25 * (u.at(i).at(j - 1) + u.at(i).at(j) + u.at(i + 1).at(j - 1) + u.at(i + 1).at(j));
                        float a = (nu * (v.at(i - 1).at(j) - 2 * v.at(i).at(j) + v.at(i + 1).at(j)) * pow(dxi, 2));
                        float b = nu * (v.at(i).at(j - 1) - 2 * v.at(i).at(j) + v.at(i).at(j + 1) * pow(dyi, 2));
                        float c = -uHere * (v.at(i + 1).at(j) - v.at(i - 1).at(j)) * 0.5 * dyi;
                        float d = -v.at(i).at(j) * (v.at(1).at(j + 1) - v.at(i).at(j - 1)) * 0.5 * dxi;
                        vs.at(i).at(j) = v.at(i).at(j) + dT * (a + b + c + d);
                  }
            }
            std::cout << "ja 2 2" << std::endl;

            // std::vector<std::vector<float>> pv = L/R;

            // for (int j=1; j < ny; j++) {
            //       for (int i=1; i < nx; i++) {

            //       }
            // }
            // float dP = ;

            // float a = (p.at(i+1).at(j) + p.at(i-1).at(j) + p.at(i).at(j+1) + p.at(i).at(j-1));
            // float b = -((rho * dx) / 16) * ((2/dT) * (u.at(i+1).at(j) - u.at(i-1).at(j) + v.at(i).at(j+1) - v.at(i).at(j-1));
            // float c = -(2/dx) * (u.at(i).at(j+1) - u.at(i).at(j-1)) * (v.at(i+1).at(j) - v.at(i-1).at(j));
            // float d = (u.at(i+1).at(j)
            // p.at(i).at(j) = (1/4) * a -((rho * dx) / 16) * ((2/dT) * (u.at(i+1).at(j) - u.at(i-1).at(j) + v.at(i).at(j+1) - v.at(i).at(j-1)) -(2/dx) * (u.at(i).at(j+1) - u.at(i).at(j-1)) * (v.at(i+1).at(j) - v.at(i-1).at(j)) - (pow(u.at(i+1).at(j) - u.at(i-1).at(j), 2) / dx) - (pow(v.at(i).at(j+1) - v.at(i).at(j-1), 2) / dx);
            // float changeInU =

            int n = 0;

            for (int j = jMin; j < jMax + 1; j++)
            {
                  for (int i = iMin; i < iMax + 1; i++)
                  {
                        n++;
                        R.at(n) = -rho / dT * ((us.at(i + 1).at(j) - us.at(i).at(j)) * dxi + (vs.at(i).at(j + 1) - vs.at(i).at(j)) * dyi);
                  }
            }

            for (int i = 0; i < L.size(); i++)
            {
                  float lTimesR = 0;
                  for (int j = 0; j < L.size(); j++)
                  {
                        lTimesR += L.at(j).at(i) * R.at(j);
                  }
                  pv.at(i) = lTimesR;
            }

            n = 0;
            std::vector<std::vector<float>> p = zeros(iMax, jMax);

            for (int j = jMin; j < jMax + 1; j++)
            {
                  for (int i = iMin; i < iMax + 1; i++)
                  {
                        n++;
                        p.at(i).at(j) = pv.at(n);
                  }
            }
      }

      // corrector step

      for (int j = jMin; j < jMax + 1; j++)
      {
            for (int i = iMin; i < iMax + 1; i++)
            {
                  u.at(i).at(j) = us.at(i).at(j) - dT / rho * (p.at(i).at(j) - p.at(i).at(j - 1)) * dxi;
            }
      }

      for (int j = jMin + 1; j < jMax + 1; j++)
      {
            for (int i = iMin; i < iMax + 1; i++)
            {
                  v.at(i).at(j) = vs.at(i).at(j) - dT / rho * (p.at(i).at(j) - p.at(i).at(j - 1)) * dxi;
            }
      }

      for (int i = 0; i < v.size(); i++)
      {
            for (int j = 0; j < v.at(i).size(); j++)
            {
                  std::cout << v.at(i).at(j);
            }
            std::cout << std::endl;
      }
      plot();
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

