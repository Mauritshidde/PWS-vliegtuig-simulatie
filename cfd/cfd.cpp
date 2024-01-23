#include "cfd.h"

#define WITHOUT_NUMPY
#include "Physics/matplotlibcpp.h"

namespace mat = matplotlibcpp;

void Cfd::createMesh()
{
    for (int i = 0; i < nz; i++)
    {
        std::vector<std::vector<MeshCube>> helper;
        std::vector<std::vector<double>> helper2;
        for (int j = 0; j < nx; j++)
        {
            std::vector<MeshCube> helperHelper;
            std::vector<double> helperHelper2;
            for (int k = 0; k < ny; k++)
            {
                helperHelper.push_back(MeshCube);
                helperHelper2.push_back(0);
            }
            helper.push_back(helperHelper);
            helper2.push_back(helperHelper2);
        }
        mesh.push_back(helper);
        divergenceVelocityField.push_back(helper2);
        gradientPressureField.push_back(helper2);
    }
}

void Cfd::setBoundaryConditions(float velocityXDirectionStart, float velocityYDirectionStart, float velocityZDirectionStart, float velocityXDirectionEnd, float velocityYDirectionEnd, float velocityZDirectionEnd)
{
    for (int i = 0; i < nz; i++)
    {
        for (int k = 0; k < ny; k++)
        {
            mesh.at(i).at(1).at(k).boundary = true;
            mesh.at(i).at(1).at(k).velocityX = velocityXDirectionStart;

            mesh.at(i).at(nx - 1).at(k).boundary = true;
            mesh.at(i).at(nx - 1).at(k).velocityX = velocityXDirectionEnd;
        }
    }

    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            mesh.at(i).at(j).at(1).boundary = true;
            mesh.at(i).at(j).at(1).velocityY = velocityYDirectionStart;

            mesh.at(i).at(j).at(ny - 1).boundary = true;
            mesh.at(i).at(j).at(ny - 1).velocityY = velocityYDirectionEnd;
        }
    }

    for (int j = 0; j < nx; j++)
    {
        for (int k = 0; k < ny; k++)
        {
            mesh.at(1).at(j).at(k).boundary = true;
            mesh.at(1).at(j).at(k).velocityZ = velocityZDirectionStart;

            mesh.at(nz - 1).at(j).at(k).boundary = true;
            mesh.at(nz - 1).at(j).at(k).velocityZ = velocityZDirectionEnd;
        }
    }
}

void Cfd::setPlaneBoundary()
{
}

void Cfd::iterativeSolver(double density, float densityX, float densityY)
{
    float newDensity = 0; // Dn(x,y)
    float x = 0;      // Dn(x+1,y)
    float y = densityX;      // Dn(x-1,y)
    float z = 0;      // Dn(x,y+1)
    float w = densityY;      // Dn(x,y-1)

    float previousNewDensity = 0; // previousDn value of (x,y) in itterative solver
    float previousx = 0;          // previousDn value of (x+1,y) in itterative solver
    float previousz = 0;          // previousDn value of (x,y+1) in itterative solver

    while (newDensity != previousNewDensity || x != previousx || z != previousz)
    {
        previousNewDensity = newDensity;
        previousx = x;
        previousz = z;

        newDensity = (density + k * ((x + y + w + z) / 4)) / (1 + k);
        x = 4 * newDensity + ((4 * newDensity + density) / 4) - (y + z + w);
        z = 4 * newDensity + ((4 * newDensity + density) / 4) - (x + y + w);
    }
}

void Cfd::calc()
{
    float tijd = 0;
    while (tijd < maxTime)
    {
        tijd += dT;
        
        // for (int i=1; i < nz-1; i++) {
            for (int j = 1; j < nx - 1; j++)
            {
                for (int k = 1; k < ny - 1; k++)
                {
                    mesh.at(1).at(j).at(k).density = iterativeSolver(mesh.at(1).at(j).at(k));
                    // mesh.medianSurroundingDensity = (mesh.at(1).at(j+1).at(k) + mesh.at(1).at(j-1).at(k) + mesh.at(1).at(j).at(k+1) + mesh.at(1).at(j).at(k-1) + 0 + 0)/4;
                }
            }
        // }

    }

    // correction

    // for (int i=1; i < nz-1; i++) {
        for (int j = 1; j < nx - 1; j++)
        {
            for (int k = 1; k < ny - 1; k++)
            {
                divergenceVelocityField.at(1).at(j).at(k) = (mesh.at(1).at(j+1).at(k).velocityX - mesh.at(1).at(j-1).at(k).velocityX + mesh.at(1).at(j).at(k+1).velocityY - mesh.at(1).at(j).at(k-1).velocityY)/2
                mesh.at(1).at(j).at(k).density = iterativeSolver(mesh.at(1).at(j).at(k));
                mesh.at(1).at(j).at(k).newPressure = ((mesh.at(1).at(j-1).at(k).pressure + mesh.at(1).at(j+1).at(k).newPressure + mesh.at(1).at(j).at(k-1).newPressure + mesh.at(1).at(j).at(k+1).newPressure) - divergenceVelocityField.at(1).at(j).at(k)) / 4;
                // mesh.medianSurroundingDensity = (mesh.at(1).at(j+1).at(k) + mesh.at(1).at(j-1).at(k) + mesh.at(1).at(j).at(k+1) + mesh.at(1).at(j).at(k-1) + 0 + 0)/4;
            }
        }
    // }

    // for (int i=1; i < nz-1; i++) {
        for (int j = 1; j < nx - 1; j++)
        {
            for (int k = 1; k < ny - 1; k++)
            {
                mesh.at(1).at(j).at(k).pressure = mesh.at(1).at(j).at(k).newPressure;
            }
        }
    // }

    // for (int i=1; i < nz-1; i++) {
        for (int j = 1; j < nx - 1; j++)
        {
            for (int k = 1; k < ny - 1; k++)
            {
                gradientPressureField.at(1).at(j).at(k).x = (mesh.at(1).at(j+1).at(k).pressure - mesh.at(1).at(j-1).at(k).pressure)/2;
                gradientPressureField.at(1).at(j).at(k).y = (mesh.at(1).at(j).at(k+1).pressure - mesh.at(1).at(j).at(k-1).pressure)/2;
                // gradientPressureField.at(1).at(j).at(k).z = (mesh.at(0+1).at(j).at(k).pressure - mesh.at(0-1).at(j).at(k).pressure)/2 
            }
        }
    // }
}

Cfd::Cfd(int setnx, int setny, int setnz)
{
    nx = setnx;
    ny = setny;
    nz = setnz;

    createMesh();
    setBoundaryConditions();
    setPlaneBoundary();
}

Cfd::~Cfd()
{
}