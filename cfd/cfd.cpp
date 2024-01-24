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
            mesh.at(i).at(0).at(k).boundary = true;
            mesh.at(i).at(0).at(k).velocityX = velocityXDirectionStart;
            mesh.at(i).at(0).at(k).pressure = velocityXDirectionStart * (rho/2.0f);

            mesh.at(i).at(nx - 1).at(k).boundary = true;
            mesh.at(i).at(nx - 1).at(k).velocityX = velocityXDirectionEnd;
            // mesh.at(i).at(nx - 1).at(k).pressure = ; // set the pressure of the boundary
        }
    }

    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            mesh.at(i).at(j).at(0).boundary = true;
            mesh.at(i).at(j).at(0).velocityY = velocityYDirectionStart;

            mesh.at(i).at(j).at(ny - 1).boundary = true;
            mesh.at(i).at(j).at(ny - 1).velocityY = velocityYDirectionEnd;
        }
    }

    for (int j = 0; j < nx; j++)
    {
        for (int k = 0; k < ny; k++)
        {
            mesh.at(0).at(j).at(k).boundary = true;
            mesh.at(0).at(j).at(k).velocityZ = velocityZDirectionStart;

            mesh.at(nz - 1).at(j).at(k).boundary = true;
            mesh.at(nz - 1).at(j).at(k).velocityZ = velocityZDirectionEnd;
        }
    }
}

void Cfd::setPlaneBoundary()
{
}

void Cfd::solvePressure(int i, int j, int k) {
    double pressure1 = mesh.at(i).at(j).at(k).pressure; // Dn(x,y) known
    double pressure2 = 0; // Dn(x+1,y)
    double pressure3 = mesh.at(i).at(j-1).at(k).pressure: // Dn(x-1,y) known
    double pressure4 = 0; // Dn(x,y+1)
    double pressure5 = mesh.at(i).at(j).at(k-1).pressure; // Dn(x,y-1) known

    for (int i=0; i < 100; i++) {
        pressure2 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure3 - pressure4 - pressure5;
        pressure4 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure2 - pressure3 - pressure5;
    }
    mesh.at(i).at(j+1).at(k).pressure = pressure2;
    mesh.at(i).at(j).at(k+1).pressure = pressure4;
}

void Cfd::solvePressureFirst(int i, int j, int k)
{
    // float newDensity = 0; // Dn(x,y)
    // float x = 0;      // Dn(x+1,y)
    // float y = densityX;      // Dn(x-1,y)
    // float z = 0;      // Dn(x,y+1)
    // float w = densityY;      // Dn(x,y-1)

    // float previousNewDensity = 0; // previousDn value of (x,y) in itterative solver
    // float previousx = 0;          // previousDn value of (x+1,y) in itterative solver
    // float previousz = 0;          // previousDn value of (x,y+1) in itterative solver

    // while (newDensity != previousNewDensity || x != previousx || z != previousz)
    // {
    //     previousNewDensity = newDensity;
    //     previousx = x;
    //     previousz = z;

    //     newDensity = (density + k * ((x + y + w + z) / 4)) / (1 + k);
    //     x = 4 * newDensity + ((4 * newDensity + density) / 4) - (y + z + w);
    //     z = 4 * newDensity + ((4 * newDensity + density) / 4) - (x + y + w);
    // }
    // std::vector<Vector3> useIndex;
    // for (int x=k-1; x <= k+1; x++) {
    //     if (mesh.at(1).at(j).at(x).updatedPressure) {
    //         useIndex.push_back
    //     }
    // }
    // for (int wx=j-1; w <= j+1; w++) {
    //     if (mesh.at(1).at(j).at(w).updatedPressure) {

    //     }
    // }

    double pressure1 = 0; ; // Dn(x,y)
    double pressure2 = 0; // Dn(x+1,y)
    double pressure3 =  mesh.at(i).at(j-1).at(k).pressure; // Dn(x-1,y)
    double pressure4 = 0; // Dn(x,y+1)
    double pressure5 =  mesh.at(i).at(j).at(k-1).pressure; // Dn(x,y-1)

    for (int i=0; i < 100; i++) {
        std::cout << pressure1 << " " << pressure2 << " " << pressure3  << " " << pressure4 << " " << pressure5 << std::endl;
        pressure1 = (pressure2 + pressure3 + pressure4 + pressure5 - divergenceVelocityField.at(i).at(j).at(k)) / 4.0f;
        pressure2 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure3 - pressure4 - pressure5;
        pressure4 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure2 - pressure3 - pressure5;
    }

    mesh.at(i).at(j).at(k).pressure = pressure1;
    mesh.at(i).at(j+1).at(k).pressure = pressure2;
    mesh.at(i).at(j).at(k+1).pressure = pressure4;
    
    
    // for (int x=k-1; x <= k+1; x++) {
    //     mesh.at(1).at(j).at(x).updatedPressure = true;
    // }
    // for (int wx=j-1; w <= j+1; w++) {
    //     mesh.at(1).at(j).at(w).updatedPressure = true;
    // }
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
                    // mesh.at(1).at(j).at(k).updatedPressure = false;
                    
                    // mesh.at(1).at(j).at(k).density = iterativeSolver(mesh.at(1).at(j).at(k));
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
                // mesh.at(1).at(j).at(k).pressure = iterativeSolver(1, j, k);
                if (j == 1) {
                    solvePressureFirst(1, j, k);
                } else {
                    solvePressure(1, j, k);
                }
                // mesh.at(1).at(j).at(k).newPressure = ((mesh.at(1).at(j-1).at(k).pressure + mesh.at(1).at(j+1).at(k).newPressure + mesh.at(1).at(j).at(k-1).newPressure + mesh.at(1).at(j).at(k+1).newPressure) - divergenceVelocityField.at(1).at(j).at(k)) / 4;
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

// void testDraw() {
//     InitWindow(0, 0, "airplane simulation");
//     ToggleFullscreen();
//     const int screenWidth = GetScreenWidth();
//     const int screenHeight = GetScreenHeight();

//     SetTargetFPS(60);

//     while (!WindowShouldClose()) {
//         BeginDrawing();

//             for (int j=0; j < 30; j++) {
//                 for (int k=0; k < 50; k++) {
//                     int val = mesh.at(1).at(j).at(k);
//                     Color col = {255, val, 255, 255};
//                     DrawRectangle(j*30, k*30, 30, 30, col);
//                     DrawLine(j*30, k*30, 30 * j, 30 * k +30, RED);
//                     DrawLine(j*30, k*30, 30 * j +30,30 * k, RED);
//                     // DrawText()
//                     // DrawLine(j*10, k*10, 10 * j +10,10 * k +10, RED);
//                     // DrawLine(j*10, k*10, 10 * j +10,10 * k +10, RED);
//                 }
//             }

//         EndDrawing();
//     }
// }