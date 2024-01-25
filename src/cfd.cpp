#include "cfd.h"

#define WITHOUT_NUMPY
#include "Physics/matplotlibcpp.h"

namespace mat = matplotlibcpp;

void Cfd::createMesh()
{
    for (int i = 0; i < nz; i++)
    {   
        std::cout << "done " << i << std::endl;
        std::vector<std::vector<MeshCube>> helper;
        std::vector<std::vector<double>> helper2;
        std::vector<std::vector<Vector3>> helper3;
        for (int j = 0; j < nx; j++)
        {
            std::vector<MeshCube> helperHelper;
            std::vector<double> helperHelper2;
            std::vector<Vector3> helperHelper3;
            for (int k = 0; k < ny; k++)
            {
                helperHelper.push_back(MeshCube());
                helperHelper2.push_back(0);
                helperHelper3.push_back({0,0,0});
            }
            helper.push_back(helperHelper);
            helper2.push_back(helperHelper2);
            helper3.push_back(helperHelper3);
        }
        mesh.push_back(helper);
        divergenceVelocityField.push_back(helper2);
        gradientPressureField.push_back(helper3);
    }
}

void Cfd::setBoundaryConditions(double velocityXDirectionStart, double velocityYDirectionStart, double velocityZDirectionStart, double velocityXDirectionEnd, double velocityYDirectionEnd, double velocityZDirectionEnd)
{
    for (int i = 0; i < nz; i++)
    {
        for (int k = 0; k < ny; k++)
        {
            mesh.at(i).at(0).at(k).boundary = true;
            mesh.at(i).at(0).at(k).velocityX = velocityXDirectionStart;
            mesh.at(i).at(0).at(k).pressure = pow(velocityXDirectionStart,2) * (rho/2.0f);

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

void Cfd::solveDensity(int i, int j, int k) {
    double currentDensity = mesh.at(i).at(j).at(k).density;
    double density1 = mesh.at(i).at(j).at(k).density; // Pn(x,y) known
    double density2 = 0; // Pn(x+1,y)
    double density3 = mesh.at(i).at(j-1).at(k).density; // Pn(x-1,y) known
    double density4 = 0; // Pn(x,y+1)
    double density5 = mesh.at(i).at(j).at(k-1).density; // Pn(x,y-1) known

    double previousDensity2, previousDensity4;

    int times = 0;
    while (times < 100 && (density2 == previousDensity2 || density4 == previousDensity4)) {
        density2 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density3 - density4 - density5;
        density4 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density2 - density3 - density5;

        times++;
    }

    mesh.at(i).at(j+1).at(k).newDensity = density2;
    mesh.at(i).at(j).at(k+1).newDensity = density4;
}

void Cfd::solveDensityFirst(int i, int j, int k) {
    double currentDensity = mesh.at(i).at(j).at(k).density;
    double density1 = 0; // Pn(x,y) known
    double density2 = 0; // Pn(x+1,y)
    double density3 = mesh.at(i).at(j-1).at(k).density; // Pn(x-1,y) known
    double density4 = 0; // Pn(x,y+1)
    double density5 = mesh.at(i).at(j).at(k-1).density; // Pn(x,y-1) known

    double previousDensity1, previousDensity2, previousDensity4;

    int times = 0;
    while (times < 100 && (density1 == previousDensity1 || density2 == previousDensity2 || density4 == previousDensity4)) {
        density1 = (currentDensity + k * ((density2 + density3 + density4 + density5) / 4.0f)) / (1.0f + k);
        density2 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density3 - density4 - density5;
        density4 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density2 - density3 - density5;

        times++;
    }
    
    mesh.at(i).at(j).at(k).newDensity = density1;
    mesh.at(i).at(j+1).at(k).newDensity = density2;
    mesh.at(i).at(j).at(k+1).newDensity = density4;
}

void Cfd::solvePressure(int i, int j, int k) {
    double pressure1 = mesh.at(i).at(j).at(k).pressure; // Pn(x,y) known
    double pressure2 = 0; // Pn(x+1,y)
    double pressure3 = mesh.at(i).at(j-1).at(k).pressure; // Pn(x-1,y) known
    double pressure4 = 0; // Pn(x,y+1)
    double pressure5 = mesh.at(i).at(j).at(k-1).pressure; // Pn(x,y-1) known

    double previousPressure2, previousPressure4;

    int times = 0;
    while (times < 100 && (pressure2 == previousPressure2 || pressure4 == previousPressure4)) {
        pressure2 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure3 - pressure4 - pressure5;
        pressure4 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure2 - pressure3 - pressure5;
        previousPressure2 = pressure2;
        previousPressure4 = pressure4;
        times++;
    }

    mesh.at(i).at(j+1).at(k).pressure = pressure2;
    mesh.at(i).at(j).at(k+1).pressure = pressure4;
}

void Cfd::solvePressureFirst(int i, int j, int k)
{
    double pressure1 = 0; ; // Dn(x,y)
    double pressure2 = 0; // Dn(x+1,y)
    double pressure3 =  mesh.at(i).at(j-1).at(k).pressure; // Dn(x-1,y)
    double pressure4 = 0; // Dn(x,y+1)
    double pressure5 =  mesh.at(i).at(j).at(k-1).pressure; // Dn(x,y-1)

    double previousPressure1, previousPressure2, previousPressure4;

    int times = 0;
    while (times < 100 && (pressure1 == previousPressure1 || pressure2 == previousPressure2 || pressure4 == previousPressure4)) {
        std::cout << pressure1 << " " << pressure2 << " " << pressure3  << " " << pressure4 << " " << pressure5 << std::endl;
        pressure1 = (pressure2 + pressure3 + pressure4 + pressure5 - divergenceVelocityField.at(i).at(j).at(k)) / 4.0f;
        pressure2 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure3 - pressure4 - pressure5;
        pressure4 = 4.0f * pressure1 + divergenceVelocityField.at(i).at(j).at(k) - pressure2 - pressure3 - pressure5;
        times++;
    }

    mesh.at(i).at(j).at(k).pressure = pressure1;
    mesh.at(i).at(j+1).at(k).pressure = pressure2;
    mesh.at(i).at(j).at(k+1).pressure = pressure4;
}

void Cfd::densityDispersion() {
    // for (int i=1; i < nz-1; i++) {
            for (int j = 1; j < nx - 1; j++)
            {
                for (int k = 1; k < ny - 1; k++)
                {
                    if (j == 1) {
                        solveDensityFirst(1, j, k);
                    } else {
                        solveDensity(1, j, k);
                    }
                    // mesh.at(1).at(j).at(k).updatedPressure = false;
                    
                    // mesh.at(1).at(j).at(k).density = iterativeSolver(mesh.at(1).at(j).at(k));
                    // mesh.medianSurroundingDensity = (mesh.at(1).at(j+1).at(k) + mesh.at(1).at(j-1).at(k) + mesh.at(1).at(j).at(k+1) + mesh.at(1).at(j).at(k-1) + 0 + 0)/4;
                }
            }
        // }

        // for (int i=1; i < nz-1; i++) {
            for (int j = 1; j < nx - 1; j++)
            {
                for (int k = 1; k < ny - 1; k++)
                {
                    mesh.at(1).at(j).at(k).density == mesh.at(1).at(j).at(k).newDensity;
                }
            }
        // }
}

void Cfd::removeDivergence() {
    // for (int i=1; i < nz-1; i++) {
        for (int j = 1; j < nx - 1; j++)
        {
            for (int k = 1; k < ny - 1; k++)
            {
                divergenceVelocityField.at(1).at(j).at(k) = (mesh.at(1).at(j+1).at(k).velocityX - mesh.at(1).at(j-1).at(k).velocityX + mesh.at(1).at(j).at(k+1).velocityY - mesh.at(1).at(j).at(k-1).velocityY)/2;
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

void Cfd::calc()
{
    double tijd = 0;
    while (tijd < maxTime)
    {
        tijd += dT;
        
        densityDispersion();

        // // for (int i=1; i < nz-1; i++) {
        //     for (int j=1; j < nx-1; j++) {
        //         for (int k=1; k < ny; k++) {
        //             double vHere = 0.25 * (mesh.at(1).at(i-1).at(j).velocityX + mesh.at(1).at(i-1).at(j+1).velocityX + mesh.at(1).at(i).at(j).velocityX + mesh.at(1).at(i).at(j+1).velocityX);
        //             float a = (nu * (u.at(i - 1).at(j) - 2 * u.at(i).at(j) + u.at(i + 1).at(j)) * pow(dxi, 2));
        //             float b = nu * (u.at(i).at(j - 1) - 2 * u.at(i).at(j) + u.at(i).at(j + 1) * pow(dyi, 2));
        //             float c = -u.at(i).at(j) * (u.at(i + 1).at(j) - u.at(i - 1).at(j)) * 0.5 * dxi;
        //             float d = -vHere * (u.at(i).at(j + 1) - u.at(i).at(j - 1)) * 0.5 * dyi;
        //             us.at(i).at(j) = u.at(i).at(j) + dT * (a + b + c + d); // nieuwe s over tijd
        //         }
        //     }
        // // }

        // for (int j = jMin; j < jMax + 1; j++)
        //     {
        //           for (int i = iMin + 1; i < iMax + 1; i++)
        //           {
        //                 std::cout << "ja 1" << std::endl;
        //                 float vHere = 0.25 * (v.at(i - 1).at(j) + v.at(i - 1).at(j + 1) + v.at(i).at(j) + v.at(i).at(j + 1));
        //                 float a = (nu * (u.at(i - 1).at(j) - 2 * u.at(i).at(j) + u.at(i + 1).at(j)) * pow(dxi, 2));
        //                 float b = nu * (u.at(i).at(j - 1) - 2 * u.at(i).at(j) + u.at(i).at(j + 1) * pow(dyi, 2));
        //                 float c = -u.at(i).at(j) * (u.at(i + 1).at(j) - u.at(i - 1).at(j)) * 0.5 * dxi;
        //                 float d = -vHere * (u.at(i).at(j + 1) - u.at(i).at(j - 1)) * 0.5 * dyi;
        //                 us.at(i).at(j) = u.at(i).at(j) + dT * (a + b + c + d); // nieuwe s over tijd
        //           }
        //     }
        //     std::cout << "ja 1 2" << std::endl;

        //     for (int j = jMin + 1; j < jMax + 1; j++)
        //     {
        //           for (int i = iMin; i < iMax + 1; i++)
        //           {
        //                 std::cout << "ja 2" << std::endl;
        //                 float uHere = 0.25 * (u.at(i).at(j - 1) + u.at(i).at(j) + u.at(i + 1).at(j - 1) + u.at(i + 1).at(j));
        //                 float a = (nu * (v.at(i - 1).at(j) - 2 * v.at(i).at(j) + v.at(i + 1).at(j)) * pow(dxi, 2));
        //                 float b = nu * (v.at(i).at(j - 1) - 2 * v.at(i).at(j) + v.at(i).at(j + 1) * pow(dyi, 2));
        //                 float c = -uHere * (v.at(i + 1).at(j) - v.at(i - 1).at(j)) * 0.5 * dyi;
        //                 float d = -v.at(i).at(j) * (v.at(1).at(j + 1) - v.at(i).at(j - 1)) * 0.5 * dxi;
        //                 vs.at(i).at(j) = v.at(i).at(j) + dT * (a + b + c + d);
        //           }
        //     }
        removeDivergence();

        Draw();
        std::cout << tijd << " " << maxTime << std::endl;
    }

    // correction

    // // for (int i=1; i < nz-1; i++) {
    //     for (int j = 1; j < nx - 1; j++)
    //     {
    //         for (int k = 1; k < ny - 1; k++)
    //         {
    //             divergenceVelocityField.at(1).at(j).at(k) = (mesh.at(1).at(j+1).at(k).velocityX - mesh.at(1).at(j-1).at(k).velocityX + mesh.at(1).at(j).at(k+1).velocityY - mesh.at(1).at(j).at(k-1).velocityY)/2;
    //             // mesh.at(1).at(j).at(k).pressure = iterativeSolver(1, j, k);
    //             if (j == 1) {
    //                 solvePressureFirst(1, j, k);
    //             } else {
    //                 solvePressure(1, j, k);
    //             }
    //             // mesh.at(1).at(j).at(k).newPressure = ((mesh.at(1).at(j-1).at(k).pressure + mesh.at(1).at(j+1).at(k).newPressure + mesh.at(1).at(j).at(k-1).newPressure + mesh.at(1).at(j).at(k+1).newPressure) - divergenceVelocityField.at(1).at(j).at(k)) / 4;
    //             // mesh.medianSurroundingDensity = (mesh.at(1).at(j+1).at(k) + mesh.at(1).at(j-1).at(k) + mesh.at(1).at(j).at(k+1) + mesh.at(1).at(j).at(k-1) + 0 + 0)/4;
    //         }
    //     }
    // // }

    // // for (int i=1; i < nz-1; i++) {
    //     for (int j = 1; j < nx - 1; j++)
    //     {
    //         for (int k = 1; k < ny - 1; k++)
    //         {
    //             mesh.at(1).at(j).at(k).pressure = mesh.at(1).at(j).at(k).newPressure;
    //         }
    //     }
    // // }

    // // for (int i=1; i < nz-1; i++) {
    //     for (int j = 1; j < nx - 1; j++)
    //     {
    //         for (int k = 1; k < ny - 1; k++)
    //         {
    //             gradientPressureField.at(1).at(j).at(k).x = (mesh.at(1).at(j+1).at(k).pressure - mesh.at(1).at(j-1).at(k).pressure)/2;
    //             gradientPressureField.at(1).at(j).at(k).y = (mesh.at(1).at(j).at(k+1).pressure - mesh.at(1).at(j).at(k-1).pressure)/2;
    //             // gradientPressureField.at(1).at(j).at(k).z = (mesh.at(0+1).at(j).at(k).pressure - mesh.at(0-1).at(j).at(k).pressure)/2 
    //         }
    //     }
    // // }
}

void Cfd::Draw() {
    while (!WindowShouldClose()) {
        BeginDrawing();

            std::cout << "start" << std::endl;
            for (int j=0; j < 100; j++) {
                for (int k=0; k < 100; k++) {
                    double val = mesh.at(1).at(j).at(k).pressure / mesh.at(1).at(1).at(1).pressure;
                    double val2 = (mesh.at(1).at(j).at(k).pressure / mesh.at(1).at(1).at(1).pressure) * 10;
                    double val3 = mesh.at(1).at(j).at(k).pressure;
                    Color col = {val3, val, val2, 255};
                    DrawRectangle(j*10, k*10, 10, 10, col);
                    DrawLine(j*10, k*10, 10 * j, 10 * k +10, RED);
                    DrawLine(j*10, k*10, 10 * j +10,10 * k, RED);
                    // DrawText()
                    // DrawLine(j*10, k*10, 10 * j +10,10 * k +10, RED);
                    // DrawLine(j*10, k*10, 10 * j +10,10 * k +10, RED);
                }
            }
        
            // for (int j=0; j< 100; j++) {
            //     for (int k=0; k < 100; k++) {
            //         std::cout << mesh.at(1).at(j).at(k).pressure << " ";
            //     }
            //     std::cout << " end " << std::endl;
            // }
            // std::cout << std::endl;
            // std::cout << std::endl;
            // std::cout << std::endl;
            // std::cout << std::endl;


        EndDrawing();
    }
}

Cfd::Cfd(int setnx, int setny, int setnz, double deltaTime, double setMaxTime, double setRho)
{
    nx = setnx;
    ny = setny;
    nz = setnz;
    dT = deltaTime;
    maxTime = setMaxTime;
    rho = setRho;

    createMesh();
    setBoundaryConditions(100,  0,  0,  0,  0,  0);
    setPlaneBoundary();
}
 
Cfd::~Cfd()
{
}
