#include "cfd.h"

#define WITHOUT_NUMPY
#include "Physics/matplotlibcpp.h"

namespace mat = matplotlibcpp;

void Cfd::createMesh()
{
    for (int i = 0; i < nz; i++)
    {   
        // std::cout << "done " << i << std::endl;
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

void Cfd::setPlaneBoundaryHelper(int startIndex, int endIndex) {
    for (int i=startIndex; i < endIndex; i++) {
        for (int j=1; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                Vector3 position;
                position.x = dx * j + startingPoint.x;
                position.y = dy * k + startingPoint.y;
                position.z = dz * i + startingPoint.z;
                Vector3 test = {0, 1, 0};
                Ray ray;
                ray.position = position;
                ray.direction = test;
                int collisions = plane.detectCollision(ray);

                // if (collisions > 0) {
                //     std::cout << collisions << std::endl;
                // }

                if (collisions % 2 != 0 && collisions > 0) {
                    mesh.at(i).at(j).at(k).boundary = true;
                }
            }
        }
        std::cout << "ja" << std::endl;
    }
}

void Cfd::setPlaneBoundary()
{
    int part1 = (int) nz/5.0f;
    int part2 = (int) 2.0f * nz/5.0f;
    int part3 = (int) 3.0f * nz/3.0f;
    int part4 = (int) 4.0f * nz/5.0f;
    std::cout << part1 << " " << part2 << " " << nz << std::endl;
    std::thread t1(&Cfd::setPlaneBoundaryHelper, this, 1, part1);
    std::thread t2(&Cfd::setPlaneBoundaryHelper, this, part1, part2);
    std::thread t3(&Cfd::setPlaneBoundaryHelper, this, part2, part3);
    std::thread t4(&Cfd::setPlaneBoundaryHelper, this, part3, part4);
    std::thread t5(&Cfd::setPlaneBoundaryHelper, this, part4, nz-1);

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
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

void Cfd::calc(double anglePitch, double angleYaw)
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

void Cfd::moveCamera() {
    float deltaTime = 0.01;
    Vector2 currentMousePos = GetMousePosition();
    
    if (IsMouseButtonDown(0))
    {
        angleYAxis += 100 * ((currentMousePos.x - previousMousePosition.x)) * deltaTime;
        angleXZAxis += 100 * ((currentMousePos.y - previousMousePosition.y)) * deltaTime;
        if (angleYAxis > 360) {
            angleYAxis -= 360;
        } else if (angleYAxis < 0) {
            angleYAxis += 360;
        }

        if (angleXZAxis > 360) {
            angleYAxis -= 360;
        } else if (angleYAxis < 0) {
            angleXZAxis += 360;
        }
    }

    if (IsKeyDown(KEY_RIGHT))
    {
        angleYAxis += 100 * deltaTime;
        if (angleYAxis > 360) {
            angleYAxis -= 360;
        }
    }
    if (IsKeyDown(KEY_LEFT))
    {
        angleYAxis -= 100 * deltaTime;
        if (angleYAxis < 0) {
            angleYAxis += 360;
        }
    }
    if (IsKeyDown(KEY_UP))
    {
        angleXZAxis += 100 * deltaTime;
        if (angleXZAxis > 360) {
            angleXZAxis -= 360;
        }
    }
    if (IsKeyDown(KEY_DOWN))
    {
        angleXZAxis -= 100 * deltaTime;
        if (angleXZAxis < 0) {
            angleXZAxis += 360;
        }
    }

    if (GetMouseWheelMove() > 0)
    {
        cameraCircleRadius += 100 * deltaTime;
        cameraPos = {0.0f, 0.0f, cameraCircleRadius};
    }
    else if (GetMouseWheelMove() < 0)
    {
        cameraCircleRadius -= 100 * deltaTime;
        cameraPos = {0.0f, 0.0f, cameraCircleRadius};
    }

    camera.position = Vector3Transform2(cameraPos, MatrixRotateXYZ2((Vector3){DEG2RAD * angleXZAxis, DEG2RAD * angleYAxis, 0}));
    previousMousePosition = currentMousePos;
}

void Cfd::Draw() {
    Vector3 position;
    position.x = 0;
    position.y = 0;
    position.z = 0;
    Vector3 test = {0, 100, 0};
    Ray ray;
    ray.position = position;
    ray.direction = test;

    int collisions = plane.detectCollision(ray);
    // std::cout << collisions << std::endl;
    // float deltaTime = GetFrameTime();
    moveCamera();
    BeginDrawing();
        ClearBackground(WHITE);
        BeginMode3D(camera);
        plane.drawModel();
        // DrawRay(ray, PINK);
    // DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f}, (Vector3){1.0f, 0.0f, 0.0f}, 0, (Vector3){0.5f, 0.5f, 0.5f}, WHITE); // 2de vector geeft aan met welke factor hij met currentangle draait

            DrawCubeWires({0,0,0}, 20, 40, 40, RED);

            // std::cout << "start" << std::endl;
            for (int i=1; i < nz; i++) {
            // std::cout << "start" << std::endl;
                
                for (int j=1; j < nx; j++) {
            // std::cout << j << std::endl;

                    for (int k=1; k < ny; k++) {
                        // double val = mesh.at(1).at(j).at(k).pressure / mesh.at(1).at(1).at(1).pressure;
                        // double val2 = (mesh.at(1).at(j).at(k).pressure / mesh.at(1).at(1).at(1).pressure) * 10;
                        // double val3 = mesh.at(1).at(j).at(k).pressure;
                        // Color col = {val3, val, val2, 255};
                        Vector3 point;
                        point.x = startingPoint.x + j * dx - 0.5 * dx;
                        point.y = startingPoint.y + k * dy - 0.5 * dy;
                        point.z = startingPoint.z + i * dz - 0.5 * dz;
                        if (mesh.at(i).at(j).at(k).boundary) {
                            // DrawCubeWires(point, dx, dy, dz, BLACK);
                            DrawCube(point, dx, dy, dz, BLACK)
                        } else {
                            // DrawCubeWires(point, dx, dy, dz, RED);
                        }
                    }
                }
            }
        EndMode3D();
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

void Cfd::run(int steps) {
    while (true) {
        Draw();
    }
    // for (int i=0; i < 360; i++) { // pitch
    //     for (int j=0; j < 360; j++) { // yaw
    //         calc(i, j);
    //     }
    // }
}

Cfd::Cfd(int setnx, int setny, int setnz, double deltaTime, double setMaxTime, double setRho)
{
    // set multithreading variables
    cores = 6;

    // set camera variables
    cameraCircleRadius = 150;
    cameraPos = {0.0f, 0.0f, cameraCircleRadius};
    cameraXYPos = {cameraPos.x, cameraPos.y};
    camera = {0};
    angleYAxis = 0;
    angleXZAxis = 0;

    camera.position = cameraPos;                  // Camera position perspective
    camera.target = (Vector3){0.0f, 0.0f, 0.0f}; // Camera looking at point  20 ?????????????? hier naar nog kijken TODO
    camera.up = (Vector3){0.0f, 30.0f, 0.0f};     // Camera up vector (rotation towards target)
    camera.fovy = 30.0f;                          // Camera field-of-view Y   effect van dit veranderen bestuderen ?????????????? TODO
    camera.projection = CAMERA_PERSPECTIVE;  /// wat doet dit TODO

    // set plane model variables
    airplane = LoadModel("models/object/airplane.obj");
    airplaneTexture = LoadTexture("models/texture/planeTextureBeter.png");
    airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;
    airplane.transform = MatrixTranslate2(0, -10.0f, 0);

    // set simulation variables
    plane.loadObjectModel();
    nx = setnx;
    ny = setny;
    nz = setnz;
    dT = deltaTime;
    maxTime = setMaxTime;
    rho = setRho;
    dx = 1.5;
    dy = 1;
    dz = 5;
    startingPoint.x = -(nx*dx)/2;
    startingPoint.y = -(ny*dy)/2 + 10;
    startingPoint.z = -(nz*dz)/2;

    // need to be replaced
    createMesh();
    setBoundaryConditions(100,  0,  0,  0,  0,  0);
    setPlaneBoundary();
}
 
Cfd::~Cfd()
{
}
