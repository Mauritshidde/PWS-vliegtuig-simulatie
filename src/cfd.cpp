#include "cfd.h"

#define WITHOUT_NUMPY
#include "Physics/matplotlibcpp.h"

namespace mat = matplotlibcpp;

void Cfd::Start() {
    
}

void Cfd::createMesh()
{
    for (int i = 0; i < nz; i++)
    {   
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

bool Cfd::getCollisionPlaneRay(Vector3 direction, Vector3 oppositeDirection, Ray ray, Ray ray2) {
    ray.direction = direction;
    ray2.direction = oppositeDirection;

    RayCollision meshHitInfo = GetRayCollisionMesh(ray, *airplane.meshes, airplane.transform);
    RayCollision meshHitInfo2 = GetRayCollisionMesh(ray2, *airplane.meshes, airplane.transform);

    if (meshHitInfo.hit && meshHitInfo2.hit) {
        return true;
    } else {
        return false;
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
                if (CheckCollisionBoxSphere(boundingBoxPlane, position, dx)) {
                        // ray.position = position;
                        // ray.direction = rayDirection;
                        // int collisions = plane.detectCollision(ray);

                        // if (collisions % 2 != 0 && collisions > 0) {
                        //     mesh.at(i).at(j).at(k).boundary = true;
                        // }
                        Ray ray;
                        Ray ray2;
                        
                        ray.position = position;
                        ray2.position = position;
                        
                        // check if the cube is inside the plane
                        if (getCollisionPlaneRay({1,0,0}, {1,0,0}, ray, ray2)) { 
                            if (getCollisionPlaneRay({0,1,0}, {0,1,0}, ray, ray2)) {
                                if (getCollisionPlaneRay({0,0,1}, {0,0,1}, ray, ray2)) {
                                    mesh.at(i).at(j).at(k).boundary = true;
                                }
                            }
                        }
                }
            }
        }
        std::cout << "ja" << std::endl;
    }
    settingPlaneBOundarys = false;

}

void Cfd::setPlaneBoundary() //222
{
    std::vector<std::thread> threads;
    cores = 10;
    int newNz = nz - 2;
    
    for (int i=0; i < cores; i++) {
        int part = 1 + (i * newNz)/cores;
        int part2 = 1 + ((i+1) * newNz)/cores;
        // std::cout << part << " " << part2 << " " << cores << " " << nz << std::endl;
        // std::thread t1(&Cfd::setPlaneBoundaryHelper, this, part, part2);
        threads.emplace_back(&Cfd::setPlaneBoundaryHelper, this, part, part2);
    }

    for (int i=0; i < cores; i++) {
        threads.at(i).join();
    }
}

// void Cfd::solveDensity(int i, int j, int k) {
//     double currentDensity = mesh.at(i).at(j).at(k).density;
//     double density1 = mesh.at(i).at(j).at(k).density; // Pn(x,y) known
//     double density2 = 0; // Pn(x+1,y)
//     double density3 = mesh.at(i).at(j-1).at(k).density; // Pn(x-1,y) known
//     double density4 = 0; // Pn(x,y+1)
//     double density5 = mesh.at(i).at(j).at(k-1).density; // Pn(x,y-1) known

//     double previousDensity2, previousDensity4;

//     int times = 0;
//     while (times < 100 && (density2 == previousDensity2 || density4 == previousDensity4)) {
//         density2 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density3 - density4 - density5;
//         density4 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density2 - density3 - density5;

//         times++;
//     }

//     mesh.at(i).at(j+1).at(k).newDensity = density2;
//     mesh.at(i).at(j).at(k+1).newDensity = density4;
// }

// void Cfd::solveDensityFirst(int i, int j, int k) {
//     double currentDensity = mesh.at(i).at(j).at(k).density;
//     double density1 = 0; // Pn(x,y) known
//     double density2 = 0; // Pn(x+1,y)
//     double density3 = mesh.at(i).at(j-1).at(k).density; // Pn(x-1,y) known
//     double density4 = 0; // Pn(x,y+1)
//     double density5 = mesh.at(i).at(j).at(k-1).density; // Pn(x,y-1) known

//     double previousDensity1, previousDensity2, previousDensity4;

//     int times = 0;
//     while (times < 100 && (density1 == previousDensity1 || density2 == previousDensity2 || density4 == previousDensity4)) {
//         density1 = (currentDensity + k * ((density2 + density3 + density4 + density5) / 4.0f)) / (1.0f + k);
//         density2 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density3 - density4 - density5;
//         density4 = (4.0f * density1 * (1.0f + k) - 4.0f * currentDensity) / k - density2 - density3 - density5;

//         times++;
//     }
    
//     mesh.at(i).at(j).at(k).newDensity = density1;
//     mesh.at(i).at(j+1).at(k).newDensity = density2;
//     mesh.at(i).at(j).at(k+1).newDensity = density4;
// }

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

// void Cfd::densityDispersion() {
//     // for (int i=1; i < nz-1; i++) {
//             for (int j = 1; j < nx - 1; j++)
//             {
//                 for (int k = 1; k < ny - 1; k++)
//                 {
//                     if (j == 1) {
//                         solveDensityFirst(1, j, k);
//                     } else {
//                         solveDensity(1, j, k);
//                     }
//                     // mesh.at(1).at(j).at(k).updatedPressure = false;
                    
//                     // mesh.at(1).at(j).at(k).density = iterativeSolver(mesh.at(1).at(j).at(k));
//                     // mesh.medianSurroundingDensity = (mesh.at(1).at(j+1).at(k) + mesh.at(1).at(j-1).at(k) + mesh.at(1).at(j).at(k+1) + mesh.at(1).at(j).at(k-1) + 0 + 0)/4;
//                 }
//             }
//         // }

//         // for (int i=1; i < nz-1; i++) {
//             for (int j = 1; j < nx - 1; j++)
//             {
//                 for (int k = 1; k < ny - 1; k++)
//                 {
//                     mesh.at(1).at(j).at(k).density == mesh.at(1).at(j).at(k).newDensity;
//                 }
//             }
//         // }
// }

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

    // for (int i=1; i < nz-1; i++) {
        for (int j = 1; j < nx - 1; j++)
        {
            for (int k = 1; k < ny - 1; k++)
            {
                mesh.at(1).at(j).at(k).velocityX -= gradientPressureField.at(1).at(j).at(k).x;
                mesh.at(1).at(j).at(k).velocityY -= gradientPressureField.at(1).at(j).at(k).y;
                // mesh.at(1).at(j).at(k).velocityZ -= gradientPressureField.at(1).at(j).at(k).z;
            }
        }
    // }
}

void Cfd::resetMesh() {
    for (int i = 1; i < nz-1; i++)
    {   
        for (int j = 1; j < nx-1; j++)
        {
            for (int k = 1; k < ny-1; k++)
            {
                mesh.at(i).at(j).at(k) = MeshCube();
            }
        }
    }
}

void Cfd::velocityMovement(float dT) {
    // for (int i=1; i < nz-1; i++) {
            for (int j=1; j < nx-1; j++) {
                for (int k=1; k < ny-1; k++) {
                    double vHere = 0.25 * (mesh.at(1).at(j-1).at(k).velocityY + mesh.at(1).at(j-1).at(k+1).velocityY + mesh.at(1).at(j).at(k).velocityY + mesh.at(1).at(j).at(k+1).velocityY);
                    float a = (nu * (mesh.at(1).at(j-1).at(k).velocityX - 2 * mesh.at(1).at(j).at(k).velocityX + mesh.at(1).at(j + 1).at(k).velocityX) * pow(dxi, 2));
                    float b = nu * (mesh.at(1).at(j).at(k - 1).velocityX - 2 * mesh.at(1).at(j).at(k).velocityX + mesh.at(1).at(j).at(k + 1).velocityX * pow(dyi, 2));
                    float c = -mesh.at(1).at(j).at(k).velocityX * (mesh.at(1).at(j +1).at(k).velocityX - mesh.at(1).at(j-1).at(k).velocityX) * 0.5 * dxi;
                    float d = -vHere * (mesh.at(1).at(j).at(k+1).velocityX - mesh.at(1).at(j).at(k-1).velocityX) * 0.5 * dyi;
                    mesh.at(1).at(j).at(k).newVelocityX = mesh.at(1).at(j).at(k).velocityX + dT * (a + b + c + d); // nieuwe s over tijd
                }
            }
        // }
    // for (int i=1; i < nz-1; i++) {
        for (int j=1; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                float vHere = 0.25 * (mesh.at(1).at(j-1).at(k).velocityX + mesh.at(1).at(j - 1).at(k+1).velocityX + mesh.at(1).at(j).at(k).velocityX + mesh.at(1).at(j).at(k+1).velocityX);
                float a = (nu * (mesh.at(1).at(j-1).at(k).velocityY - 2 * mesh.at(1).at(j).at(k).velocityY + mesh.at(1).at(j+1).at(k).velocityY) * pow(dxi, 2));
                float b = nu * (mesh.at(1).at(j).at(k-1).velocityY - 2 * mesh.at(1).at(j).at(k).velocityY + mesh.at(1).at(j).at(k+1).velocityY * pow(dyi, 2));
                float c = -mesh.at(1).at(j).at(k).velocityY * (mesh.at(1).at(j + 1).at(k).velocityY - mesh.at(1).at(j -1).at(k).velocityY) * 0.5 * dxi;
                float d = -vHere * (mesh.at(1).at(j).at(k + 1).velocityY - mesh.at(1).at(j).at(k-1).velocityY) * 0.5 * dyi;
                mesh.at(1).at(j).at(k).newVelocityY = mesh.at(1).at(j).at(k).velocityY + dT * (a + b + c + d); // nieuwe s over tijd
            }
        }
    // }

    // for (int i=1; i < nz-1; i++) {
        for (int j=1; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                mesh.at(1).at(j).at(k).velocityX = mesh.at(1).at(j).at(k).newVelocityX;
                mesh.at(1).at(j).at(k).velocityY = mesh.at(1).at(j).at(k).newVelocityY;
                // std::cout << "newe x" << mesh.at(1).at(j).at(k).newVelocityX << " y " << mesh.at(1).at(j).at(k).newVelocityY << std::endl;
            }
        }
    // }
}

Vector2 Cfd::calc(double anglePitch, double angleYaw)
{
    float cl, cd;
    double tijd = 0;
    while (tijd < maxTime)
    {
        tijd += dT;
        
        // TODO the movement of the velocity and pressure NOTE density is constant

        velocityMovement(dT);
        removeDivergence();


        // for (int j = jMin + 1; j < jMax + 1; j++)
        // {
        //     for (int i = iMin; i < iMax + 1; i++)
        //     {
        //         std::cout << "ja 2" << std::endl;
        //         float uHere = 0.25 * (u.at(i).at(j - 1) + u.at(i).at(j) + u.at(i + 1).at(j - 1) + u.at(i + 1).at(j));
        //         float a = (nu * (v.at(i - 1).at(j) - 2 * v.at(i).at(j) + v.at(i + 1).at(j)) * pow(dxi, 2));
        //         float b = nu * (v.at(i).at(j - 1) - 2 * v.at(i).at(j) + v.at(i).at(j + 1) * pow(dyi, 2));
        //         float c = -uHere * (v.at(i + 1).at(j) - v.at(i - 1).at(j)) * 0.5 * dyi;
        //         float d = -v.at(i).at(j) * (v.at(1).at(j + 1) - v.at(i).at(j - 1)) * 0.5 * dxi;
        //         vs.at(i).at(j) = v.at(i).at(j) + dT * (a + b + c + d);
        //     }
        // }

        if (drawing) {
            Draw();
        }
        std::cout << tijd << " " << maxTime << std::endl;
    }

    // TODO correction fase
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

    // TODO use pressure field to calculate lift and drag
    return {cl, cd};
}

void Cfd::moveCamera(float deltaTime) {
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

void Cfd::drawVelocityVectors() {
    for (int i=1; i < 2; i++) {
        for (int j=1; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                Vector3 point;
                point.x = startingPoint.x + j * dx - 0.5 * dx;
                point.y = startingPoint.y + k * dy - 0.5 * dy;
                point.z = startingPoint.z + i * dz - 0.5 * dz;
                if (mesh.at(i).at(j).at(k).boundary) {
                    // DrawCubeWires(point, dx, dy, dz, BLACK);
                    DrawCube(point, dx, dy, dz, BLACK);
                } else {
                    
                    float velocityX = mesh.at(i).at(j).at(k).velocityX;
                    float velocityY = mesh.at(i).at(j).at(k).velocityY;
                    float velocityZ = mesh.at(i).at(j).at(k).velocityZ;
                    float velocity = sqrt(pow(velocityX,2) + pow(velocityY,2) + pow(velocityZ,2));
                    
                    // TODO make vector whichs color depends on velocity
                    double val = (velocity / 200.0f) *30;
                    double val2 = (velocity / 500.0f);
                    // double val2 = (mesh.at(1).at(j).at(k).pressure / mesh.at(1).at(1).at(1).pressure) * 10;
                    // double val3 = mesh.at(1).at(j).at(k).pressure;
                    Color velocityColor = {255, val2, val, 255};
                    
                    Vector3 velocityDirection = {velocityX,velocityY,velocityZ};
                    velocityDirection = Vector3Normalize2(velocityDirection);
                    velocityDirection.x = (velocityDirection.x * 0.5 * dx + point.x) *10;
                    velocityDirection.y = (velocityDirection.y * 0.5 * dy + point.y) *10;
                    velocityDirection.z = (velocityDirection.z * 0.5 * dz + point.z) *10;
                    // std::cout << point.x << "  x " << velocityDirection.x << std::endl;
                    // std::cout << point.y << " y " << velocityDirection.y << std::endl;
                    // std::cout << point.z << " z " << velocityDirection.z << std::endl;
                    DrawLine3D(point, velocityDirection, velocityColor); //111
                    // DrawLine3D(point, {point.x, point.y, point.z+dz}, BLUE);
                    // DrawCubeWires(point, dx, dy, dz, RED);
                    std::cout << velocity << " ";
                }
            }
            std::cout  << std::endl;
            // std::cout << "velocity start " << mesh.at(1).at(0).at(0).velocityX << std::endl; 
        }
    }
            std::cout  << std::endl;
            std::cout  << std::endl;
            std::cout  << std::endl;
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

    moveCamera(GetFrameTime());
    BeginDrawing();
        ClearBackground(WHITE);
        BeginMode3D(camera);
            // DrawModel(airplane, {0,0,0}, 1.0f, RED);
            // DrawBoundingBox(boundingBoxPlane, ORANGE);
            DrawPoint3D(boundingBoxPlaneMax, PINK);
            DrawPoint3D(boundingBoxPlaneMin, BLACK);
            // std::cout << boundingBoxPlaneMax.x << " "  << boundingBoxPlaneMax.y << " " << boundingBoxPlaneMax.z << std::endl;
            drawVelocityVectors();
        EndMode3D();
    EndDrawing();
}

void Cfd::run(int steps) { //333
    double stepsize = 360.0f/steps;
    std::vector<std::vector<Vector2>> cfdResults;
    for (double i=0; i < 360; i+=stepsize) { // pitch
        std::vector<Vector2> cfdResultsHelper;
        for (double j=0; j < 360; j+=stepsize) { // yaw
            airplane.transform = MatrixRotateXYZ2((Vector3){DEG2RAD * i, DEG2RAD * j, DEG2RAD * 0});
            resetMesh();
            // setPlaneBoundary();
            Vector2 consts = calc(i, j);
            cfdResultsHelper.push_back({consts.x, consts.y});
        }
    }
}

Cfd::Cfd(int setnx, int setny, int setnz, double deltaTime, double setMaxTime, double setRho, bool drawingEnabled)
{   
    float Re = 100; // Reynolds number
    nu = 1 / Re;

    // set multithreading variables
    cores = 6;
    settingPlaneBOundarys = false;

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
    airplane.transform = MatrixTranslate2(0, 0.0f, 0);

    // set simulation variables
    plane.loadObjectModel();

    boundingBoxPlane = GetModelBoundingBox(airplane);
    Vector3 boundingBoxPlaneMin = boundingBoxPlane.min;
    Vector3 boundingBoxPlaneMax = boundingBoxPlane.max;

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

    drawing = drawingEnabled;
    if (!drawingEnabled) {
        CloseWindow();
    }

    // functions for generating the mesh
    createMesh();
    setBoundaryConditions(100,  0,  0,  0,  0,  0);
}
 
Cfd::~Cfd()
{
}
