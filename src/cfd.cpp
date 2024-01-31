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
        divergenceVelocityScalarField.push_back(helper2);
        gradientPressureField.push_back(helper3);
        divergenceVelocityField.push_back(helper3);
        divergenceFreeField.push_back(helper3);
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
            // mesh.at(i).at(0).at(k).pressure = pow(velocityXDirectionStart,2) * (rho/2.0f);

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

void Cfd::solvePressure(int i, int j, int k)
{
    if (!mesh.at(i).at(j).at(k).boundary) {
        std::vector<Vector3> toUse;
        if (!mesh.at(i).at(j+1).at(k).boundary) {
            toUse.push_back({i,j+1,k});
        }
        if (!mesh.at(i).at(j-1).at(k).boundary) {
            toUse.push_back({i,j-1,k});
        }
        if (!mesh.at(i).at(j).at(k+1).boundary) {
            toUse.push_back({i,j,k+1});
        }
        if (!mesh.at(i).at(j).at(k-1).boundary) {
            toUse.push_back({i,j,k-1});
        }
        // if (!mesh.at(i+1).at(j).at(k).boundary) {
        //     toCalc.push_back({i+1,j,k});
        // }
        // if (!mesh.at(i-1).at(j).at(k).boundary) {
        //     toCalc.push_back({i-1,j,k});
        // }
        Vector3 indexes;
        std::vector<Vector3> toCalc;
        std::vector<double> values;
        std::vector<double> previousValues;
        for (int l=0; l < toUse.size(); l++) {
            indexes = toUse.at(l);
            if (!mesh.at(indexes.x).at(indexes.y).at(indexes.z).pressureChanged) {
                toCalc.push_back(indexes);
                values.push_back(0);
                previousValues.push_back(mesh.at(indexes.x).at(indexes.y).at(indexes.z).pressure);
            }
        }

        double pressure1 = 0; // Dn(x,y)

        int times = 0;
        bool end = false;
        while (times < 100 && !end) {
            if (!mesh.at(indexes.x).at(indexes.y).at(indexes.z).pressureChanged) {
                double val = 0;
                for (int l=0; l < toUse.size(); l++) {
                    indexes = toUse.at(l);
                    val += mesh.at(indexes.x).at(indexes.y).at(indexes.z).pressure;
                }
                pressure1 = (val - divergenceVelocityScalarField.at(i).at(j).at(k)) / 4.0f;
            }   

            for (int l=0; l < toCalc.size(); l++) {
                double val = 0;
                for (int m=0; m < toUse.size(); m++) {
                    indexes = toUse.at(m);
                    val += mesh.at(indexes.x).at(indexes.y).at(indexes.z).pressure;
                }
                values.at(l) = 4.0f * pressure1 + divergenceVelocityScalarField.at(i).at(j).at(k) - val;
            }

            for (int l=0; l < toCalc.size(); l++) {
                end = true;
                if (values.at(l) != previousValues.at(l)) {
                    end = false;
                }
            } 

            times++;
        }

        for (int l=0; l < toCalc.size(); l++) {
            indexes = toUse.at(l);
            mesh.at(indexes.x).at(indexes.y).at(indexes.z).pressureChanged = true;
        }
    }
}

void Cfd::removeDivergence() {
    // for (int i=1; i < nz-1; i++) {
        for (int j = 2; j < nx - 2; j++)
        {
            for (int k = 2; k < ny - 2; k++)
            { // TODO check if boundary for the complete code of the removeDivergence function
                if (!mesh.at(1).at(j).at(k).boundary) {
                    divergenceVelocityScalarField.at(1).at(j).at(k) = (mesh.at(1).at(j+1).at(k).velocityX - mesh.at(1).at(j-1).at(k).velocityX + mesh.at(1).at(j).at(k+1).velocityY - mesh.at(1).at(j).at(k-1).velocityY)/2;
                    divergenceVelocityField.at(1).at(j).at(k).x = divergenceVelocityScalarField.at(1).at(j).at(k);
                    divergenceVelocityField.at(1).at(j).at(k).y = divergenceVelocityScalarField.at(1).at(j).at(k);
                    divergenceVelocityField.at(1).at(j).at(k).z = divergenceVelocityScalarField.at(1).at(j).at(k);
                } else {
                    divergenceVelocityScalarField.at(1).at(j).at(k) = 0;
                    divergenceVelocityField.at(1).at(j).at(k) = {0,0,0};
                }
                solvePressure(1,j,k);
            }
        }
    // }

    // for (int i=1; i < nz-1; i++) {
        for (int j = 2; j < nx - 2; j++)
        {
            for (int k = 2; k < ny - 2; k++)
            {
                mesh.at(1).at(j).at(k).pressureChanged = false;
                if (!mesh.at(1).at(j).at(k).boundary) {
                    gradientPressureField.at(1).at(j).at(k).x = (mesh.at(1).at(j+1).at(k).pressure - mesh.at(1).at(j-1).at(k).pressure)/2;
                    gradientPressureField.at(1).at(j).at(k).y = (mesh.at(1).at(j).at(k+1).pressure - mesh.at(1).at(j).at(k-1).pressure)/2;
                    gradientPressureField.at(1).at(j).at(k).z = (mesh.at(0+1).at(j).at(k).pressure - mesh.at(2-1).at(j).at(k).pressure)/2; 
                } else {
                    gradientPressureField.at(1).at(j).at(k) = {0,0,0};
                }
            }
        }
    // }



    // for (int i=1; i < nz-1; i++) {
        for (int j = 2; j < nx - 2; j++)
        {
            for (int k = 2; k < ny - 2; k++)
            {
                if (!mesh.at(1).at(j).at(k).boundary) {
                    divergenceFreeField.at(1).at(j).at(k).x = divergenceVelocityField.at(1).at(j).at(k).x - gradientPressureField.at(1).at(j).at(k).x;
                    divergenceFreeField.at(1).at(j).at(k).y = divergenceVelocityField.at(1).at(j).at(k).y - gradientPressureField.at(1).at(j).at(k).y;
                    divergenceFreeField.at(1).at(j).at(k).z = divergenceVelocityField.at(1).at(j).at(k).z - gradientPressureField.at(1).at(j).at(k).z;

                    mesh.at(1).at(j).at(k).velocityX += (mesh.at(1).at(j).at(k).velocityX - divergenceFreeField.at(1).at(j).at(k).x)/dx *(dT/10);
                    mesh.at(1).at(j).at(k).velocityY += (mesh.at(1).at(j).at(k).velocityY - divergenceFreeField.at(1).at(j).at(k).y)/dy *(dT/10);
                    mesh.at(1).at(j).at(k).velocityZ += (mesh.at(1).at(j).at(k).velocityZ - divergenceFreeField.at(1).at(j).at(k).z)/dz *(dT/10);
                } else {
                    divergenceFreeField.at(1).at(j).at(k) = {0,0,0};
                }
                // TODO solve velocity by calculating the divergence velocity back
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
                divergenceVelocityField.at(i).at(j).at(k) = {0,0,0};
            }
        }
    }
}

void Cfd::velocityMovement(float dT) {
    std::vector<std::vector<std::vector<Vector3>>> tempVelocity;
    for (int i=0; i < nz; i++) {
        std::vector<std::vector<Vector3>> helper;
        for (int j=0; j < nx; j++) {
            std::vector<Vector3> helperHelper;
            for (int k=0; k < ny; k++) {
                helperHelper.push_back({0,0,0});
            }
            helper.push_back(helperHelper);
        }
        tempVelocity.push_back(helper);
    }

    // for (int i=1; i < nz-1; i++) {
        for (int j=1; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                if (!mesh.at(1).at(j).at(k).boundary) {
                    // double duDt = 1(mesh.at())
                    double duDt = -(mesh.at(1).at(j).at(k).velocityX * (mesh.at(1).at(j).at(k).velocityX - mesh.at(1).at(j-1).at(k).velocityX) / dx +
                            mesh.at(1).at(j).at(k).velocityY * (mesh.at(1).at(j).at(k).velocityX - mesh.at(1).at(j).at(k-1).velocityX) / dy + 0) /dx;
                            // mesh.at(1).at(j).at(k).velocityZ * (mesh.at(1).at(j).at(k).velocityX - mesh.at(2-1).at(j).at(k).velocityX) / dz) / dx;

                    double dvDt = -(mesh.at(1).at(j).at(k).velocityX * (mesh.at(1).at(j).at(k).velocityY - mesh.at(1).at(j-1).at(k).velocityY) / dx +
                            mesh.at(1).at(j).at(k).velocityY * (mesh.at(1).at(j).at(k).velocityY - mesh.at(1).at(j).at(k-1).velocityY) / dy + 0) /dy;
                            // mesh.at(1).at(j).at(k).velocityZ * (mesh.at(1).at(j).at(k).velocityY - mesh.at(2-1).at(j).at(k).velocityY) / dz) / dy;

                    // double dwDt = -(mesh.at(1).at(j).at(k).velocityX * (mesh.at(1).at(j).at(k).velocityZ - mesh.at(1).at(j-1).at(k).velocityZ) / dx +
                    //         mesh.at(1).at(j).at(k).velocityY * (mesh.at(1).at(j).at(k).velocityZ - mesh.at(1).at(j).at(k-1).velocityZ) / dy +
                    //         mesh.at(1).at(j).at(k).velocityZ * (mesh.at(1).at(j).at(k).velocityZ - mesh.at(2-1).at(j).at(k).velocityZ) / dz) / dz;


                    tempVelocity.at(1).at(j).at(k).x = mesh.at(1).at(j).at(k).velocityX + duDt * dT;
                    tempVelocity.at(1).at(j).at(k).y = mesh.at(1).at(j).at(k).velocityY + dvDt * dT;
                    // tempVelocity.at(1).at(j).at(k).z = mesh.at(1).at(j).at(k).velocityZ + dwDt * dT;
                }
            }
        }
    // } 

    // for (int i=1; i < nz-1; i++) {
        for (int j=1; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                if (!mesh.at(1).at(j).at(k).boundary) {
                    mesh.at(1).at(j).at(k).velocityX = tempVelocity.at(1).at(j).at(k).x;
                    mesh.at(1).at(j).at(k).velocityY = tempVelocity.at(1).at(j).at(k).y;
                    // mesh.at(1).at(j).at(k).velocityZ = tempU.at(1).at(j).at(k).z;
                }
            }
        }
    // }

}

Vector3 Cfd::getNetPressureOnPlane() {
    Vector3 netPressure = {0,0,0};

    for (int i=1; i < nz-1; i++) {
        for (int j=1; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                if (mesh.at(1).at(j).at(k).boundary) {
                    if (!mesh.at(i).at(j+1).at(k).boundary) {
                        netPressure.x += mesh.at(1).at(j+1).at(k).pressure * dy * dz;
                    }
                    if (!mesh.at(i).at(j-1).at(k).boundary) {
                        netPressure.x -= mesh.at(1).at(j-1).at(k).pressure * dy * dz;
                    }
                    if (!mesh.at(i).at(j).at(k+1).boundary) {
                        netPressure.y += mesh.at(1).at(j).at(k+1).pressure * dx * dz;
                    }
                    if (!mesh.at(i).at(j).at(k-1).boundary) {
                        netPressure.y -= mesh.at(1).at(j).at(k-1).pressure * dx * dz;
                    }
                    if (!mesh.at(i+1).at(j).at(k).boundary) {
                        netPressure.z += mesh.at(i+1).at(j).at(k).pressure * dx * dy;
                    }
                    if (!mesh.at(i-1).at(j).at(k).boundary) {
                        netPressure.z -= mesh.at(i-1).at(j).at(k).pressure * dx * dy;
                    }
                }
            }
        }
    }

    return netPressure;
}

Vector2 Cfd::calc(double anglePitch, double angleYaw)
{
    float cl, cd;
    double tijd = 0;
    while (tijd < maxTime)
    {
        tijd += dT;
        
        // TODO the movement of the pressure NOTE density is constant

        velocityMovement(dT);
        // removeDivergence();

        if (drawing) {
            Draw();
        }
        std::cout << tijd << " " << maxTime << std::endl;
    }

    // TODO correction fase
    // correction 

    Vector3 forces = getNetPressureOnPlane();
    // TODO the 100 is the starting velocity of the boudnary on the left
    cl = forces.y / (rho * pow(100 ,2) * 0.5);
    cd = forces.x / (rho * pow(100 ,2) * 0.5);
    // float cz = forces.z / (rho * pow(100 ,2) * 0.5);

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
    for (int i=1; i < nz-1; i+=(nz/30)) {
        for (int j=1; j < nx-1; j+=(nx/(nx/2))) {
            for (int k=1; k < ny-1; k+=(ny/(ny/4))) {
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
                    
                    double val = (velocity / 200.0f) *300;
                    double val2 = (velocity / 500.0f);
                    double val3 = velocity * 180;
                    
                    Color velocityColor = {255, val2, val3, 255};
                    
                    Vector3 velocityDirection = {velocityX,velocityY,velocityZ};
                    velocityDirection = Vector3Normalize2(velocityDirection);
                    velocityDirection.x = (velocityDirection.x * 0.5 * dx + point.x);
                    velocityDirection.y = (velocityDirection.y * 0.5 * dy + point.y);
                    velocityDirection.z = (velocityDirection.z * 0.5 * dz + point.z);
                    // std::cout << point.x << "  x " << velocityDirection.x << std::endl;
                    // std::cout << point.y << " y " << velocityDirection.y << std::endl;
                    // std::cout << point.z << " z " << velocityDirection.z << std::endl;
                    DrawLine3D(point, velocityDirection, velocityColor); //111
                    // DrawLine3D(point, {point.x, point.y, point.z+dz}, BLUE);
                    // DrawCubeWires(point, dx, dy, dz, RED);
                    // std::cout << velocity << " ";
                }
            }
            // std::cout  << std::endl;
            // std::cout << "velocity start " << mesh.at(1).at(0).at(0).velocityX << std::endl; 
        }
    }
            // std::cout  << std::endl;
            // std::cout  << std::endl;
            // std::cout  << std::endl;
}

void Cfd::draw2DGrid() {
    for (int i=1; i < 2; i++) {
        for (int j=0; j < nx-1; j++) {
            for (int k=1; k < ny-1; k++) {
                Vector3 point;
                point.x = j * dx - 0.5 * dx;
                point.y = k * dy - 0.5 * dy;
                point.z = startingPoint.z + i * dz - 0.5 * dz;

                if (mesh.at(i).at(j).at(k).boundary) {
                    DrawRectangle(point.x*4, point.y*4, dx*4, dy*4, BLACK);

                } else {
                    float velocityX = mesh.at(i).at(j).at(k).velocityX;
                    float velocityY = mesh.at(i).at(j).at(k).velocityY;
                    float velocityZ = mesh.at(i).at(j).at(k).velocityZ;
                    float velocity = sqrt(pow(velocityX,2) + pow(velocityY,2) + pow(velocityZ,2));
                    
                    double val = (velocity / 200.0f) *300;
                    double val2 = (velocity / 500.0f);
                    double val3 = velocity * 180;
                    
                    Color velocityColor = {255, val2, val3, 255};
                    DrawRectangle(point.x*4, point.y*4, dx*4, dy*4, velocityColor);
                    std::cout << velocity << " ";
                }
            }
            std::cout  << std::endl;
        }
    }
    std::cout  << std::endl;
    std::cout  << std::endl;
    std::cout  << std::endl;
}

void Cfd::Draw() {
    moveCamera(GetFrameTime());
    BeginDrawing();
        ClearBackground(WHITE);
        draw2DGrid();
        // BeginMode3D(camera);
        //     drawVelocityVectors();
        // EndMode3D();
    EndDrawing();
}

void Cfd::run(int steps, double stepsizePitch, double stepsizeYaw) { //333
    double stepsize = 360.0f/steps;
    std::vector<std::vector<Vector2>> cfdResults;
    for (double i=0; i < 360; i+=stepsize) { // pitch
        std::vector<Vector2> cfdResultsHelper;
        for (double j=0; j < 360; j+=stepsize) { // yaw
            airplane.transform = MatrixRotateXYZ2((Vector3){DEG2RAD * i, DEG2RAD * j, DEG2RAD * 0});
            resetMesh();
            // setPlaneBoundary();
                for (int j=1; j < ny-1; j++) {
                    mesh.at(1).at(nx/2).at(j).boundary = true;
                }
                    mesh.at(1).at(nx/2).at(ny/2-1).boundary = false;
                    mesh.at(1).at(nx/2).at(ny/2-2).boundary = false;
                    mesh.at(1).at(nx/2).at(ny/2+1).boundary = false;
                    mesh.at(1).at(nx/2).at(ny/2+2).boundary = false;

            Vector2 consts = calc(i, j);
            cfdResultsHelper.push_back({consts.x, consts.y});
        }
    }

    std::vector<Vector2> cfdResultsPitch, cfdResultsYaw;
    for (double i=0; i < 360; i+=stepsizePitch) { // pitch
        airplane.transform = MatrixRotateXYZ2((Vector3){DEG2RAD * i, DEG2RAD * 0, DEG2RAD * 0});
        resetMesh();
        setPlaneBoundary();
        Vector2 consts = calc(i, 0);
        cfdResultsPitch.push_back({consts.x, consts.y});
    }

    for (double i=0; i < 360; i+=stepsizeYaw) {
        airplane.transform = MatrixRotateXYZ2((Vector3){DEG2RAD * 0, DEG2RAD * i, DEG2RAD * 0});
        resetMesh();
        setPlaneBoundary();

        Vector2 consts = calc(0, i);
        cfdResultsYaw.push_back({consts.x, consts.y});
    }

    createLiftFiles(&cfdResults, &cfdResultsPitch, &cfdResultsYaw);
}

Cfd::Cfd(int setnx, int setny, int setnz, double deltaTime, double setMaxTime, double setRho, bool drawingEnabled)
{   
    Re = 100;
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
    setBoundaryConditions(10,  0,  0,  0,  0,  0);
}
 
Cfd::~Cfd()
{
}
