#include "cfd.h"
#include <cuda_runtime.h>

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
        tempVelocity.push_back(helper3);
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
            // mesh.at(i).at(0).at(k).pressure = pow(velocityXDirectionStart,2) * (rho/2.0f); // quess for starting pressure

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

__device__ bool getCollisionPlaneRay(Vector3 direction, Vector3 oppositeDirection, Ray ray, Ray ray2) {
    ray.direction = direction;
    ray2.direction = oppositeDirection;

    // RayCollision meshHitInfo = GetRayCollisionMesh(ray, *airplane.meshes, airplane.transform);
    // RayCollision meshHitInfo2 = GetRayCollisionMesh(ray2, *airplane.meshes, airplane.transform);
    
    // if (meshHitInfo.hit && meshHitInfo2.hit) {
    //     return true;
    // } else {
    //     return false;
    // }
}

__global__ void setPlaneBoundaryHelper(MeshCube *mesh, Vector3 startingPoint, int dx, int dy, int dz, int i, int j) {
    Vector3 position;
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    position.x = dx * j + startingPoint.x;
    position.y = dy * k + startingPoint.y;
    position.z = dz * i + startingPoint.z;
    // if (CheckCollisionBoxSphere(boundingBoxPlane, position, dx)) {
        Ray ray;
        Ray ray2;
        
        ray.position = position;
        ray2.position = position;
        
        // check if the cube is inside the plane
        if (getCollisionPlaneRay({1,0,0}, {1,0,0}, ray, ray2)) { 
            if (getCollisionPlaneRay({0,1,0}, {0,1,0}, ray, ray2)) {
                if (getCollisionPlaneRay({0,0,1}, {0,0,1}, ray, ray2)) {
                    mesh[k].boundary = true;
                }
            }
        }
    // }
}

void Cfd::setPlaneBoundary() //222
{
    int N2 = mesh.at(0).at(0).size();
    // int threadsPerBlock = 256;
    // int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    size_t size = sizeof(MeshCube) * N2;
    for (int i=1; i < nz-1; i++) {
        for (int j=1; j < nx-1; j++) {
        
            MeshCube *array = mesh.at(i).at(j).data();
            MeshCube *array_p;

            cudaMalloc(&array_p, size);
            cudaMemcpy(array_p, array, size, cudaMemcpyHostToDevice);
            
            setPlaneBoundaryHelper<<<grid_size, block_size>>>(array_p, startingPoint, dx, dy, dz, i, j);

            cudaMemcpy(array, array_p, size, cudaMemcpyDeviceToHost);
        }
    }
    std::cout << "done this" << std::endl;
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

__global__ void velocityMovementHelper(int N, int M, int B, MeshCube *mesh, float dT, double rho, double dx, double dy, double dz) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    int k = index / (M * B);
    int j = (index % (M * B)) / B;
    int i = (index % B);

    if (i < N * M * B && (i > 0 && i < N-1) && (j > 0 && j < M-1) && (k > 0 && k < B-1)) {
        if (!mesh[index].boundary) {
            // int index = z * (nx * ny) + x * ny + y;
            double vx = mesh[index].velocityX;
            double vy = mesh[index].velocityY;
            double vz = mesh[index].velocityZ;
            double duDt = -(vx * (vx - mesh[k + i*M * B + (j-1) * B].velocityX) / dx +
                    vy * (vx - mesh[(k-1) + i*M * B + j * B].velocityX) / dy + 
                    vz * (vx - mesh[k + (i-1)*M * B + j * B].velocityX) / dz) / dx;

            double dvDt = -(vx * (vy - mesh[k + i*M * B + (j-1) * B].velocityY) / dx +
                    vy * (vy - mesh[(k-1) + i*M * B + j * B].velocityY) / dy + 
                    vz * (vy - mesh[k + (i-1)*M * B + j * B].velocityY) / dz) / dy;

            double dwDt = -(vx * (vz - mesh[k + i*M * B + (j-1) * B].velocityZ) / dx +
                    vy * (vz - mesh[(k-1) + i*M * B + j * B].velocityZ) / dy +
                    vz * (vz - mesh[k + (i-1)*M * B + j * B].velocityZ) / dz) / dz;

            mesh[index].tempVelocity.x = mesh[index].velocityX + duDt * dT;
            mesh[index].tempVelocity.y = mesh[index].velocityY + dvDt * dT;
            mesh[index].tempVelocity.z = mesh[index].velocityZ + dwDt * dT;
            mesh[index].pressure = (duDt / dT + dvDt / dT + dwDt / dT) * rho; 
        }
    }
}

__global__ void velocityMovementUpdater(int N, int M, int B, MeshCube *mesh) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    int k = index / (M * B);
    int j = (index % (M * B)) / B;
    int i = (index % B);

    if (i < N * M * B && (i > 0 && i < N-1) && (j > 0 && j < M-1) && (k > 0 && k < B-1)) {
        // int index = z * (nx * ny) + x * ny + y;
        if (!mesh[index].boundary) {
            mesh[k + i*M * B + j * B].velocityX = mesh[k + i*M * B + j * B].tempVelocity.x;
            mesh[k + i*M * B + j * B].velocityX = mesh[k + i*M * B + j * B].tempVelocity.y;
            mesh[k + i*M * B + j * B].velocityX = mesh[k + i*M * B + j * B].tempVelocity.z; 
        }
    }
}

// void Cfd::velocityMovement(float dT) {
//     // std::cout << "start " << std::endl;
//     // int N = mesh.at(0).at(0).size();
//     // int M = mesh.at(0).size();
//     // int B = mesh.size();
    

//     // auto start = std::chrono::system_clock::now();
//     velocityMovementHelper<<<grid_size, block_size>>>(N, M, B, array_p, dT, rho, dx, dy, dz);
//     // auto end = std::chrono::system_clock::now();
//     // std::chrono::duration<double> elapsed_seconds = end-start;
//     // std::cout << elapsed_seconds.count() << std::endl;
// }

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
    std::vector<std::vector<std::vector<double>>> *diffuseV;
    
    for (int i=0; i < N; i++) {
        
        for (int j=0; j < M; j++) {
            for (int k=0; k < B; k++) {
                mesh_array[k + i*M * B + j * B] = mesh.at(i).at(j).at(k);
            }
        }
    } 

    std::cout << "start loop for getting pressure and velocity" << std::endl;
    // setting gpu memory
    MeshCube *array_p;

    cudaMalloc(&array_p, N * M * B * sizeof(MeshCube));
    cudaMemcpy(array_p, mesh_array, N * M * B * sizeof(MeshCube), cudaMemcpyHostToDevice);
    
    while (tijd < maxTime)
    {
        auto start = std::chrono::system_clock::now();
        tijd += dT;
        
        // TODO the movement of the pressure NOTE density is constant

        std::vector<std::thread> threads;
        int newNz = nz - 2;

        // std::cout << "start movement" << std::endl;
        // velocityMovement(dT);

        velocityMovementHelper<<<grid_size, block_size>>>(N, M, B, array_p, dT, rho, dx, dy, dz);
        velocityMovementUpdater<<<grid_size, block_size>>>(N, M, B, array_p);
        


        // std::cout << "end movement" << std::endl;

        // for (int i=1; i < nz-1; i++) {
        //     for (int j=1; j < nx-1; j++) {
        //         for (int k=1; k < ny-1; k++) {
        //             if (!mesh_array[k + i*M * B + j * B].boundary) {
        //                 mesh_array[k + i*M * B + j * B].velocityX = mesh_array[k + i*M * B + j * B].tempVelocity.x;
        //                 mesh_array[k + i*M * B + j * B].velocityX = mesh_array[k + i*M * B + j * B].tempVelocity.y;
        //                 mesh_array[k + i*M * B + j * B].velocityX = mesh_array[k + i*M * B + j * B].tempVelocity.z;
        //             }
        //         }
        //     }
        // }

        if (drawing) {
            // for (int i=0; i < N; i++) {
            //     for (int j=0; j < M; j++) {
            //         for (int k=0; k < B; k++) {
            //             mesh.at(i).at(j).at(k) = mesh_array[k + i*M * B + j * B];
            //         }
            //     }
            // } 
            Draw();
        }
        
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << elapsed_seconds.count() << std::endl;
        std::cout << tijd << std::endl;
    }

    // deleting array_p and writing gpu memory to cpu
    cudaMemcpy(mesh_array, array_p, N * M * B * sizeof(MeshCube), cudaMemcpyDeviceToHost);
    cudaFree(array_p);

    std::cout << "done with loop getting pressure and velocity" << maxTime << std::endl;
    for (int i=0; i < N; i++) {
        for (int j=0; j < M; j++) {
            for (int k=0; k < B; k++) {
                mesh.at(i).at(j).at(k) = mesh_array[k + i*M * B + j * B];
            }
        }
    } 
    // TODO correction fase
    // correction 

    Vector3 forces = getNetPressureOnPlane();
    // TODO the 100 is the starting velocity of the boudnary on the left
    cl = forces.y / (rho * pow(10 ,2) * 0.5);
    cd = forces.x / (rho * pow(10 ,2) * 0.5);
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
                    
                    Color velocityColor = {val, val2, val3, 255};
                    
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
    for (int j=0; j < nx-1; j++) {
        for (int k=1; k < ny-1; k++) {
            Vector3 point;
            point.x = j * dx - 0.5 * dx;
            point.y = k * dy - 0.5 * dy;
            point.z = startingPoint.z + dz - 0.5 * dz;

            if (mesh_array[k + 1*M * B + j * B].boundary) {
                DrawRectangle(point.x*4, point.y*4, dx*4, dy*4, BLACK);

            } else {
                float velocityX = mesh_array[k + 1*M * B + j * B].velocityX;
                float velocityY = mesh_array[k + 1*M * B + j * B].velocityY;
                float velocityZ = mesh_array[k + 1*M * B + j * B].velocityZ;
                float velocity = sqrt(pow(velocityX,2) + pow(velocityY,2) + pow(velocityZ,2));
                
                double val = (velocity / 200.0f) *300;
                double val2 = (velocity / 500.0f);
                double val3 = velocity * 180;
                
                Color velocityColor = {255, val2, val3, 255};
                DrawRectangle(point.x*4, point.y*4, dx*4, dy*4, velocityColor);
                // std::cout << mesh.at(1).at(j).at(k).pressure << " ";
                // std::cout << velocity << " ";
            }
        }
        // std::cout  << std::endl;
    }
    // std::cout  << std::endl;
    // std::cout  << std::endl;
    // std::cout  << std::endl;
}

void Cfd::Draw() {
    moveCamera(GetFrameTime());
    BeginDrawing();
        ClearBackground(WHITE);
        if (drawing3D) {
            BeginMode3D(camera);
                drawVelocityVectors();
            EndMode3D();
        } else {
            draw2DGrid();
        }
    EndDrawing();
}

void Cfd::run(int steps, double stepsizePitch, double stepsizeYaw) { //333
    N = mesh.at(0).at(0).size();
    M = mesh.at(0).size();
    B = mesh.size();
    
    grid_size = ((N * M * B + 255) / block_size);
    mesh_array = (MeshCube*)malloc(N * M * B * sizeof(MeshCube));

    double stepsize = 360.0f/steps;
    std::vector<std::vector<Vector2>> cfdResults;
    for (double i=0; i <= 360; i+=stepsize) { // pitch
        std::vector<Vector2> cfdResultsHelper;
        for (double j=0; j <= 360; j+=stepsize) { // yaw
            airplane.transform = MatrixRotateXYZ2((Vector3){DEG2RAD * i, DEG2RAD * j, DEG2RAD * 0});
            resetMesh();
            // setPlaneBoundary();
                // for (int j=1; j < ny-1; j++) {
                //     mesh.at(1).at(nx/2).at(j).boundary = true;
                // }
                //     mesh.at(1).at(nx/2).at(ny/2-1).boundary = false;
                //     mesh.at(1).at(nx/2).at(ny/2-2).boundary = false;
                //     mesh.at(1).at(nx/2).at(ny/2+1).boundary = false;
                //     mesh.at(1).at(nx/2).at(ny/2+2).boundary = false;

            Vector2 consts = calc(i, j);
            cfdResultsHelper.push_back(consts);
        }
    }
    std::cout << "done" << std::endl;
    std::vector<Vector2> cfdResultsPitch, cfdResultsYaw;
    for (double i=0; i <= 360; i+=stepsizePitch) { // pitch
        airplane.transform = MatrixRotateXYZ2((Vector3){DEG2RAD * i, DEG2RAD * 0, DEG2RAD * 0});
        resetMesh();
        setPlaneBoundary();
        Vector2 consts = calc(i, 0);
        cfdResultsPitch.push_back({consts.x, consts.y});
    }
    std::cout << "done2" << std::endl;

    for (double i=0; i <= 360; i+=stepsizeYaw) {
        airplane.transform = MatrixRotateXYZ2((Vector3){DEG2RAD * 0, DEG2RAD * i, DEG2RAD * 0});
        resetMesh();
        setPlaneBoundary();

        Vector2 consts = calc(0, i);
        cfdResultsYaw.push_back({consts.x, consts.y});
    }
    std::cout << "done3" << std::endl;

    createLiftFiles(&cfdResults, &cfdResultsPitch, &cfdResultsYaw);
    if (drawing) {
        CloseWindow();
    }
    
    free(mesh_array);
    std::cout << "cfd-program completed calculating cl and cd over pitch and yaw and exited succesfully";
}

Cfd::Cfd(int setnx, int setny, int setnz, double deltaTime, double setMaxTime, double setRho, bool drawingEnabled, bool draw3D)
{   
    Re = 100;
    nu = 1 / Re;

    // set variables for gpu
    block_size = 2048;

    // set multithreading variables
    cores = 12;
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
    dx = 2;
    dy = 2;
    dz = 5;
    startingPoint.x = -(nx*dx)/2;
    startingPoint.y = -(ny*dy)/2 + 10;
    startingPoint.z = -(nz*dz)/2;
    
    drawing3D = draw3D;
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
