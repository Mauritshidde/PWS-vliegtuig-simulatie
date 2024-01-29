#pragma once
#include <vector>
#include <raylib.h>
#include <math.h>
#include <bits/stdc++.h>
#include "Physics/MeshCube.h"
#include "Physics/ModelLoader.h"
#include "extra/raymath2.h"
#include "ui/loadingScreen.h"
// #include "matrix.h"

class Cfd
{
private:
    // camera variables
    Vector2 cameraYZPos;
    Vector3 cameraPos;
    Vector2 cameraXYPos;
    Vector2 previousMousePosition;
    Camera camera;
    float cameraCircleRadius;
    float angleYAxis;
    float angleXZAxis;

    // drawing variables
    Model airplane;
    Texture airplaneTexture;
    bool drawing;
    BoundingBox boundingBoxPlane;
    Vector3 boundingBoxPlaneMin;
    Vector3 boundingBoxPlaneMax;

    // multi threading
    int cores;
    bool settingPlaneBOundarys;

    // mesh variables
    int nx;   // amount of cells in x direction // steps in x direction
    int ny;   // amount of cells in y direction // steps in y direction
    int nz;
    double dx; // step size in x
    double dy; // step size in y
    double dz; // step size in z
    double dxi; // 1/dx
    double dyi; // 1/dy
    double dzi; // 1/dz
    Vector3 startingPoint; // starting point of the grid
    FluidDynamicsModel plane;

    // calculation consts
    double k; // const for amount of change of density in each time step
    double nu; // const
    double Re; // reynolds number const
    double rho; // density const
    double maxTime; // max time for the program untill the program should quit
    double dT; // time steps

    std::vector<std::vector<std::vector<MeshCube>>> mesh;
    std::vector<std::vector<std::vector<double>>> divergenceVelocityScalarField;
    std::vector<std::vector<std::vector<Vector3>>> gradientPressureField, divergenceVelocityField, divergenceFreeField;

    // functions for the creation of the grid
    void Start();
    void createMesh();
    void resetMesh();
    void setBoundaryConditions(double velocityXDirectionStart, double velocityYDirectionStart, double velocityZDirectionStart, double velocityXDirectionEnd, double velocityYDirectionEnd, double velocityZDirectionEnd);

    // functions for setting the plane boundary
    void setPlaneBoundaryHelper(int startIndex, int endIndex);
    bool getCollisionPlaneRay(Vector3 direction, Vector3 oppositeDirection, Ray ray, Ray ray2);
    void detectColission();
    void setPlaneBoundary(); // make parts of the plane part of the boundary conditions 

    // functions for calculating the movement of the fluid
    void removeDivergence();
    void velocityMovement(float dT);
    void solvePressure(int i, int j, int k);
    void solvePressureFirst(int i, int j, int k);
    Vector2 calc(double anglePitch, double angleYaw);    

    // graphics functions (these are optionally when running the cfd) 
    void moveCamera(float deltaTime);
    void drawVelocityVectors();
    void Draw();
public:
    void run(int steps);

    Cfd(int setnx = 20, int setny = 10, int setnz = 60, double deltaTime = 0.1, double setMaxTime = 1000, double setRho = 1.293, bool drawingEnabled = true);
    ~Cfd();
};