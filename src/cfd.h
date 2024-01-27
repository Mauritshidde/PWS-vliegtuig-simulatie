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

    // multi threading
    int cores;
    bool done;
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
    std::vector<std::vector<std::vector<double>>> divergenceVelocityField;
    std::vector<std::vector<std::vector<Vector3>>> gradientPressureField;

    // functions for the creation of the grid
    void Start();
    void createMesh();
    void resetMesh();
    void setBoundaryConditions(double velocityXDirectionStart, double velocityYDirectionStart, double velocityZDirectionStart, double velocityXDirectionEnd, double velocityYDirectionEnd, double velocityZDirectionEnd);

    // functions for setting the plane boundary
    void setPlaneBoundaryHelper(int startIndex, int endIndex);
    void detectColission();
    void setPlaneBoundary(); // make parts of the plane part of the boundary conditions 

    // functions for calculating the movement of the fluid
    void densityDispersion();
    void removeDivergence();
    void solveDensity(int i, int j, int k);
    void solveDensityFirst(int i, int j, int k);
    void solvePressure(int i, int j, int k);
    void solvePressureFirst(int i, int j, int k);
    void calc(double anglePitch, double angleYaw);    

    // graphics functions (these are optionally when running the cfd) 
    void moveCamera();
    void Draw();
public:
    void run(int steps);

    Cfd(int setnx = 90, int setny = 60, int setnz = 80, double deltaTime = 0.1, double setMaxTime = 1000, double setRho = 1.293);
    ~Cfd();
};