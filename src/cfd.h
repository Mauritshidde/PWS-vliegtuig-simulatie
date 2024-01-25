#pragma once
#include <vector>
#include <raylib.h>
#include <math.h>
#include <bits/stdc++.h>
#include "Physics/MeshCube.h"
// #include "matrix.h"

class Cfd
{
private:
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

    void createMesh();
    void setBoundaryConditions(double velocityXDirectionStart, double velocityYDirectionStart, double velocityZDirectionStart, double velocityXDirectionEnd, double velocityYDirectionEnd, double velocityZDirectionEnd);
    void setPlaneBoundary(); // make parts of the plane part of the boundary conditions
    void iterativeSolver(double density);
    void densityDispersion();
    void removeDivergence();
    
    void solveDensity(int i, int j, int k);
    void solveDensityFirst(int i, int j, int k);
    void solvePressure(int i, int j, int k);
    void solvePressureFirst(int i, int j, int k);
public:
    void calc();    
    void Draw();
    
    Cfd(int setnx = 100, int setny = 100, int setnz = 100, double deltaTime = 0.1, double setMaxTime = 1000);
    ~Cfd();
};