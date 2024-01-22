#pragma once
#include <vector>
#include <raylib.h>
#include <math.h>
#include <bits/stdc++.h>
#include "MeshCube.h"
// #include "matrix.h"

class Cfd
{
private:
    // mesh variables
    int nx;   // amount of cells in x direction // steps in x direction
    int ny;   // amount of cells in y direction // steps in y direction
    int nz;
    float dx; // step size in x
    float dy; // step size in y
    float dz; // step size in z
    float dxi; // 1/dx
    float dyi; // 1/dy
    float dzi; // 1/dz
    
    // calculation consts
    float nu; // const
    float Re; // reynolds number const
    float rho; // density const
    float maxTime; // max time for the program untill the program should quit
    float dT; // time steps

    std::vector<std::vector<std::vector<MeshCube>>> mesh;

    void createMesh();
    void setBoundaryConditions(float velocityXDirectionStart, float velocityYDirectionStart, float velocityZDirectionStart, float velocityXDirectionEnd, float velocityYDirectionEnd, float velocityZDirectionEnd);
    void setPlaneBoundary(); // make parts of the plane part of the boundary conditions
public:
    void calc();    

    Cfd(int setnx = 100, int setny = 100, int setnz = 100);
    ~Cfd();
};