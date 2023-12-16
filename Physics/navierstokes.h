#pragma once
#include <vector>
#include <raylib.h>
#include <math.h>
#include <bits/stdc++.h>
#include "matrix.h"

#define WITHOUT_NUMPY
#include "matplotlibcpp.h"
namespace mat = matplotlibcpp;

class NavierStokes
{
private:
      int nx;   // amount of cells in x direction
      int ny;   // amount of cells in y direction
      int iMin; // minimaal 1
      int iMax;
      int jMin;
      int jMax;
      int Lx;
      int Ly;

      float uBottom, uTop, vLeft, vRight;
      float dx;
      float dy;
      float dxi;
      float dyi;
      float nu;
      float Re;
      float rho;
      float detOfL;
      float maxTime;
      float dT;

      std::vector<float> R, pv;
      std::vector<std::vector<float>> L, invL;
      std::vector<std::vector<float>> us, vs;
      std::vector<float> boundaryMinU, boundaryMaxU, boundaryMinV, boundaryMaxV;
      std::vector<std::vector<float>> v, u, p;

public:
      NavierStokes();
      ~NavierStokes();

      void generateVectors();
      void createMesh();
      void boundaryConditions();
      std::vector<std::vector<float>> calcL(std::vector<std::vector<float>> L);
      void plot();
      void calc();
};