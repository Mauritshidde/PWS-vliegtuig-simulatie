#pragma once
#include <raylib.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

class FluidDynamicsModel
{
private:
    std::vector<Vector3> v, vn;
    std::vector<Vector2> vt;
    std::vector<std::vector<std::string>> f;
    std::vector<std::vector<Vector3>> objectTrianglePoints;

    void readObjectModelFile(std::string *as, std::string *bs, std::string *cs, float *a, float *b, float *c);

public:
    FluidDynamicsModel();
    ~FluidDynamicsModel();

    int detectCollision(Ray ray);
    void loadObjectModel();
    void drawModel();

    void modelToCubes();
};