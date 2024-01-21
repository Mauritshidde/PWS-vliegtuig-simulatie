#include "MeshCube.h"

MeshCube::MeshCube(float vx, float vy, float vz, float setPressure, bool bound)
{
    velocityXDirection = vx;
    velocityYDirection = vy;
    velocityZDirection = vz;
    pressure = setPressure;
    boundary = bound;
}

MeshCube::~MeshCube()
{
}