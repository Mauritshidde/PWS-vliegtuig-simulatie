#include "MeshCube.h"

MeshCube::MeshCube(float vx, float vy, float vz, float setPressure, float setDensity, bool bound, bool updatedPressure)
{
    velocityX = vx;
    velocityY = vy;
    velocityZ = vz;
    pressure = setPressure;
    density = setDensity;
    boundary = bound;
    updatedPressure = updatedPressure;
}

MeshCube::~MeshCube()
{
}