#include "MeshCube.h"

MeshCube::MeshCube(double vx, double vy, double vz, double setPressure, double setDensity, bool bound, bool updatedPressure)
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