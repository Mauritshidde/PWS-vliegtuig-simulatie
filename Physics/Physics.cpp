#pragma once
#include <raylib.h>
#include "forceVector.cpp"
#include <vector>
#include <math.h>

class Physics
{
private:
    /* data */
public:
    // add function to use 3d force vectors
    // force vector additions
    forceVector addForces(std::vector<forceVector> inputVector);
    forceVector createForceVector(std::vector<Vector3> force, std::vector<Vector3> location);
    float distanceBetweenPoints(Vector3 point1, Vector3 point2);
    float calcTorque(std::vector<forceVector> forces, Vector3 centerOfMass);
    Physics(/* args */);
    ~Physics();
};

forceVector addForces(std::vector<forceVector> inputVector)
{
    forceVector sumVector = forceVector();
    for (int i = 0; i < inputVector.size(); i++)
    {
        sumVector.force.x += inputVector.at(i).force.x;
        sumVector.force.y += inputVector.at(i).force.y;
        sumVector.force.z += inputVector.at(i).force.z;
    }
    return sumVector;
}

forceVector Physics::createForceVector(std::vector<Vector3> force, std::vector<Vector3> location)
{
}

float Physics::distanceBetweenPoints(Vector3 point1, Vector3 point2)
{
    float xComponent = pow(point1.x - point2.x, 2);
    float yComponent = pow(point1.y - point2.y, 2);
    float zComponent = pow(point1.z - point2.z, 2);
    float lenght = sqrt(xComponent + yComponent + zComponent);

    return lenght;
}

float Physics::calcTorque(std::vector<forceVector> forces, Vector3 centerOfMass) // calcTorque for the x-axis
{
    float torqueSum = 0;
    for (int i = 0; i < forces.size(); i++)
    {
        float distance = distanceBetweenPoints(forces.at(i).location, centerOfMass);
        torqueSum += distance;
        // calc distance/lenght between forces.at(i) and centerOfMass if forces.at(i) on the left of centerOfMass its negative
        // add this value to torqueSum
    }

    if (torqueSum > 0)
    {
        // keep goin straight forward/ without rotation
    }
    else if (torqueSum < 0)
    {
        // go right if positive and go left if negatife didn't find formula yet
    }
}

Physics::Physics(/* args */)
{
}

Physics::~Physics()
{
}

int main()
{

    return 0;
}
