#pragma once
#include <raylib.h>
#include "physicsvector.cpp"
#include <vector>
#include <math.h>

class Physics
{
private:
    /* data */
public:
    // add function to use 3d force vectors
    // force vector additions
    physicsVector translateForce(float force, Vector3 rotation); // create a vector3 force in the direction of rotation, don't know how yet
    physicsVector addForces(std::vector<physicsVector> inputVector);
    // physicsVector createPhysicsVector(std::vector<Vector3> force, std::vector<Vector3> location);
    float distanceBetweenPoints(Vector3 point1, Vector3 point2);
    float calcTorque(std::vector<physicsVector> forces, Vector3 centerOfMass);
    Physics(/* args */);
    ~Physics();
};

physicsVector addForces(std::vector<physicsVector> inputVector)
{
    physicsVector sumVector = physicsVector({0, 0, 0}, {0, 0, 0});
    for (int i = 0; i < inputVector.size(); i++)
    {
        sumVector.components.x += inputVector.at(i).components.x;
        sumVector.components.y += inputVector.at(i).components.y;
        sumVector.components.z += inputVector.at(i).components.z;
    }
    return sumVector;
}

float Physics::distanceBetweenPoints(Vector3 point1, Vector3 point2)
{
    float xComponent = pow(point1.x - point2.x, 2);
    float yComponent = pow(point1.y - point2.y, 2);
    float zComponent = pow(point1.z - point2.z, 2);
    float length = sqrt(xComponent + yComponent + zComponent);

    return length;
}

float Physics::calcTorque(std::vector<physicsVector> forces, Vector3 centerOfMass) // calcTorque for the x-axis
{
    Vector3 torqueSum = {0, 0, 0};
    for (int i = 0; i < forces.size(); i++)
    {
        float distance = distanceBetweenPoints(forces.at(i).location, centerOfMass);
        torqueSum.x += forces.at(i).components.x * distance;
        torqueSum.y += forces.at(i).components.y * distance;
        torqueSum.z += forces.at(i).components.z * distance;
        // calc distance/length between forces.at(i) and centerOfMass if forces.at(i) on the left of centerOfMass its negative
        // add this value to torqueSum
    }

    if (torqueSum.x > 0)
    {
        // keep goin straight forward/ without rotation
    }
    else if (torqueSum.x < 0)
    {
        // go right if positive and go left if negative didn't find formula yet
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
