#pragma once
#include <raylib.h>
#include "physicsvector.h"
#include <vector>
#include <math.h>
#define g 9.81

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
    Vector3 crossProduct(Vector3 vec1, Vector3 vec2);
    float distanceBetweenPoints(Vector3 point1, Vector3 point2);
    float calcTorque(std::vector<physicsVector> forces, Vector3 centerOfMass);
    float calcHypot(Vector3 components);
    Vector3 calcAcceleration(physicsVector force, float mass);
    Vector3 calcDeltaV(float deltaTime, Vector3 acceleration);
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

Vector3 vectorAddition(Vector3 vec1, Vector3 vec2)
{
    Vector3 result;
    result.x = vec1.x + vec2.x;
    result.y = vec1.y + vec2.y;
    result.z = vec1.z + vec2.z;
    return result;
}

Vector3 vectorSubtraction(Vector3 vec1, Vector3 vec2)
{
    Vector3 result;
    result.x = vec1.x - vec2.x;
    result.y = vec1.y - vec2.y;
    result.z = vec1.z - vec2.z;
    return result;
}

Vector3 Physics::crossProduct(Vector3 vec1, Vector3 vec2)
{
    Vector3 result;

    result.x = vec1.y * vec2.z - vec1.z * vec2.y;
    result.y = vec1.z * vec2.x - vec1.x * vec2.z;
    result.z = vec1.x * vec2.y - vec1.y * vec2.x;

    return result;
}

float Physics::distanceBetweenPoints(Vector3 point1, Vector3 point2)
{
    float xComponent = pow(point1.x - point2.x, 2);
    float yComponent = pow(point1.y - point2.y, 2);
    float zComponent = pow(point1.z - point2.z, 2);
    float length = sqrt(xComponent + yComponent + zComponent);

    return length;
}

float Physics::calcTorque(std::vector<physicsVector> forces, Vector3 centerOfMass)
{
    Vector3 torqueSum = {0, 0, 0};
    for (int i = 0; i < forces.size(); i++)
    {
        Vector3 distance = vectorSubtraction(centerOfMass, forces.at(i).location);
        Vector3 extraTorque = crossProduct(forces.at(i).components, distance);
        torqueSum = vectorAddition(torqueSum, extraTorque);
    }
    // for rpm rotatie ook inersia nodig
    // if (torqueSum.x > 0)
    // {
    //     // keep goin straight forward/ without rotation
    // }
    // else if (torqueSum.x < 0)
    // {
    //     // go right if positive and go left if negative didn't find formula yet
    // }
}

float Physics::calcHypot(Vector3 components)
{
    return sqrt(pow(components.x, 2) + pow(components.y, 2) + pow(components.z, 2));
}

Vector3 Physics::calcAcceleration(physicsVector force, float mass) // newton: F = m・a ------ a = F/m
{
    Vector3 acceleration;
    acceleration.x = force.components.x / mass;
    acceleration.y = force.components.y / mass;
    acceleration.z = force.components.z / mass;
    return acceleration;
}

Vector3 Physics::calcDeltaV(float deltaTime, Vector3 acceleration)
{
    Vector3 deltaV;
    deltaV.x = acceleration.x * deltaTime;
    deltaV.y = acceleration.y * deltaTime;
    deltaV.z = acceleration.z * deltaTime;
    return deltaV;
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
