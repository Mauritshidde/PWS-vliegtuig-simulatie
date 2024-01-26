#include "Physics.h"

physicsVector Physics::addForces(std::vector<physicsVector> inputVector)
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

Vector3 Physics::vectorAddition(Vector3 vec1, Vector3 vec2)
{
    Vector3 result;
    result.x = vec1.x + vec2.x;
    result.y = vec1.y + vec2.y;
    result.z = vec1.z + vec2.z;
    return result;
}

Vector3 Physics::vectorSubtraction(Vector3 vec1, Vector3 vec2)
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

Vector3 Physics::calcTorque(std::vector<physicsVector> forces, Vector3 centerOfMass)
{
    Vector3 torqueTotal = {0, 0, 0}; // resulting torque for pitch yaw and roll
    float torquePitch, torqueYaw, torqueRoll;
    Vector3 distance;
    Vector3 currentForce;
    for (int i = 0; i < forces.size(); i++)
    {
        currentForce = forces.at(i).components;
        distance = vectorSubtraction(centerOfMass, forces.at(i).location);

        torquePitch = distance.y * currentForce.y + distance.z * currentForce.z;
          torqueYaw = distance.x * currentForce.x + distance.z * currentForce.z;
         torqueRoll = distance.x * currentForce.x + distance.y * currentForce.y;

        torqueTotal.x += torquePitch;
        torqueTotal.y += torqueYaw;
        torqueTotal.z += torqueRoll;
    }
    return torqueTotal;
}

Vector3 Physics::calcAngularAcceleration(std::vector<physicsVector> forces, float mass, Vector3 centerOfMass, Vector3 momentOfInertia)
{
    Vector3 angularAcceleration, torque;
    torque = calcTorque(forces, centerOfMass);
    std::cout << "torque xyz" << torque.x << " " << torque.y << " " << torque.z << " \n";
    std::cout << "moment xyz" << momentOfInertia.x << " " << momentOfInertia.y << " " << momentOfInertia.z << " \n";
    angularAcceleration.x = torque.x / momentOfInertia.x;
    angularAcceleration.y = torque.y / momentOfInertia.y;
    angularAcceleration.z = torque.z / momentOfInertia.z;
    return angularAcceleration;
}

float Physics::calcHypot(Vector3 components)
{
    return sqrt(pow(components.x, 2) + pow(components.y, 2) + pow(components.z, 2));
}

Vector3 Physics::calcAcceleration(std::vector<physicsVector> forces, float mass) // newton: F = m・a ------ a = F/m
{
    Vector3 acceleration;
    for (int i = 0; i < forces.size(); i++)
    {
        acceleration.x += forces.at(i).components.x / mass;
        acceleration.y += forces.at(i).components.y / mass;
        acceleration.z += forces.at(i).components.z / mass;
    }
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

Vector3 Physics::calcForceGravity(float mass)
{
    Vector3 ForceG = {0, 0, 0};
    ForceG.y = -mass * gAcceleration;
    return ForceG;
}

Physics::Physics(/* args */)
{
}

Physics::~Physics()
{
}
