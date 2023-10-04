#pragma once
#include <raylib.h>
#include "forceVector.cpp"
#include <vector>
class Physics
{
private:
    /* data */
public:
    // add function to use 3d force vectors
    // force vector additions
    forceVector addForces(std::vector<forceVector> inputVector);
    forceVector createForceVector(std::vector<Vector3> force, std::vector<Vector3> location);
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

float Physics::calcTorque(std::vector<forceVector> forces, Vector3 centerOfMass) // calcTorque for the x-axis
{
    float torqueSum = 0;
    for (int i=0; i < forces.size(); i++) 
    {
        // calc distance/lenght between forces.at(i) and centerOfMass if forces.at(i) on the left of centerOfMass its negative
        // add this value to torqueSum
    }

    if (torqueSum == 0) {
        //keep goin straight forward/ without rotation
    } else 
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
