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
