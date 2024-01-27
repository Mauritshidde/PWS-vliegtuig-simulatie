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