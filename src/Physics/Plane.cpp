#include "Plane.h"

Plane::Plane(float givenMass, Vector3 startingPos)
{
      mass = givenMass;
      speedInDirections = {0, 0, 0};
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;
      // centerOfLiftWingL = calcCenterOfLiftWing();
      // centerOfLiftWingR = calcCenterOfLiftWing();
}

Plane::~Plane()
{
}

Vector3 Plane::calcLift()
{
      // TODO lift formula
}

Vector3 Plane::calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth)
{
      // TODO lift formula
}