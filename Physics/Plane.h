#include <raylib.h>


class Plane {
private:
      /* data */
public:
      Plane(float givenMass, Vector3 startingPos);
      ~Plane();
      Vector3 calcLift();
      Vector3 calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth);
      
      Vector3 pos, centerOfMass, centerOfLiftWingR, centerOfLiftWingL;
      float speed;
      float anglePitch, angleYaw, angleRoll;
      float mass;
};

Plane::Plane(float givenMass, Vector3 startingPos) {
      mass = givenMass;
      speed = 0;
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;
      // centerOfLiftWingL = calcCenterOfLiftWing();
      // centerOfLiftWingR = calcCenterOfLiftWing();
}

Plane::~Plane() {

}

Vector3 Plane::calcLift() {
      // TODO lift formula
}

Vector3 Plane::calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth) {
      // TODO lift formula
}