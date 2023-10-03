#include <raylib.h>


class Plane {
private:
      /* data */
public:
      Vector3 pos;
      float speed;
      float anglePitch, angleYaw, angleRoll;
      float mass;
      
      
      Plane();
      ~Plane();
};

Plane::Plane() {
      speed = 0;
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;
}