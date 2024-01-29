#pragma once
#include <raylib.h>
class physicsVector
{
private:
      /* data */
public:
      Vector3 components;
      Vector3 location;
      bool isAbsoluteForce;
      physicsVector(Vector3 setComponents, Vector3 setLocation);
      ~physicsVector();
};