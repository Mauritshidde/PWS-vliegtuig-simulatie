#pragma once
#include <raylib.h>
#include <vector>
#include "physicsvector.h"
class Plane
{
private:
      /* data */
      Model airplane;
      Texture2D airplaneTexture;

      float rotationMultiplier;

public:
      Plane(float givenMass = 10000, Vector3 startingPos = {0, 0, 0}, float givenrRotationMultiplier = 10);
      ~Plane();
      Vector3 calcLift();
      Vector3 calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth);

      void Start();
      void Draw();
      void Update(float deltaTime);

      float currentEngineTrust = 0.0f; // in newton
      float maxEngineTrust = 116000;   // in newton
      // float leftMotorForce;
      // float rightMotorForce;
      Vector3 pos, centerOfMass, centerOfLiftWingR, centerOfLiftWingL;
      Vector3 speedInDirections;
      float totalSpeed;
      float anglePitch, angleYaw, angleRoll;
      float mass; // in kg
};
