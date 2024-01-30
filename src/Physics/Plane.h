#pragma once
#include <raylib.h>
#include <vector>
#include <string>

#include "Physics/Physics.h"
#include "physicsvector.h"
#include "../liftFileCode/readFile.h"

class Plane
{
private:
      Model airplane;
      Texture2D airplaneTexture;
      LiftFileReader files;

      float rotationMultiplier;
      std::string liftFileName;

public:
      Plane(std::string fileName = "Boeing737", float startVelocity = 100, float rho = 1.225);
      ~Plane();
      void Draw();
      void Update(float deltaTime, float rho);
      Vector3 calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth);
      Vector2 getConsts(float pitch, float yaw, bool useYaw, bool usePitch);

      void calcLift(float rho);
      void evaluateForces(std::vector<physicsVector> forces);
      void updateVel(float deltaTime);
      void updateAngularVel(float deltaTime);
      void updateRotation(float deltaTime);

      void reduceAngleDegrees();
      Vector3 reduceAngleDegrees(Vector3 angle);
      void rotatePoints();
      void rotateVector();


      Physics planePhysics;

      float currentEngineTrust = 0.0f; // in newton
      float maxEngineTrust = -116000;   // in newton
      float engineOffset; //distance of the engine to the center of mass
      Vector2 consts;
      float cl;
      float cd;
      float lift;
      float drag;
      float leftMotorForce;
      float rightMotorForce;
      Vector3 leftMotorDirectionPoint, rightMotorDirectionPoint;
      Vector3 externalPos, pos, centerOfMass, centerOfLiftWingR, centerOfLiftWingL;
      Vector3 velocity, angularVelocity, acceleration, angularAcceleration;
      Vector3 leftMotorThrustDirection, leftMotorThrust,rightMotorThrust;
      Vector3 liftDirection, liftForce, dragDirection, dragForce;
      std::vector<physicsVector> forces;
      std::vector<physicsVector*> internalForces, relativeForces;
      
      Vector3 momentOfInertia;

      float speed;         // in m/s
      float wingArea;         // surface area of wing in m2
      float planeFrontalArea; // the surface area of the face of the plane that is parellel to the velocity direction
      // float totalSpeed;
      float anglePitch, angleYaw, angleRoll;
      float previousAnglePitch, previousAngleYaw, previousAngleRoll;
      float mass; // in kg
      float height;
};
