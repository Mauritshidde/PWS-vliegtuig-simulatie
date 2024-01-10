#include "Plane.h"

#define RAYMATH_IMPLEMENTATION
#include "../../include/modules/raymath.h"

Plane::Plane(float givenMass, Vector3 startingPos, float givenrRotationMultiplier)
{
      mass = givenMass;
      speedInDirections = {0, 0, 0};
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;
      rotationMultiplier = givenrRotationMultiplier;
      angleUpdated = false;
      Vector2 consts = getConsts(anglePitch, angleYaw, true, true);
      
      cl = consts.x;
      cd = consts.y;
      // centerOfLiftWingL = calcCenterOfLiftWing();
      // centerOfLiftWingR = calcCenterOfLiftWing();
}

Plane::~Plane()
{
}

Vector2 Plane::getConsts(float pitch, float yaw, bool usePitch, bool useYaw) {
      if (!usePitch && !useYaw) {
            usePitch = true;
            useYaw = true;
      }
      Vector2 consts = getConstFromLiftFile(pitch, yaw, useYaw, usePitch);
      return consts;
}

Vector3 Plane::calcLift()
{
      // lift = cl * rho * pow(velocity, 2) * wingArea * 0.5;
      //Drag = cd * rho * pow(velocity, 2) * planeArea * 0.5
      // TODO lift formula
}

Vector3 Plane::calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth)
{
      // TODO lift formula
}

void Plane::Start()
{
      airplaneTexture = LoadTexture("models/texture/skyboxtexture.png");
      airplane = LoadModel("models/object/plane.obj");
      // airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;
}

void Plane::Draw()
{
      DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f}, (Vector3){1.0f, 0.0f, 0.0f}, 0, (Vector3){0.5f, 0.5f, 0.5f}, RED); // 2de vector geeft aan met welke factor hij met currentangle draait
}

void Plane::Update(float deltaTime)
{
      if (IsKeyDown(KEY_W))
      {
            anglePitch += rotationMultiplier * deltaTime;
            if (anglePitch > 360)
            {
                  anglePitch -= 360;
            }
            angleUpdated = true;
      }
      else if (IsKeyDown(KEY_S))
      {
            anglePitch -= rotationMultiplier * deltaTime;
            if (anglePitch < 0)
            {
                  anglePitch += 360;
            }
            angleUpdated = true;

      }

      if (IsKeyDown(KEY_A))
      {
            angleYaw += rotationMultiplier * deltaTime;
            if (angleYaw > 360)
            {
                  angleYaw -= 360;
            }
            angleUpdated = true;
      }
      else if (IsKeyDown(KEY_D))
      {
            angleYaw -= rotationMultiplier * deltaTime;
            if (angleYaw < 0)
            {
                  angleYaw += 360;
            }
            angleUpdated = true;
      }

      if (IsKeyDown(KEY_Q))
      {
            angleRoll += rotationMultiplier * deltaTime;
            if (angleRoll > 360)
            {
                  angleRoll -= 360;
            }
            angleUpdated = true;
      }
      else if (IsKeyDown(KEY_E))
      {
            angleRoll -= rotationMultiplier * deltaTime;
            if (angleRoll < 0)
            {
                  angleRoll += 360;
            }
            angleUpdated = true;
      }

      if (angleUpdated) {
            airplane.transform = MatrixRotateXYZ((Vector3){DEG2RAD * anglePitch, DEG2RAD * angleYaw, DEG2RAD * angleRoll});
            getConsts(anglePitch, angleYaw, true, true);
      }
}