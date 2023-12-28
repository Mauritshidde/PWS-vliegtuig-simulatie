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

void Plane::Start() {
      airplaneTexture = LoadTexture("models/texture/skyboxtexture.png");
      airplane = LoadModel("models/object/plane.obj");
      airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;
}

void Plane::Draw() {
      DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f}, (Vector3){1.0f, 0.0f, 0.0f}, 0, (Vector3){0.5f, 0.5f, 0.5f}, RED); // 2de vector geeft aan met welke factor hij met currentangle draait 
}

void Plane::Update(float deltaTime) {
      if (IsKeyDown(KEY_W)) {
            anglePitch += rotationMultiplier * deltaTime;
            if (anglePitch > 360) {
                  anglePitch -= 360;
            }
      } else if (IsKeyDown(KEY_S)) {
            anglePitch -= rotationMultiplier * deltaTime;
            if (anglePitch < 0) {
                  anglePitch += 360;
            }
      }

      if (IsKeyDown(KEY_A)) {
            angleYaw += rotationMultiplier * deltaTime;
            if (angleYaw > 360) {
                  angleYaw -= 360;
            }
      } else if (IsKeyDown(KEY_D)) {
            angleYaw -= rotationMultiplier * deltaTime;
            if (angleYaw < 0) {
                  angleYaw += 360;
            }
      }

      if (IsKeyDown(KEY_Q)) {
            angleRoll += rotationMultiplier * deltaTime;
            if (angleRoll > 360) {
                  angleRoll -= 360;
            }
      } else if (IsKeyDown(KEY_E)) {
            angleRoll -= rotationMultiplier * deltaTime;
            if (angleRoll < 0) {
                  angleRoll += 360;
            }
      }
      airplane.transform = MatrixRotateXYZ((Vector3){ DEG2RAD*anglePitch, DEG2RAD*angleYaw, DEG2RAD*angleRoll});
}