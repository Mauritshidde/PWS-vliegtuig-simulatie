#include "Plane.h"

#define RAYMATH_IMPLEMENTATION
#include "../../include/modules/raymath.h"

Plane::Plane(float givenMass, Vector3 startingPos, float givenRotationMultiplier, float startVelocity, float rho)
{
      mass = givenMass;
      speedInDirections = {0, 0, 0};
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;

      previousAnglePitch = anglePitch;
      previousAngleYaw = angleYaw;
      
      rotationMultiplier = givenRotationMultiplier;

      velocity = startVelocity; // m/s
      wingArea = 10; // surface area of the wing in m2
      planeFrontalArea = 10;

      liftFileName = "Boeing737";
      files = LiftFileReader(liftFileName);
      Vector2 consts = getConsts(anglePitch, angleYaw, true, true);
      
      cl = consts.x;
      cd = consts.y;
      
      
      calcLift(rho); // set lift and drag
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
      Vector2 consts = files.getConstFromLiftFile(pitch, yaw, useYaw, usePitch);
      return consts;
}

void Plane::calcLift(float rho)
{
      lift = cl * pow(velocity, 2) * wingArea * 0.5;
      drag = cd * rho * pow(velocity, 2) * planeFrontalArea * 0.5;
      // TODO lift formula
}

Vector3 Plane::calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth)
{
      // TODO lift formula
}

void Plane::Start()
{
      airplaneTexture = LoadTexture("models/texture/planeTextureBeter.png");
      airplane = LoadModel("models/object/airplane.obj");
      airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;
}

void Plane::Draw()
{
      DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f}, (Vector3){1.0f, 0.0f, 0.0f}, 0, (Vector3){0.5f, 0.5f, 0.5f}, WHITE); // 2de vector geeft aan met welke factor hij met currentangle draait

}

void Plane::Update(float deltaTime, float rho)
{
      if (IsKeyDown(KEY_W))
      {
            anglePitch += rotationMultiplier * deltaTime;
            if (anglePitch > 360)
            {
                  anglePitch -= 360;
            }
      }
      else if (IsKeyDown(KEY_S))
      {
            anglePitch -= rotationMultiplier * deltaTime;
            if (anglePitch < 0)
            {
                  anglePitch += 360;
            }
      }

      if (IsKeyDown(KEY_A))
      {
            angleYaw += rotationMultiplier * deltaTime;
            if (angleYaw > 360)
            {
                  angleYaw -= 360;
            }
      }
      else if (IsKeyDown(KEY_D))
      {
            angleYaw -= rotationMultiplier * deltaTime;
            if (angleYaw < 0)
            {
                  angleYaw += 360;
            }
      }

      if (IsKeyDown(KEY_Q))
      {
            angleRoll += rotationMultiplier * deltaTime;
            if (angleRoll > 360)
            {
                  angleRoll -= 360;
            }
      }
      else if (IsKeyDown(KEY_E))
      {
            angleRoll -= rotationMultiplier * deltaTime;
            if (angleRoll < 0)
            {
                  angleRoll += 360;
            }
      }
      
      if (previousAnglePitch != anglePitch || previousAngleYaw != angleYaw || previousAngleRoll != angleRoll) {
            // std::cout << angleYaw << " " << anglePitch << std::endl;
            airplane.transform = MatrixRotateXYZ((Vector3){DEG2RAD * anglePitch, DEG2RAD * angleYaw, DEG2RAD * angleRoll});
            getConsts(anglePitch, angleYaw, true, true);
      }

      previousAnglePitch = anglePitch;
      previousAngleYaw = angleYaw;
      previousAngleRoll = angleRoll;

      calcLift(rho);
}