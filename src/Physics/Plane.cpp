#include "Plane.h"

#define RAYMATH_IMPLEMENTATION
#include "../../include/modules/raymath.h"

//moments of inertia in kg meter^2
#define xMoment 24675886.91
#define yMoment 44877574.54
#define zMoment 67384152.71

Plane::Plane(std::string planeName, float startVelocity, float rho)
{
      std::ifstream f("planes/planeData.json");
      nlohmann::json planeData = nlohmann::json::parse(f);
      f.close();

      planePhysics = Physics();
      velocity = {0, 0, 0};
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;
      momentOfInertia = {xMoment, yMoment, zMoment};

      previousAnglePitch = anglePitch;
      previousAngleYaw = angleYaw;

      rotationMultiplier = 10; // multiplier for the speed of rotating the plane using WASD

      speed = startVelocity; // m/s
      speed = 257;
      planeFrontalArea = 10;

      wingArea = planeData["Planes"][planeName]["wing area"].get<float>(); // surface area of the wing in m2
      mass = planeData["Planes"][planeName]["maximal mass"].get<float>();

      liftFileName = planeName;

      files = LiftFileReader(liftFileName);
      Vector2 consts = getConsts(anglePitch, angleYaw, true, true);

      cl = consts.x;
      cd = consts.y;
      calcLift(rho); // set lift and drag
      physicsVector fG = physicsVector(planePhysics.calcForceGravity(mass), {0, 0, 0}); 
      forces.push_back(fG);
}

Plane::~Plane()
{
}

Vector2 Plane::getConsts(float pitch, float yaw, bool useYaw, bool usePitch)
{
      if (!usePitch && !useYaw)
      {
            usePitch = true;
            useYaw = true;
      }
      Vector2 consts = files.getConstFromLiftFile(pitch, yaw, useYaw, usePitch);
      return consts;
}

void Plane::calcLift(float rho)
{
      lift = cl * rho * pow(speed, 2) * wingArea * 0.5;
      drag = cd * rho * pow(speed, 2) * planeFrontalArea * 0.5;
}

Vector3 Plane::calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth)
{
      // TODO lift formula
      return {0, 0, 0};
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

      if (previousAnglePitch != anglePitch || previousAngleYaw != angleYaw || previousAngleRoll != angleRoll)
      {
            airplane.transform = MatrixRotateXYZ((Vector3){DEG2RAD * anglePitch, DEG2RAD * angleYaw, DEG2RAD * angleRoll});
            getConsts(anglePitch, angleYaw, false, true);
      }

      previousAnglePitch = anglePitch;
      previousAngleYaw = angleYaw;
      previousAngleRoll = angleRoll;

      calcLift(rho);
      
      evaluateForces(forces);
      updateVel(deltaTime);

      // std::cout << "speed: " << velocity << " lift: " << lift << " mass: " << 9.81 * mass << " Drag: " << drag << " pitch: " << anglePitch << " yaw: " << angleYaw << std::endl;
}

void Plane::evaluateForces(std::vector<physicsVector> forces)
{
     angularAcceleration = planePhysics.calcAngularAcceleration(forces, mass, centerOfMass, momentOfInertia);
     acceleration = planePhysics.calcAcceleration(forces, mass);
}

void Plane::updateVel(float deltaTime)
{
      Vector3 deltaVelocity = planePhysics.calcDeltaV(deltaTime, acceleration);
      velocity.x += deltaVelocity.x;
      velocity.y += deltaVelocity.y;
      velocity.z += deltaVelocity.z;
}
