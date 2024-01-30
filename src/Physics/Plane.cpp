#include "Plane.h"

#define RAYMATH_IMPLEMENTATION
#include "../../include/modules/raymath.h"

// moments of inertia in kg meter^2
#define xMoment 24675886.91
#define yMoment 44877574.54
#define zMoment 67384152.71

Plane::Plane(std::string planeName, float startVelocity, float rho)
{
      std::ifstream f("planes/planeData.json");
      nlohmann::json planeData = nlohmann::json::parse(f);
      f.close();

      airplaneTexture = LoadTexture("models/texture/planeTextureBeter.png");
      airplane = LoadModel("models/object/airplane.obj");
      airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;

      planePhysics = Physics();
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;
      momentOfInertia = {xMoment, yMoment, zMoment};

      previousAnglePitch = anglePitch;
      previousAngleYaw = angleYaw;
      previousAngleRoll = angleRoll;

      rotationMultiplier = 10; // multiplier for the speed of rotating the plane using WASD

      speed = startVelocity; // m/s
      speed = 257;
      
      // wingArea = planeData["Planes"][planeName]["wing area"].get<float>(); // surface area of the wing in m2
      mass = planeData["Planes"][planeName]["maximal mass"].get<float>();
      centerOfMass = {0, 0, 0};

      liftFileName = planeName;

      files = LiftFileReader(liftFileName);
      consts = getConsts(anglePitch, angleYaw, true, true);

      cl = consts.x;
      cd = consts.y;
      velocity = {0, 0, 0};
      acceleration = {0, 0, 0};
      angularVelocity = {0, 0, 0};
      angularAcceleration = {0, 0, 0};
      externalPos = {0, 0, 0};
      // std::cout << " xVel " << velocity.x << " yVel " << velocity.y << " zVel " << velocity.z << "\n";
      // std::cout << " xAccel " << acceleration.x << " yAccel " << acceleration.y << " zAccel " << acceleration.z << "\n";
      // std::cout << " xangvel " << angularAcceleration.x << " yangvel " << angularAcceleration.y << " zangvel " << angularAcceleration.z << "\n";
      // std::cout << " xangAccel " << angularAcceleration.x << " yangAccel " << angularAcceleration.y << " zangAccel " << angularAcceleration.z << "\n";

      engineOffset = 14;
      leftMotorThrust = {0,0, -maxEngineTrust};
      rightMotorThrust = {0,0, -maxEngineTrust};

      // physicsVector fG = physicsVector(planePhysics.calcForceGravity(mass), {centerOfMass.x + 10000, centerOfMass.y, centerOfMass.z});
      // leftMotorThrustDirection = Vector3Transform(leftMotorThrust, MatrixRotateXYZ((Vector3){DEG2RAD * angleYaw, DEG2RAD * anglePitch, DEG2RAD * angleRoll}));
      // Vector3 engineStartPosition = {centerOfMass.x - engineOffset, centerOfMass.y, centerOfMass.z};
      // Vector3 engineDirectionVec = planePhysics.vectorAddition(engineStartPosition, {0, 0, maxEngineTrust});
      // physicsVector forceRightMotor = physicsVector({0, 0, maxEngineTrust} , {centerOfMass.x + engineOffset, centerOfMass.y, centerOfMass.z});
      physicsVector forceLeftMotor = physicsVector(leftMotorThrust , {centerOfMass.x - engineOffset, centerOfMass.y, centerOfMass.z}); //TODO add variable engine thrust
      physicsVector forceRightMotor = physicsVector(rightMotorThrust , {centerOfMass.x + engineOffset, centerOfMass.y, centerOfMass.z}); //TODO add variable engine thrust
      
      // leftMotorDirectionPoint = planePhysics.vectorAddition(forceLeftMotor.location, forceLeftMotor.components);
      // rightMotorDirectionPoint = planePhysics.vectorAddition(forceRightMotor.location, {Vector3Normalize(forceRightMotor.components)});
      forces.push_back(forceLeftMotor);
      forces.push_back(forceRightMotor);
      // forces.push_back(forceRightMotor);
      // forces.push_back(fG);
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
      speed = planePhysics.calcHypot(velocity);
      if (speed)
      {
            std::cout << " xvel " << velocity.x << " yvel " << velocity.y << " zvel " << velocity.z << "\n";
            float angleYZ = atan(velocity.y / velocity.z) * 360; // pitch
            float angleXZ = atan(velocity.x / velocity.z) * 360; // yaw
            // float angleYX = atan(velocity.y / velocity.x) * 360; // roll
            std::cout << " YZ " << angleYZ << " xz " << angleXZ << " \n "; // << angleYX << "\n";
            Vector3 angleOfAttack = {angleYZ - anglePitch, angleXZ - angleYaw, 0}; //, angleYX - angleRoll};
            std::cout << " aoa.x " << angleOfAttack.x << " aoa.y " << angleOfAttack.y << " \n "; 
            angleOfAttack = reduceAngleDegrees(angleOfAttack);
            std::cout << " aoa.x " << angleOfAttack.x << " aoa.y " << angleOfAttack.y << " \n "; 
            consts = getConsts(angleOfAttack.x, angleOfAttack.y, (angleOfAttack.x != 0 ), (angleOfAttack.y != 0));
            cl = consts.x;
            cd = consts.y;
            std::cout << " cl" << cl << " cd " << cd << " \n";
      }
      lift = cl * rho * pow(speed, 2) * 0.5;
      drag = cd * rho * pow(speed, 2) * 0.5;
      dragDirection = Vector3Negate(Vector3Normalize(velocity));
      dragForce = Vector3Scale(dragDirection, drag);
      std::cout << " xdrag " << dragForce.x << " ydrag " << dragForce.y << " zdrag " << dragForce.z << "\n";
}

Vector3 Plane::calcCenterOfLiftWing(Vector3 startOfWing, Vector3 endOfWing, float startWingWidth, float endWingWidth)
{
      // TODO lift formula
      return {0, 0, 0};
}

void Plane::Draw()
{
      // 2de vector geeft aan met welke factor hij met currentangle draait
      std::cout << " ?? " << forces.at(0).components.x*100 << " " << forces.at(0).components.y*100 << " " << forces.at(0).components.z *100 << std::endl;
      DrawLine3D({0,0,0}, {forces.at(0).components.x/pow(10, 3), forces.at(0).components.y/pow(10, 3), forces.at(0).components.z /pow(10, 3)}, RED);
      DrawModelEx(airplane, externalPos, (Vector3){1.0f, 0.0f, 0.0f}, 0, (Vector3){0.5f, 0.5f, 0.5f}, WHITE); 
}

void Plane::Update(float deltaTime, float rho)
{
      if (IsKeyDown(KEY_W))
      {
            anglePitch += rotationMultiplier * deltaTime;
      }
      else if (IsKeyDown(KEY_S))
      {
            anglePitch -= rotationMultiplier * deltaTime;
      }

      if (IsKeyDown(KEY_A))
      {
            angleYaw += rotationMultiplier * deltaTime;
      }
      else if (IsKeyDown(KEY_D))
      {
            angleYaw -= rotationMultiplier * deltaTime;
      }

      if (IsKeyDown(KEY_Q))
      {
            angleRoll += rotationMultiplier * deltaTime;
      }
      else if (IsKeyDown(KEY_E))
      {
            angleRoll -= rotationMultiplier * deltaTime;
      }
      reduceAngleDegrees();
      if (previousAnglePitch != anglePitch || previousAngleYaw != angleYaw || previousAngleRoll != angleRoll)
      {
            airplane.transform = MatrixRotateXYZ((Vector3){DEG2RAD * anglePitch, DEG2RAD * angleYaw, DEG2RAD * angleRoll});
            rotatePoints();
      }

      previousAnglePitch = anglePitch;
      previousAngleYaw = angleYaw;
      previousAngleRoll = angleRoll;

      calcLift(rho);

      rotateVector();
      evaluateForces(forces);
      updateVel(deltaTime);
      updateAngularVel(deltaTime);
      updateRotation(deltaTime);
      externalPos = planePhysics.moveWithVelocity(externalPos, velocity, deltaTime); //TODO update pos in torque

      // physicsVector forceLeftMotor = physicsVector(leftMotorThrustDirection , {centerOfMass.x - engineOffset, centerOfMass.y, centerOfMass.z}); //TODO add variable engine thrust

      //test prints
      // for (int i = 0; i < forces.size(); i++)
      // {
      //       std::cout << " xf " << forces.at(i).components.x << " yf " << forces.at(i).components.y << " zf " << forces.at(i).components.z << "\n";
      //       std::cout << " xloc " << forces.at(i).location.x << " yloc " << forces.at(i).location.y << " zloc " << forces.at(i).location.z << "\n";
      // }
      // std::cout << "speed: " << velocity << " lift: " << lift << " mass: " << 9.81 * mass << " Drag: " << drag << " pitch: " << anglePitch << " yaw: " << angleYaw << std::endl;
}

void Plane::evaluateForces(std::vector<physicsVector> forces)
{
      angularAcceleration = planePhysics.calcAngularAcceleration(forces, mass, centerOfMass, momentOfInertia, (Vector3){anglePitch, angleYaw, angleRoll});
      acceleration = planePhysics.calcAcceleration(forces, mass);
      // std::cout << " xAccel " << acceleration.x << " yAccel " << acceleration.y << " zAccel " << acceleration.z << "\n";
}

void Plane::updateVel(float deltaTime)
{
      Vector3 deltaVelocity = planePhysics.calcDeltaV(deltaTime, acceleration);
      velocity.x += deltaVelocity.x;
      velocity.y += deltaVelocity.y;
      velocity.z += deltaVelocity.z;
      // std::cout << " xVel " << velocity.x << " yVel " << velocity.y << " zVel " << velocity.z << "\n";
}

void Plane::updateAngularVel(float deltaTime)
{
      Vector3 deltaAngularV = planePhysics.calcDeltaV(deltaTime, angularAcceleration);
      angularVelocity.x += deltaAngularV.x;
      angularVelocity.y += deltaAngularV.y;
      angularVelocity.z += deltaAngularV.z;
      // std::cout << " xDanV " << deltaAngularV.x << " yDanV " << deltaAngularV.y << " zDanV " << deltaAngularV.z << "\n";
      // std::cout << " xAn " << angularAcceleration.x << " yAn " << angularAcceleration.y << " zAn " << angularAcceleration.z << "\n";
      // std::cout << " xanvel " << angularVelocity.x << " yanvel " << angularVelocity.y << " zanvel " << angularVelocity.z << "\n";
}

void Plane::updateRotation(float deltaTime)
{
      anglePitch += angularVelocity.x * deltaTime * 360;
      angleYaw += angularVelocity.y * deltaTime * 360;
      angleRoll += angularVelocity.z * deltaTime * 360;
      reduceAngleDegrees();
}

void Plane::reduceAngleDegrees() // 0 < Angle < 360
{
      if (anglePitch > 360)
      {
            anglePitch -= 360;
      }
      else if (anglePitch < 0)
      {
            anglePitch += 360;
      }
      if (angleYaw > 360)
      {
            angleYaw -= 360;
      }
      else if (angleYaw < 0)
      {
            angleYaw += 360;
      }
      if (angleRoll > 360)
      {
            angleRoll -= 360;
      }
      else if (angleRoll < 0)
      {
            angleRoll += 360;
      }
}

Vector3 Plane::reduceAngleDegrees(Vector3 angle) // 0 < Angle < 360
{
      while (angle.x > 360)
      {
            angle.x -= 360;
      }
      while (angle.x < 0)
      {
            angle.x += 360;
      }
      while (angle.y > 360)
      {
            angle.y -= 360;
      }
      while (angle.y < 0)
      {
            angle.y += 360;
      }
      while (angle.z > 360)
      {
            angle.z -= 360;
      }
      while (angle.z < 0)
      {
            angle.z += 360;
      }
      return angle;
}

void Plane::rotatePoints()
{
      leftMotorDirectionPoint = Vector3Transform(leftMotorDirectionPoint, MatrixRotateXYZ((Vector3){DEG2RAD *anglePitch, DEG2RAD *angleYaw, DEG2RAD *angleRoll}));
}

void Plane::rotateVector()
{
      // leftMotorThrustDirection = Vector3Transform(leftMotorThrust, MatrixRotateXYZ((Vector3){DEG2RAD * angleYaw, DEG2RAD * anglePitch, DEG2RAD * angleRoll}));
      // rotatePoints();
      forces.at(0).components = Vector3Transform(leftMotorThrust, MatrixRotateXYZ((Vector3){DEG2RAD * anglePitch, DEG2RAD * angleYaw, DEG2RAD * angleRoll}));
      forces.at(1).components = Vector3Transform(rightMotorThrust, MatrixRotateXYZ((Vector3){DEG2RAD * anglePitch, DEG2RAD * angleYaw, DEG2RAD * angleRoll}));
      // forces.at(0).components = planePhysics.vectorSubtraction(leftMotorDirectionPoint, forces.at(0).location);
      // std::cout << sqrt(pow(forces.at(0).components.x, 2) + pow(forces.at(0).components.y, 2) + pow(forces.at(0).components.z, 2)) << " hhhhh" << std::endl;
      // std::cout << forces.at(0).components.x << " " << forces.at(0).components.y << " " << forces.at(0).components.z << " hhhhh" << std::endl;
}