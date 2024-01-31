#include "skyDecoration.h"

#define RAYMATH_IMPLEMENTATION
#include "../include/modules/raymath.h"

skyDecoration::skyDecoration()
{
      SetRandomSeed(GetTime());
      globalPosition = {0, 0, 0};
      birds = LoadModel("models/object/birds.obj");
      anglePitch = 0;
      angleYaw = 0;
      angleRoll = 0;
      birds.transform = MatrixRotateXYZ((Vector3){DEG2RAD * anglePitch, DEG2RAD * angleYaw, DEG2RAD * angleRoll});
     
      int amountOfDecorations = GetRandomValue(15, 20);
      Vector3 randomPos;
      for (int i = 0; i < amountOfDecorations; i++)
      {
            randomPos = {static_cast<float>(GetRandomValue(-300, 300)), static_cast<float>(GetRandomValue(-90, 90)), static_cast<float>(GetRandomValue(-300, 300))};
            allDecorations.push_back(randomPos);
      }
}

skyDecoration::~skyDecoration()
{
}

void skyDecoration::draw()
{
      if (isOutOfBounds()) globalPosition = Vector3Scale(Vector3Negate(Vector3Normalize(globalPosition)), skyboxRadius);
      
      for (int i = 0; i < allDecorations.size(); i++)
      {
            relativePosition.x = allDecorations.at(i).x + globalPosition.x;
            relativePosition.y = allDecorations.at(i).y + globalPosition.y;
            relativePosition.z = allDecorations.at(i).z + globalPosition.z;
            DrawModel(birds, relativePosition, 10.0f, BLACK);
      }      
}

bool skyDecoration::isOutOfBounds()
{
      if (globalPosition.x > skyboxRadius)
      {
            return true;
      }
      else if (globalPosition.x < -skyboxRadius)
      {
            return true;
      }
      if (globalPosition.y > skyboxRadius)
      {
            return true;
      }
      else if (globalPosition.y < -skyboxRadius)
      {
            return true;
      }
      if (globalPosition.z > skyboxRadius)
      {
            return true;
      }
      else if (globalPosition.z < -skyboxRadius)
      {
            return true;
      } 
      return false;
}
