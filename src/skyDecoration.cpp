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
     
      int amountOfDecorations = GetRandomValue(5, 12);
      Vector3 randomPos;
      for (int i = 0; i < amountOfDecorations; i++)
      {
            randomPos = {static_cast<float>(GetRandomValue(-200, 200)), static_cast<float>(GetRandomValue(-50, 50)), static_cast<float>(GetRandomValue(-200, 200))};
            allDecorations.push_back(randomPos);
      }
}

skyDecoration::~skyDecoration()
{
}

void skyDecoration::draw()
{
      Vector3 absolutePosition;
      for (int i = 0; i < allDecorations.size(); i++)
      {
            absolutePosition.x = allDecorations.at(i).x + globalPosition.x;
            absolutePosition.y = allDecorations.at(i).y + globalPosition.y;
            absolutePosition.z = allDecorations.at(i).z + globalPosition.z;
            DrawModel(birds, absolutePosition, 10.0f, BLACK);
      }      
}
