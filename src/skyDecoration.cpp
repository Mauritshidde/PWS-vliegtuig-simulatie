#include "skyDecoration.h"

#define RAYMATH_IMPLEMENTATION
#include "../include/modules/raymath.h"

skyDecoration::skyDecoration()
{
      SetRandomSeed(GetTime());
      birds = LoadModel("models/object/birds.obj");
      anglePitch = 90;
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
      for (int i = 0; i < allDecorations.size(); i++)
      {
      DrawModel(birds, allDecorations.at(i), 10.0f, BLACK);
      }      
}
