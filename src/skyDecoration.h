#pragma once
#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "matrix.h"
#include "Physics/Physics.h"

class skyDecoration
{
private:
      Model birds;
      Model clouds;
      float anglePitch, angleYaw, angleRoll;
      std::vector<Vector3> allDecorations;
      Vector3 relativePosition;
      int skyboxRadius = 400;
public:
      Physics physics = Physics();
      Vector3 globalPosition;
      skyDecoration() = default;
      skyDecoration(Vector3 position);
      ~skyDecoration();
      void draw();
      bool isOutOfBounds();
};
