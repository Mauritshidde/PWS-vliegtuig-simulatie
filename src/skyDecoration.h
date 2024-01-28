#pragma once
#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "matrix.h"

class skyDecoration
{
private:
      Model birds;
      Model clouds;
      float anglePitch, angleYaw, angleRoll;
      std::vector<Vector3> allDecorations;
public:
      skyDecoration();
      ~skyDecoration();
      void draw();
};
