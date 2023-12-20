#pragma once 
#include "navierstokes.h"

int main(int argc, char const *argv[])
{
      NavierStokes fluidSimulation;
      fluidSimulation = NavierStokes();
      fluidSimulation.calc();
      return 0;
}