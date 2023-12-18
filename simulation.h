#pragma once
#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#define RAYGUI_IMPLEMENTATION
#include "modules/raygui.h"

#define RAYMATH_IMPLEMENTATION
#include "modules/raymath.h"

#include "gui/simulationGui.h"
#include <json/json.h>

#include "Physics/ModelLoader.h"
#include "Physics/Plane.h"
#include "ui/menu.h"

// tijdelijke plaats voor variablen die bij een andere class horen

class RunSimulation
{
private:
    FluidDynamicsModel plane; // use this class instead of the model class for plane, because an error occurs when the model is loaded form the obj file
    Model airplane;
    Model skybox;

    Menu mainMenu;
    Plane planePhysicsModel;

    // Slider testtest;
    // Button testtest2;

    Vector2 previousMousePosition;

    Texture2D airplaneTexture;
    Texture2D skyboxTexture;
    Vector2 cameraYZPos;
    Vector3 cameraPos;
    Vector2 cameraXYPos;
    Camera mainCamera;

    int renderWidth;
    int renderHeight;

    float adiabaticIndex; // adiabatic index of air
    float gasConstant; // the gas constant of air
    float temperature; // the temperature
    float speedOfSound; // the speed of sound in the simulation

    float angleYAxis;
    float angleXZAxis;
    float cameraCircleRadius;

    float maxAngle;
    float currentPitchAngle;
    float currentYawAngle;
    float currentRollAngle;

    const char *minText = "0";
    const char *maxText = "360";

    bool notOnGUI(Vector2 mousePosition);
    void moveCamera(float deltaTime);
    void Start(int screenHeight, int screenWidth);
    void Update(float deltaTime);
    void Render();
public:
    RunSimulation();
    ~RunSimulation();

    void run();
};