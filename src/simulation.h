#pragma once
#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <thread>

// #include "../include/gui/simulationGui.h"
// #include <json/json.h>
#include "Physics/ModelLoader.h"
#include "Physics/Plane.h"
// #include "ui/loadingScreen.h"

class RunSimulation
{
private:
    Plane plane;
    // FluidDynamicsModel plane; // use this class instead of the model class for plane, because an error occurs when the model is loaded form the obj file
    // Model airplane;
    Model skybox;

    Plane planePhysicsModel;

    Vector2 previousMousePosition;

    // Texture2D airplaneTexture;
    Texture2D skyboxTexture;
    Vector2 cameraYZPos;
    Vector3 cameraPos;
    Vector2 cameraXYPos;
    Camera mainCamera;

    int renderWidth;
    int renderHeight;

    float adiabaticIndex; // adiabatic index of air
    float gasConstant;    // the gas constant of air
    float temperature;    // the temperature
    float speedOfSound;   // the speed of sound in the simulation

    float angleYAxis;
    float angleXZAxis;
    float cameraCircleRadius;
    float cameraRotationMultiplier;
    float cameraZoomMultiplier;

    float rotationMultiplier;
    float maxAngle;

    float airDensity; // ρ = P/(R×T)
    float height;

    const char *minText = "0";
    const char *maxText = "360";

    bool notOnGUI(Vector2 mousePosition);
    void moveCamera(float deltaTime);
    void loadSimulationVars();
    void Update(float deltaTime);
    void Render();
    void currentAirDensity();

public:
    RunSimulation();
    ~RunSimulation();

    void Start(int screenHeight, int screenWidth, bool *isLoading);
    bool loading = true;

    void run();
};