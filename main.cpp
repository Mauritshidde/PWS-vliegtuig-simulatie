#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#define RAYGUI_IMPLEMENTATION
#include "modules/raygui.h"
#include "gui/simulationGui.h"
#include <json/json.h>
#include "Physics/ModelLoader.h"
#include "Physics/Plane.h"
#include "ui/menu.h"

// tijdelijke plaats voor variablen die bij een andere class horen
float maxAirspeed = 1; // defined by mach number has to be lower than 1; this speed is given in mach whilst the speed in most other parts of the code is in m/s

bool Button002Pressed = false;

class RunSimulation
{
private:
    FluidDynamicsModel plane; // use this class instead of the model class for plane, because an error occurs when the model is loaded form the obj file
    Model airplane;
    Model skybox;

    Menu mainMenu;
    Plane planePhysicsModel;

    Slider testtest;
    Button testtest2;

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

RunSimulation::RunSimulation()
{
}

RunSimulation::~RunSimulation()
{
}

void RunSimulation::Start(int screenWidth, int screenHeight)
{
    plane.loadObjectModel();

    renderWidth = GetRenderWidth();
    renderHeight = GetRenderHeight();

    mainMenu = Menu(screenWidth, screenHeight);

    angleYAxis = 0;
    angleXZAxis = 0;
    cameraCircleRadius = 120;

    cameraPos = {0.0f, 0.0f, -120.0f};
    cameraXYPos = {cameraPos.x, cameraPos.y};
    mainCamera = {0};

    mainCamera.position = cameraPos;                  // Camera position perspective
    mainCamera.target = (Vector3){0.0f, 20.0f, 0.0f}; // Camera looking at point
    mainCamera.up = (Vector3){0.0f, 10.0f, 0.0f};     // Camera up vector (rotation towards target)
    mainCamera.fovy = 30.0f;                          // Camera field-of-view Y
    mainCamera.projection = CAMERA_PERSPECTIVE;

    skybox = LoadModel("models/object/skybox.obj");
    airplane = LoadModel("models/object/plane.obj");
    airplaneTexture = LoadTexture("models/texture/skyboxtexture.png");
    skyboxTexture = LoadTexture("models/texture/skyboxtexture.png");
    airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;
    skybox.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = skyboxTexture;
    
    planePhysicsModel = Plane(41145, {0,0,0});

    testtest = Slider(0, planePhysicsModel.maxEngineTrust, renderWidth - (renderWidth / 8) + (renderWidth / 38), renderHeight / 3.6, (renderWidth / 8) - (renderWidth / 38) * 2, renderHeight / 54);
}

bool RunSimulation::notOnGUI(Vector2 mousePosition)
{
    if (mousePosition.x >= renderWidth - (renderWidth / 8))
    {
        return false;
    }
    else
    {
        return true;
    }
}

void RunSimulation::moveCamera(float deltaTime)
{
    Vector2 currentMousePos = GetMousePosition();

    if (notOnGUI(currentMousePos))
    {
        if (IsMouseButtonDown(0))
        {
            angleYAxis += ((currentMousePos.x - previousMousePosition.x)) * deltaTime;
            angleXZAxis += ((currentMousePos.y - previousMousePosition.y)) * deltaTime;
        }

        if (IsKeyDown(KEY_RIGHT))
        {
            angleYAxis += deltaTime;
        }
        if (IsKeyDown(KEY_LEFT))
        {
            angleYAxis -= deltaTime;
        }
        if (IsKeyDown(KEY_UP))
        {
            angleXZAxis += deltaTime;
        }
        if (IsKeyDown(KEY_DOWN))
        {
            angleXZAxis -= deltaTime;
        }

        if (GetMouseWheelMove() > 0)
        {
            cameraCircleRadius += 100 * deltaTime;
        }
        else if (GetMouseWheelMove() < 0)
        {
            cameraCircleRadius -= 100 * deltaTime;
        }

        float x, y, z;
        x = cameraCircleRadius * sin(angleYAxis) * cos(angleXZAxis);
        y = cameraCircleRadius * sin(angleYAxis) * sin(angleXZAxis);
        z = cameraCircleRadius * cos(angleYAxis);

        mainCamera.position.x = x;
        mainCamera.position.y = y;
        mainCamera.position.z = z;

        previousMousePosition = currentMousePos;
    }
}

void RunSimulation::Update(float deltaTime)
{
    if (mainMenu.startScreen) {
        mainMenu.Update(GetScreenWidth(), GetScreenHeight());
    } else {
        speedOfSound = sqrt(adiabaticIndex * gasConstant * temperature);

        if (planePhysicsModel.totalSpeed/speedOfSound <= 0.8) {
            planePhysicsModel.totalSpeed = 0.8 * speedOfSound;
        }
        // first value updates over time
        // after that value updates by gui or key inputs
        moveCamera(deltaTime);
    }
    // if (IsMouseButtonPressed(0)) {
    // } else if (IsMouseButtonPressed(1)) {

    // }
}

void RunSimulation::Render()
{
    if (mainMenu.startScreen) {
        mainMenu.Draw(GetScreenWidth(), GetScreenHeight());
    } else {
        Rectangle rec = {20, 40, 200, 150};
        Rectangle panelContentRec = {0, 0, 340, 340};
        Rectangle panelView = {0};
        Vector2 panelScroll = {99, -20};
        Rectangle sliderRec = {renderWidth - 240, 40, 200, 150};
        BeginDrawing();
            ClearBackground(BLACK);

            BeginMode3D(mainCamera);
                DrawGrid(10, 10.0f);
            EndMode3D();

            BeginMode3D(mainCamera);
                DrawModel(skybox, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f, skybox.materials->maps->color);
                // DrawModel(airplane, (Vector3){0.0f, 0.0f, 0.0f}, 0.5f, BLACK);
                plane.drawModel();
                DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){0.0f, 100.0f, 0.0f}, RED);
                DrawGrid(10, 10.0f);
            EndMode3D();
            
            GuiPanel((Rectangle){renderWidth - (renderWidth / 8), 0, (renderWidth / 8), renderHeight}, NULL);
            GuiSlider((Rectangle){renderWidth - (renderWidth / 8) + (renderWidth / 38), renderHeight / 3.6, (renderWidth / 8) - (renderWidth / 38) * 2, renderHeight / 54}, NULL, NULL, &planePhysicsModel.maxEngineTrust, 0, planePhysicsModel.currentEngineTrust);

            testtest.DrawSlider();
            testtest2.DrawButton();
            // Button002Pressed = GuiButton((Rectangle){ 824, 288, 120, 24 }, "SAMPLE TEXT");
            // GuiLayoutName();
            // GuiGroupBox((Rectangle){ 66, 24, 276, 312 }, "STANDARD");
            // GuiSlider((Rectangle){ 96, 48, 216, 16 }, TextFormat("%0.f", value), NULL, &value, 0.0f, 1000.0f);
        EndDrawing();
    }
}

void RunSimulation::run()
{
    InitWindow(0, 0, "airplane simulation");
    ToggleFullscreen();
    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();

    SetTargetFPS(60);
    Start(screenWidth, screenHeight);

    while (!WindowShouldClose())
    {
        float deltaTime = GetFrameTime();
        Update(deltaTime);
        Render();
    }
}

int main()
{
    RunSimulation simulatie;
    simulatie.run();

    return 0;
}