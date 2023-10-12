#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>

#define RAYGUI_IMPLEMENTATION
#include "modules/raygui.h"
#include "gui_layout_name.h"

class RunSimulation
{
private:
    Model airplane;
    Model skybox;

    Vector2 previousMousePosition;

    Texture2D airplaneTexture;
    Texture2D skyboxTexture;
    Vector2 cameraYZPos;
    Vector3 cameraPos;
    Vector2 cameraXYPos;
    Camera mainCamera;

    int renderWidth;
    int renderHeight;

    float angleYAxis = 0;
    float angleXZAxis = 0;
    float cameraCircleRadius = 120;
public:
    float test;
    RunSimulation(/* args */);
    ~RunSimulation();
    void moveCamera(float deltaTime);
    void Start();
    void Update(float deltaTime);
    void Render();
    void run();
};

RunSimulation::RunSimulation(/* args */)
{
}

RunSimulation::~RunSimulation()
{
}



float zCirclePosCam(float x, float radius)
{
    float z;

    z = sqrt(pow(radius, 2) - pow(x, 2));

    return z;
}

void RunSimulation::Start()
{
    // GuiLoadStyle("terminal.rgs");
    renderWidth = GetRenderWidth();
    renderHeight = GetRenderHeight();
    test = 3;
    airplane = LoadModel("tinker.obj");
    airplaneTexture = LoadTexture("Untitled1485_20230104061358.png");
    // airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;

    cameraPos = {0.0f, 0.0f, -120.0f};
    cameraXYPos = {cameraPos.x, cameraPos.y};
    mainCamera = {0};

    mainCamera.position = cameraPos;                  // Camera position perspective
    mainCamera.target = (Vector3){0.0f, 20.0f, 0.0f}; // Camera looking at point
    mainCamera.up = (Vector3){0.0f, 10.0f, 0.0f};     // Camera up vector (rotation towards target)
    mainCamera.fovy = 30.0f;                          // Camera field-of-view Y
    mainCamera.projection = CAMERA_PERSPECTIVE;

    skybox = LoadModel("skybox.obj");
    skyboxTexture = LoadTexture("skyboxtexture.png");
    skybox.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = skyboxTexture;
}

void RunSimulation::moveCamera(float deltaTime) {
    Vector2 currentMousePos = GetMousePosition();

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

void RunSimulation::Update(float deltaTime)
{
    moveCamera(deltaTime);
}

float value = 0.5f;
void RunSimulation::Render()
{
    Rectangle rec = { 20, 40, 200, 150 };
    Rectangle panelContentRec = {0, 0, 340, 340 };
    Rectangle panelView = { 0 };
    Vector2 panelScroll = { 99, -20 };
    Rectangle sliderRec = {renderWidth - 240, 40, 200, 150};
    BeginDrawing();
    ClearBackground(BLACK);

        BeginMode3D(mainCamera);
            DrawModel(skybox, (Vector3){0.0f,0.0f,0.0f}, 1.0f, skybox.materials->maps->color);
            DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f }, (Vector3){180.0f, 0.0f, .0f }, 270.0f, (Vector3){0.4f,0.4f,0.4f}, GRAY);
            DrawLine3D((Vector3){0.0f, 0.0f, 0.0f }, (Vector3){0.0f, 100.0f, 0.0f }, RED);  
            DrawGrid(10, 10.0f);
        EndMode3D();
        GuiLayoutName();
        GuiGroupBox((Rectangle){ 66, 24, 276, 312 }, "STANDARD");
        GuiSlider((Rectangle){ 96, 48, 216, 16 }, TextFormat("%0.f", value), NULL, &value, 0.0f, 1000.0f);
    
    EndDrawing();
}

void RunSimulation::run() 
{
    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();
    InitWindow(screenWidth, screenHeight, "airplane simulation");

    SetTargetFPS(60);

    Start();

    while (!WindowShouldClose())
    {
        float deltaTime = GetFrameTime();
        Update(deltaTime);
        Render();
    }
}

int main() {
    RunSimulation simulatie;
    simulatie.run();
    
    return 0;
}