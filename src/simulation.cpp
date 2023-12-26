#include "simulation.h"

#define RAYGUI_IMPLEMENTATION
#include "../include/modules/raygui.h"

#define RAYMATH_IMPLEMENTATION
#include "../include/modules/raymath.h"

RunSimulation::RunSimulation()
{
}

RunSimulation::~RunSimulation()
{
}

void RunSimulation::Start(int screenWidth, int screenHeight)
{
    plane.loadObjectModel();

    maxAngle = 360;
    currentPitchAngle = 0;
    currentYawAngle = 0;
    currentRollAngle = 0;

    renderWidth = GetRenderWidth();
    renderHeight = GetRenderHeight();

    // mainMenu = Menu(screenWidth, screenHeight);

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

    // testtest = Slider(0, planePhysicsModel.maxEngineTrust, renderWidth - (renderWidth / 8) + (renderWidth / 38), renderHeight / 3.6, (renderWidth / 8) - (renderWidth / 38) * 2, renderHeight / 54);
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
    // if (mainMenu.startScreen) {
    //     mainMenu.Update(GetScreenWidth(), GetScreenHeight());
    // } else {
        speedOfSound = sqrt(adiabaticIndex * gasConstant * temperature);

        if (planePhysicsModel.totalSpeed/speedOfSound <= 0.8) {
            planePhysicsModel.totalSpeed = 0.8 * speedOfSound;
        }
        // first value updates over time
        // after that value updates by gui or key inputs
        moveCamera(deltaTime);
    // }
    // if (IsMouseButtonPressed(0)) {
    // } else if (IsMouseButtonPressed(1)) {

    // }
}

void RunSimulation::Render()
{
    // if (mainMenu.startScreen) {
    //     mainMenu.Draw(GetScreenWidth(), GetScreenHeight());
    // } else {
        Rectangle sliderRec = {renderWidth - 240, 40, 200, 150};
        // Rectangle rec = {20, 40, 200, 150};
        // Rectangle panelContentRec = {0, 0, 340, 340};
        // Rectangle panelView = {0};
        // Vector2 panelScroll = {99, -20};
        BeginDrawing();
            ClearBackground(BLACK);

            BeginMode3D(mainCamera);
                DrawGrid(10, 10.0f);
            EndMode3D();

            BeginMode3D(mainCamera);
                DrawModel(skybox, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f, skybox.materials->maps->color);
                
                // airplane.transform = MatrixRotateXYZ((Vector3){ DEG2RAD*currentPitchAngle, DEG2RAD*currentYawAngle, DEG2RAD*currentRollAngle});
                DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f}, (Vector3){1.0f, 0.0f, 0.0f}, currentPitchAngle, (Vector3){0.5f, 0.5f, 0.5f}, RED); // 2de vector geeft aan met welke factor hij met currentangle draait 
                // DrawModel(airplane, (Vector3){0.0f, 0.0f, 0.0f}, 0.5f, BLACK);
                // plane.drawModel();
                DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){0.0f, 100.0f, 0.0f}, RED);
                DrawGrid(10, 10.0f);
            EndMode3D();
            
            Rectangle guiPanelSize = (Rectangle){renderWidth - (renderWidth / 8), 0, (renderWidth / 8), renderHeight};
            GuiPanel(guiPanelSize, NULL);
            // GuiSlider((Rectangle){renderWidth - (renderWidth / 8) + (renderWidth / 38), renderHeight / 3.6, (renderWidth / 8) - (renderWidth / 38) * 2, renderHeight / 54}, NULL, NULL, &planePhysicsModel.maxEngineTrust, 0, planePhysicsModel.currentEngineTrust);
            
            float guiTriangleWidth = guiPanelSize.width - guiPanelSize.x/25;
            float panelSliderWidthDifference = guiPanelSize.width - guiTriangleWidth; // differnce in width of the panel and the gui slider rectangle
            // float guiTriangleHeight = guiPanelSize; 

            GuiSlider((Rectangle){guiPanelSize.x + 0.5 * panelSliderWidthDifference, renderHeight / 3.6, guiTriangleWidth, renderHeight / 54}, minText, maxText, &currentPitchAngle, 0, maxAngle);
            GuiSlider((Rectangle){guiPanelSize.x + 0.5 * panelSliderWidthDifference, renderHeight / 5, guiTriangleWidth, renderHeight / 50}, minText, maxText, &currentYawAngle, 0, maxAngle);
            GuiSlider((Rectangle){guiPanelSize.x + 0.5 * panelSliderWidthDifference, renderHeight / 6.4, guiTriangleWidth, renderHeight / 46}, minText, maxText, &currentRollAngle, 0, maxAngle);
            // GuiSlider(Rectangle bounds, const char *textLeft, const char *textRight, float *value, float minValue, float maxValue);
            
            // testtest.DrawSlider();
            // testtest2.DrawButton();
            // Button002Pressed = GuiButton((Rectangle){ 824, 288, 120, 24 }, "SAMPLE TEXT");
            // GuiLayoutName();
            // GuiGroupBox((Rectangle){ 66, 24, 276, 312 }, "STANDARD");
            // GuiSlider((Rectangle){ 96, 48, 216, 16 }, TextFormat("%0.f", value), NULL, &value, 0.0f, 1000.0f);
        EndDrawing();
    // }
}

void RunSimulation::run()
{
    // InitWindow(0, 0, "airplane simulation");
    if (!IsWindowFullscreen()) {
        ToggleFullscreen();

    }
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