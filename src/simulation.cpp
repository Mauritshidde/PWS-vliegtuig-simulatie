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
    // plane.loadObjectModel();
    rotationMultiplier = 10;
    maxAngle = 360;
    // currentPitchAngle = 0;
    // currentYawAngle = 0;
    // currentRollAngle = 0;

    renderWidth = GetRenderWidth();
    renderHeight = GetRenderHeight();

    angleYAxis = 0;
    angleXZAxis = 0;
    cameraCircleRadius = 120;
    cameraRotationMultiplier = 10;
    cameraZoomMultiplier = 30;
    
    cameraPos = {0.0f, 0.0f, cameraCircleRadius};
    cameraXYPos = {cameraPos.x, cameraPos.y};
    mainCamera = {0};

    mainCamera.position = cameraPos;                  // Camera position perspective
    mainCamera.target = (Vector3){0.0f, 0.0f, 0.0f}; // Camera looking at point  20 ?????????????? hier naar nog kijken TODO
    mainCamera.up = (Vector3){0.0f, 10.0f, 0.0f};     // Camera up vector (rotation towards target)
    mainCamera.fovy = 30.0f;                          // Camera field-of-view Y   effect van dit veranderen bestuderen ?????????????? TODO
    mainCamera.projection = CAMERA_PERSPECTIVE;  /// wat doet dit TODO

    skybox = LoadModel("models/object/skybox.obj");
    // airplane = LoadModel("models/object/plane.obj");
    // airplaneTexture = LoadTexture("models/texture/skyboxtexture.png");
    skyboxTexture = LoadTexture("models/texture/skyboxtexture.png");
    // airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;
    skybox.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = skyboxTexture;

    planePhysicsModel = Plane(41145, {0, 0, 0});

    plane.Start();
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
            angleYAxis += cameraRotationMultiplier * ((currentMousePos.x - previousMousePosition.x)) * deltaTime;
            angleXZAxis += cameraRotationMultiplier * ((currentMousePos.y - previousMousePosition.y)) * deltaTime;
            if (angleYAxis > 360) {
                angleYAxis -= 360;
            } else if (angleYAxis < 0) {
                angleYAxis += 360;
            }

            if (angleXZAxis > 360) {
                angleYAxis -= 360;
            } else if (angleYAxis < 0) {
                angleXZAxis += 360;
            }
        }

        if (IsKeyDown(KEY_RIGHT))
        {
            angleYAxis += cameraRotationMultiplier * deltaTime;
            if (angleYAxis > 360) {
                angleYAxis -= 360;
            }
        }
        if (IsKeyDown(KEY_LEFT))
        {
            angleYAxis -= cameraRotationMultiplier * deltaTime;
            if (angleYAxis < 0) {
                angleYAxis += 360;
            }
        }
        if (IsKeyDown(KEY_UP))
        {
            angleXZAxis += cameraRotationMultiplier * deltaTime;
            if (angleXZAxis > 360) {
                angleXZAxis -= 360;
            }
        }
        if (IsKeyDown(KEY_DOWN))
        {
            angleXZAxis -= cameraRotationMultiplier * deltaTime;
            if (angleXZAxis < 0) {
                angleXZAxis += 360;
            }
        }

        if (GetMouseWheelMove() > 0)
        {
            cameraCircleRadius += cameraZoomMultiplier * deltaTime;
            cameraPos = {0.0f, 0.0f, cameraCircleRadius};
        }
        else if (GetMouseWheelMove() < 0)
        {
            cameraCircleRadius -= cameraZoomMultiplier * deltaTime;
            cameraPos = {0.0f, 0.0f, cameraCircleRadius};
        }

        // float x, y, z;
        // x = cameraCircleRadius * sin(angleYAxis) * cos(angleXZAxis);
        // y = cameraCircleRadius * sin(angleYAxis) * sin(angleXZAxis);
        // z = cameraCircleRadius * cos(angleYAxis);

        // mainCamera.position.x = x;
        // mainCamera.position.y = y;
        // mainCamera.position.z = z;


        mainCamera.position = Vector3Transform(cameraPos, MatrixRotateXYZ((Vector3){DEG2RAD * angleXZAxis, DEG2RAD * angleYAxis, 0}));
        previousMousePosition = currentMousePos;
    }
}

void RunSimulation::Update(float deltaTime)
{
    speedOfSound = sqrt(adiabaticIndex * gasConstant * temperature);

    if (planePhysicsModel.totalSpeed / speedOfSound <= 0.8)
    {
        planePhysicsModel.totalSpeed = 0.8 * speedOfSound;
    }
    // first value updates over time
    // after that value updates by gui or key inputs

    plane.Update(deltaTime);

    moveCamera(deltaTime);
}

void RunSimulation::Render()
{
    Rectangle sliderRec = {renderWidth - 240, 40, 200, 150};

    BeginDrawing();
        ClearBackground(BLACK);

        BeginMode3D(mainCamera);
            DrawModel(skybox, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f, skybox.materials->maps->color);
            plane.Draw();
            // airplane.transform = MatrixRotateXYZ((Vector3){ DEG2RAD*plane.anglePitch, DEG2RAD*plane.angleYaw, DEG2RAD*plane.angleRoll});
            // DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f}, (Vector3){1.0f, 0.0f, 0.0f}, 0, (Vector3){0.5f, 0.5f, 0.5f}, RED); // 2de vector geeft aan met welke factor hij met currentangle draait
            // DrawModel(airplane, (Vector3){0.0f, 0.0f, 0.0f}, 0.5f, BLACK);
            // plane.drawModel();
            DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){0.0f, 100.0f, 0.0f}, RED);
            DrawGrid(10, 10.0f);
        EndMode3D();

        Rectangle guiPanelSize = (Rectangle){renderWidth - (renderWidth / 8), 0, (renderWidth / 8), renderHeight};
        GuiPanel(guiPanelSize, NULL);
        // GuiSlider((Rectangle){renderWidth - (renderWidth / 8) + (renderWidth / 38), renderHeight / 3.6, (renderWidth / 8) - (renderWidth / 38) * 2, renderHeight / 54}, NULL, NULL, &planePhysicsModel.maxEngineTrust, 0, planePhysicsModel.currentEngineTrust);

        float guiTriangleWidth = guiPanelSize.width - guiPanelSize.x / 25;
        float panelSliderWidthDifference = guiPanelSize.width - guiTriangleWidth; // differnce in width of the panel and the gui slider rectangle
        // float guiTriangleHeight = guiPanelSize;

        GuiSlider((Rectangle){guiPanelSize.x + 0.5 * panelSliderWidthDifference, renderHeight / 3.6, guiTriangleWidth, renderHeight / 54}, minText, maxText, &plane.anglePitch, 0, maxAngle);
        GuiSlider((Rectangle){guiPanelSize.x + 0.5 * panelSliderWidthDifference, renderHeight / 5, guiTriangleWidth, renderHeight / 50}, minText, maxText, &plane.angleYaw, 0, maxAngle);
        GuiSlider((Rectangle){guiPanelSize.x + 0.5 * panelSliderWidthDifference, renderHeight / 6.4, guiTriangleWidth, renderHeight / 46}, minText, maxText, &plane.angleRoll, 0, maxAngle);
        // GuiSlider(Rectangle bounds, const char *textLeft, const char *textRight, float *value, float minValue, float maxValue);

    EndDrawing();
}

void RunSimulation::run()
{
    // InitWindow(0, 0, "airplane simulation");
    if (!IsWindowFullscreen())
    {
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