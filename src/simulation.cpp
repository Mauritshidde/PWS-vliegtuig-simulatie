#include "simulation.h"

#define RAYGUI_IMPLEMENTATION
#include "../include/modules/raygui.h"

#define RAYMATH_IMPLEMENTATION
#include "../include/modules/raymath.h"

namespace mat = matplotlibcpp;

RunSimulation::RunSimulation(std::string setFileName)
{
    fileName = setFileName;
}

RunSimulation::~RunSimulation()
{
}

void RunSimulation::Start(int screenWidth, int screenHeight)
{
    rotationMultiplier = 10;
    maxAngle = 360;
    rho = 1.225; // kg / m3

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
    skyboxTexture = LoadTexture("models/texture/skyboxtexture.png");
    skybox.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = skyboxTexture;

    plane = Plane(fileName, 100);
    
    plotXRange = linspace(0, 10, 101);
    timeElapsed = 0;
    // Vector2 aeroConsts;
    // for (float x = 0; x < plotXRange.size(); x++)
    // {
    //     aeroConsts = plane.getConsts(x*10, 0, false, true);
    //     plotYValues.push_back(aeroConsts.x);
    // }
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

        mainCamera.position = Vector3Transform(cameraPos, MatrixRotateXYZ((Vector3){DEG2RAD * angleXZAxis, DEG2RAD * angleYAxis, 0}));
        previousMousePosition = currentMousePos;
    }
}

void RunSimulation::Update(float deltaTime)
{
    speedOfSound = sqrt(adiabaticIndex * gasConstant * temperature);

    if (plane.speed / speedOfSound <= 0.8)
    {
        plane.speed = 0.8 * speedOfSound; //cap the speed of the plane at mach 0.8
    }
    // first value updates over time
    // after that value updates by gui or key inputs
    
    plane.Update(deltaTime, rho);

    moveCamera(deltaTime);
}

void RunSimulation::Render()
{
    Rectangle sliderRec = {renderWidth - 240, 40, 200, 150};

    BeginDrawing();
        ClearBackground(BLACK);
        DrawFPS(800,600);

        BeginMode3D(mainCamera);
            DrawModel(skybox, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f, skybox.materials->maps->color);
            plane.Draw();
            DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){0.0f, 100.0f, 0.0f}, RED);
            DrawGrid(10, 10.0f);
        EndMode3D();

        Rectangle guiPanelSize = (Rectangle){renderWidth - (renderWidth / 8), 0, (renderWidth / 8), renderHeight};
        GuiPanel(guiPanelSize, NULL);

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

    Start(screenWidth, screenHeight);

    while (!WindowShouldClose())
    {
        float deltaTime = GetFrameTime();
        Update(deltaTime);
        Render();
        std::cout << GetFPS() << std::endl;
        if (plotYValues.size() < plotXRange.size()) 
        {
            if (timeElapsed + deltaTime > 0.1)
            {
                timeElapsed = 0;
                std::cout << plotYValues.size() << " size \n";
                plotYValues.push_back(plane.angleYaw);
            }
            else
            {
                timeElapsed += deltaTime;
            }
        }
    }
    mat::plot(plotXRange, plotYValues, "-o");
    mat::save("test.pdf");
}