#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>

// globally used variables like the airplane model
Vector2 previousMousePosition;
Model model;
Texture2D texture;
Vector2 cameraYZPos;
Vector3 cameraPos = {0.0f, 0.0f, -120.0f};
Vector2 cameraXYPos = {cameraPos.x, cameraPos.y};
Camera camera = {0};

float p = 0;
float l = 0;
float radius = 120;

float zCirclePosCam(float x, float radius)
{
    float z;

    z = sqrt(pow(radius, 2) - pow(x, 2));

    return z;
}

void Start()
{
    model = LoadModel("tinker.obj");
    // model = LoadModel("tinker.obj");
    texture = LoadTexture("Untitled1485_20230104061358.png");
}

void Update(float deltaTime)
{
    Vector2 currentMousePos = GetMousePosition();
    float x, y, z;

    if (IsMouseButtonDown(0))
    {
        p += ((currentMousePos.x - previousMousePosition.x)) * deltaTime;
        l += ((currentMousePos.y - previousMousePosition.y)) * deltaTime;
    }

    if (IsKeyDown(KEY_RIGHT))
    {
        p += deltaTime;
    }
    if (IsKeyDown(KEY_LEFT))
    {
        p -= deltaTime;
    }
    if (IsKeyDown(KEY_UP))
    {
        l += deltaTime;
    }
    if (IsKeyDown(KEY_DOWN))
    {
        l -= deltaTime;
    }

    if (GetMouseWheelMove() > 0)
    {
        radius += 100 * deltaTime;
    }
    else if (GetMouseWheelMove() < 0)
    {
        radius -= 100 * deltaTime;
    }

    x = radius * sin(p) * cos(l);
    y = radius * sin(p) * sin(l);
    z = radius * cos(p);

    std::cout << GetMouseWheelMove() << std::endl;
    
    camera.position.x = x;
    camera.position.y = y;
    camera.position.z = z;

    previousMousePosition = currentMousePos;
}

void Render()
{
    // render model
    BeginDrawing();
    ClearBackground(BLACK);

        BeginMode3D(camera);
            DrawModelEx(model, (Vector3){0.0f, 0.0f, 0.0f }, (Vector3){180.0f, 0.0f, .0f }, 270.0f, (Vector3){0.4f,0.4f,0.4f}, WHITE);
            DrawLine3D((Vector3){0.0f, 0.0f, 0.0f }, (Vector3){0.0f, 100.0f, 0.0f }, RED);  
            DrawGrid(10, 10.0f);
        EndMode3D();
    
    EndDrawing();
}

int main() {
    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();
    InitWindow(screenWidth, screenHeight, "airplane simulation");

    camera.position = cameraPos;                  // Camera position perspective
    camera.target = (Vector3){0.0f, 20.0f, 0.0f}; // Camera looking at point
    camera.up = (Vector3){0.0f, 10.0f, 0.0f};     // Camera up vector (rotation towards target)
    camera.fovy = 30.0f;                          // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;

    SetTargetFPS(60);

    Start();

    while (!WindowShouldClose())
    {
        float deltaTime = GetFrameTime();
        Update(deltaTime);
        Render();
    }

    return 0;
}