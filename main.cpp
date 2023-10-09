#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>

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

    float angleYAxis = 0;
    float angleXZAxis = 0;
    float cameraCircleRadius = 120;
public:
    RunSimulation(/* args */);
    ~RunSimulation();
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
    airplane = LoadModel("tinker.obj");
    airplaneTexture = LoadTexture("Untitled1485_20230104061358.png");
    airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;

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

void RunSimulation::Update(float deltaTime)
{
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

void RunSimulation::Render()
{
    // render model
    BeginDrawing();
    ClearBackground(BLACK);

        BeginMode3D(mainCamera);
            DrawModel(skybox, (Vector3){0.0f,0.0f,0.0f}, 1.0f, WHITE);
            DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f }, (Vector3){180.0f, 0.0f, .0f }, 270.0f, (Vector3){0.4f,0.4f,0.4f}, WHITE);
            DrawLine3D((Vector3){0.0f, 0.0f, 0.0f }, (Vector3){0.0f, 100.0f, 0.0f }, RED);  
            DrawGrid(10, 10.0f);
        EndMode3D();
    
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