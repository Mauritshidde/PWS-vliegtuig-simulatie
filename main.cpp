#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>

// globally used variables like the airplane model
Vector2 previousMousePosition;
Model model;
Texture2D texture;
Vector2 cameraYZPos;
Vector3 cameraPos = { 0.0f, 0.0f, -120.0f };
Vector2 cameraXYPos = {cameraPos.x, cameraPos.y};
Camera camera = { 0 };

int turnDirection = 1;
int yMultiplier = 1;

 
float yCirclePos(float x, float radius) {
    float y;

    y = sqrt(pow(radius, 2) - pow(x, 2));

    return y;
}

void Start() {
    model = LoadModel("tinker.obj"); 
    texture = LoadTexture("Untitled1485_20230104061358.png"); 
}

void Update() {
    Vector2 currentMousePos = GetMousePosition();
    if (IsMouseButtonDown(0)) {
        // std::cout << "neew" << std::endl;
        if (previousMousePosition.x != currentMousePos.x) {
            // camera.target = (Vector3){ 0.0f, 100.0f, 0.0f };      // Camera looking at point
            float x = cameraPos.x;
            float xDisplacement = turnDirection * (currentMousePos.x-previousMousePosition.x);
            x += xDisplacement;
            if (x > 100) {
                x = 100 - (x-100);
                cameraPos.x = x;
                turnDirection = -turnDirection;
                yMultiplier = -yMultiplier;
            } else if (x < -100) {
                x = -100 - (x+100);
                cameraPos.x = x;
                turnDirection = -turnDirection;
                yMultiplier = -yMultiplier;
            } else {
                cameraPos.x = x + xDisplacement; 

            }
                cameraPos.y = yMultiplier * yCirclePos(cameraPos.x, 100);
                camera.position = cameraPos;

            std::cout << cameraPos.x << " " << yMultiplier*yCirclePos(x, 100) << " " << yMultiplier << " " << turnDirection << std::endl;
        }
    }
    previousMousePosition = currentMousePos;
}

void Render() {
    // render model
    BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);

            DrawModel(model, (Vector3){0.0f, 0.0f, 0.0f }, 0.4f, WHITE);
            DrawGrid(10, 10.0f);

        EndMode3D();
    
    EndDrawing();
}

int main() {
    // load model
    // cameraPos.y = -120.0f;

    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();
    InitWindow(screenWidth, screenHeight, "airplane simulation");

    camera.position = cameraPos;// Camera position perspective
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 30.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;

    SetTargetFPS(60);

    Start();

    while (!WindowShouldClose()) { 
        Update();
        Render();
        
    }

    

    return 0;
}