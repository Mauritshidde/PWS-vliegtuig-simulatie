#include <raylib.h>
#include <iostream>
#include <vector>

// globally used variables like the airplane model
Vector2 previousMousePosition;
Model model;
Texture2D texture;

Camera camera = { 0 };
 

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
            std::cout << previousMousePosition.x << " " << currentMousePos.x << std::endl;
        }
    }
    previousMousePosition = currentMousePos;
}

void Render() {
    // render model
    BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);

            DrawModel(model, (Vector3){ 0.0f, -8.0f, 0.0f }, 0.4f, WHITE);
            DrawGrid(10, 10.0f);

        EndMode3D();
    
    EndDrawing();
}

int main() {
    // load model

    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();
    InitWindow(screenWidth, screenHeight, "airplane simulation");

    camera.position = (Vector3){ 0.0f, 50.0f, -120.0f };// Camera position perspective
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