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
int turnDirection2 = 1;
int zMultiplier = 1;
int zMultiplier2 = 1;
float p = 0;
float l = 0;
float radius = 120;

 
float zCirclePos(float x, float radius) {
    float z;

    z = sqrt(pow(radius, 2) - pow(x, 2));

    return z;
}

void Start() {
    model = LoadModel("tinker.obj"); 
    texture = LoadTexture("Untitled1485_20230104061358.png"); 
}

void Update(float deltaTime) {
    Vector2 currentMousePos = GetMousePosition();
    if (IsMouseButtonDown(0)) {
        // std::cout << "neew" << std::endl;
        // if (previousMousePosition.x != currentMousePos.x) {
        //     // camera.target = (Vector3){ 0.0f, 100.0f, 0.0f };      // Camera looking at point
        //     float x = cameraPos.x;
        //     float xDisplacement = turnDirection * (currentMousePos.x-previousMousePosition.x);
        //     x += xDisplacement;
        //     if (x > 100) {
        //         x = 100 - (x-100);
        //         cameraPos.x = x;
        //         turnDirection = -turnDirection;
        //         zMultiplier = -zMultiplier;
        //     } else if (x < -100) {
        //         x = -100 - (x+100);
        //         cameraPos.x = x;
        //         turnDirection = -turnDirection;
        //         zMultiplier = -zMultiplier;
        //     } else {
        //         cameraPos.x = x + xDisplacement; 

        //     }
        //         cameraPos.z = zMultiplier * zCirclePos(cameraPos.x, 100);
        //         camera.position = cameraPos;

        //     std::cout << cameraPos.x << " " << zMultiplier*zCirclePos(x, 100) << " " << zMultiplier << " " << turnDirection << std::endl;
        // }

        // angle between cameraPos and z axis
        // float vecLenght = sqrt(pow(cameraPos.x, 2) + pow(cameraPos.y, 2) + pow(cameraPos.z, 2));
        // float zLenght = pow(vecLenght.z, 2);
        // float vecvec = cameraPos.x * 0 +
        // pow(cos(), -1)
        
        
        // if (previousMousePosition.y != currentMousePos.y) {
        //     // camera.target = (Vector3){ 0.0f, 100.0f, 0.0f };      // Camera looking at point
        //     float y = cameraPos.y;
        //     float yDisplacement = turnDirection2 * (currentMousePos.y-previousMousePosition.y);
        //     y += yDisplacement;
        //     if (y > 100) {
        //         y = 100 - (y-100);
        //         cameraPos.y = y;
        //         turnDirection2 = -turnDirection2;
        //         zMultiplier2 = -zMultiplier2;
        //     } else if (y < -100) {
        //         y = -100 - (y+100);
        //         cameraPos.y = y;
        //         turnDirection2 = -turnDirection2;
        //         zMultiplier2 = -zMultiplier2;
        //     } else {
        //         cameraPos.y = y + yDisplacement; 

        //     }
        //         cameraPos.z = zMultiplier2 * zCirclePos(cameraPos.y, 100);
        //         camera.position = cameraPos;

        //     std::cout << cameraPos.y << " " << zMultiplier2*zCirclePos(y, 100) << " " << zMultiplier2 << " " << turnDirection2 << std::endl;
        // }
    }
    
    float x, y, z;

    if (IsMouseButtonDown(0)) {
        p += ((currentMousePos.x-previousMousePosition.x)) * deltaTime;
        l += ((currentMousePos.y-previousMousePosition.y))* deltaTime;
        
    }
    if (GetMouseWheelMove() > 0) {
        radius += 100 * deltaTime;
    } else if (GetMouseWheelMove() < 0) {
        radius -= 100 * deltaTime;
    }


    x = radius *  sin(p) * cos(l);
    y = radius * sin(p) * sin(l);
    z = radius * cos(p);
    std::cout << GetMouseWheelMove() << std::endl;
    camera.position.x = x;
    camera.position.y = y;
    camera.position.z = z;


    previousMousePosition = currentMousePos;
}

void Render() {
    // render model
    BeginDrawing();
        ClearBackground(BLACK);

        BeginMode3D(camera);

            // DrawModel(model, (Vector3){0.0f, 0.0f, 0.0f }, 0.4f, WHITE);
            DrawModelEx(model, (Vector3){0.0f, 0.0f, 0.0f }, (Vector3){180.0f, 0.0f, .0f }, 270.0f, (Vector3){0.4f,0.4f,0.4f}, WHITE);
            // DrawLine3D((Vector3){0.0f, 0.0f, 0.0f }, (Vector3){100.0f, 0.0f, 0.0f }, RED);  
            DrawLine3D((Vector3){0.0f, 0.0f, 0.0f }, (Vector3){0.0f, 100.0f, 0.0f }, RED);  
            // DrawLine3D((Vector3){0.0f, 0.0f, 0.0f }, (Vector3){0.0f, 0.0f, 100.0f }, RED);  
            DrawGrid(10, 10.0f);
    // void DrawModel(Model model, Vector3 position, float scale, Color tint);               // Draw a model (with texture if set)
    // void DrawModelEx(Model model, Vector3 position, Vector3 rotationAxis, float rotationAngle, Vector3 scale, Color tint); // Draw a model with extended parameters
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
    camera.target = (Vector3){ 0.0f, 20.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 10.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 30.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;

    SetTargetFPS(60);

    Start();

    while (!WindowShouldClose()) { 
        float deltaTime = GetFrameTime();
        Update(deltaTime);
        Render();
        
    }

    

    return 0;
}