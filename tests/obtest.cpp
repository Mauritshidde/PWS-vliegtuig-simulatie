#include "raylib.h"
#include "fstream"
#include <vector>
#include <iostream>
#include <string>
//------------------------------------------------------------------------------------
// Program main entry point
//------------------------------------------------------------------------------------
std::vector<Vector3> v, vn;
std::vector<Vector2> vt;
std::vector<std::vector<std::string>> f;
std::vector<std::vector<Vector3>> data;

void read() {
    double a, b, c;
    std::string as, bs, cs;
    std::string v2;
    std::ifstream fin("skybox.txt");

    while(fin >> v2){
        if (v2 == "f"){
            fin >> as >> bs >> cs;
            f.push_back({as, bs, cs});
        } else if (v2 == "v") {
            fin >> a >> b >> c;
            v.push_back((Vector3){a, b, c});
        } else if (v2 == "vt") {
            fin >> a >> b;
            vt.push_back((Vector2){a, b});
        } else if (v2 == "vn") {
            fin >> a >> b >> c;
            vn.push_back((Vector3){a, b, c});
        }
        // std::cout << a << " " << b << " " <<  c << " " << std::endl;
    }
    
    for (int i=0; i < f.size(); i++) {
        std::vector<Vector3> data2;
        for (int k=0; k < 3; k++) {
            char z;
            int index = 0;
            int j = 0;
            while (z != '/') {
                z = f.at(i).at(k).at(j);
                if (z != '/') {
                    index += z - '0';
                }
                j++;
            }
            data2.push_back({v.at(index).x, v.at(index).y, v.at(index).z});
        }
        data.push_back(data2);
    }
    
    // for (int i=0; i < f.size(); i++) {
        
    //     data2.push_back(v.at(int(f.at(i).at(0))));
    //     data2.push_back(v.at(int(f.at(i).at(1))));
    //     data2.push_back(v.at(int(f.at(i).at(2))));
    // }
}

int main(void)
{
    read();
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1920;
    const int screenHeight = 1080;

    InitWindow(screenWidth, screenHeight, "raylib [core] example - 3d camera free");

    Camera3D camera = { 0 };
    camera.position = (Vector3){ 10.0f, 10.0f, 10.0f }; // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type

    Vector3 cubePosition = { 0.0f, 0.0f, 0.0f };
    
    std::cout << data.size() << std::endl;
    SetTargetFPS(60);                   // Set our game to run at 60 frames-per-second

    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        UpdateCamera(&camera, CAMERA_FREE);

        if (IsKeyDown('Z')) camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };

        BeginDrawing();

            ClearBackground(BLACK);

            BeginMode3D(camera);
                if (data.size() >= 1) {
                    for (int i=0; i < data.size(); i++) {
                        DrawTriangle3D(data.at(i).at(0), data.at(i).at(1), data.at(i).at(2), RED);
                        // DrawTriangle3D( {10, 10, 10}, {20, 20, 20}, {40, 30, -10}, RED);

                        // DrawLine3D(data.at(i).at(0), data.at(i).at(1), RED);
                    }
                }
                DrawGrid(10, 1.0f);
                
            EndMode3D();

            // DrawRectangle( 10, 10, 320, 133, Fade(SKYBLUE, 0.5f));
            // DrawRectangleLines( 10, 10, 320, 133, BLUE);

            DrawText("Free camera default controls:", 20, 20, 10, BLACK);
            DrawText("- Mouse Wheel to Zoom in-out", 40, 40, 10, DARKGRAY);
            DrawText("- Mouse Wheel Pressed to Pan", 40, 60, 10, DARKGRAY);
            DrawText("- Alt + Mouse Wheel Pressed to Rotate", 40, 80, 10, DARKGRAY);
            DrawText("- Z to zoom to (0, 0, 0)", 40, 120, 10, DARKGRAY);

        EndDrawing();
    }

    CloseWindow();        // Close window and OpenGL context

    return 0;
}