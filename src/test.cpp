#include <iostream>
#include <raylib.h>
#include <vector>

bool startScreen = true;

void Start() {

}
const char *text = "START";
void Render() {
    BeginDrawing();
        ClearBackground(WHITE);

        DrawText(text, 100, 100, 100, BLACK);
    EndDrawing();
}

void Update() {

}

int main() {
    InitWindow(0, 0, "airplane simulation");
    ToggleFullscreen();
    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();

    SetTargetFPS(60);
    Start();

    while (!WindowShouldClose())
    {
        float deltaTime = GetFrameTime();
        Update();
        Render();
    }


    return 0;
}