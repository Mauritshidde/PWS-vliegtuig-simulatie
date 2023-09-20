#include <raylib.h>
#include <iostream>
#include <vectors>

// globally used variables like the airplane model

void Start() {
    
}

void Update() {

}

void Render() {
    // render model
    BeginDrawing();
    
    EndDrawing();
}

int main() {
    // load model

    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();

    InitWindow(screenWidth, screenHeight, "airplane simulation");

    

    return 0;
}