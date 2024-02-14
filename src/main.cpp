#include "simulation.h"
#include "ui/menu.h"
#include "cfd.h"
#include "ui/loadingScreen.h"
// #include "ui/cfdOptionsMenu.h"

#include <raylib.h>
#include <iostream>

int main()
{
    bool running = true;
    Menu menu;

    InitWindow(0, 0, "airplane simulation");
    ToggleFullscreen();
    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();
    // SetConfigFlags(FLAG_MSAA_4X_HINT);
    SetConfigFlags(FLAG_FULLSCREEN_MODE | FLAG_WINDOW_RESIZABLE);
    // createLiftFiles(200, 0.01);
    SetTargetFPS(60);
    menu = Menu(screenWidth, screenHeight);
    while (!WindowShouldClose() && running)
    {
        float deltaTime = GetFrameTime();
        menu.Update(screenWidth, screenHeight);
        menu.Draw(screenWidth, screenHeight);

        if (!menu.startScreen)
        {
            if (menu.buttonPressed == 0)
            {
                loadingScreen(100, screenWidth, screenHeight);
                RunSimulation simulatie("Boeing737");
                simulatie.run();
                running = false;
            }
            else if (menu.buttonPressed == 1)
            {
                // CfdMenu cfdMenu = CfdMenu(screenWidth, screenHeight);
                // cfdMenu.Draw(screenWidth, screenHeight);
                Cfd cfd(90, 100, 100, 0.1, 300, 1.293, false, false);
                cfd.run(20, 1, 3);
                // cfd.Draw();
                running = false;
                // start cdf program for calculating cl and cd;
            }
            else if (menu.buttonPressed == 2)
            {
                // load application for loading a model // selecting parts from the model;
                running = false;
            }
            else if (menu.buttonPressed == 3)
            {
                // open settings
                running = false;
            }
            else if (menu.buttonPressed == 4)
            {
                running = false;
            }
        }
    }

    return 0;
}