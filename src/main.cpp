#include "simulation.h"
#include "ui/menu.h"
#include "cfd.h"
#include "liftFileCode/createFile.h"
#include "ui/loadingScreen.h"

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
                Cfd cfd;
                cfd.run(30, 5, 5);
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