#include "simulation.h"
#include "ui/menu.h"
#include "navierstokes.h"
#include "liftFileCode/createFile.h"

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
    // createLiftFiles(200, 0.01);
    while (!WindowShouldClose() && running)
    {
        float deltaTime = GetFrameTime();
        menu.Update(screenWidth, screenHeight);
        menu.Draw(screenWidth, screenHeight);

        if (!menu.startScreen)
        {
            if (menu.buttonPressed == 0)
            {
                RunSimulation simulatie;
                simulatie.run();
                running = false;
            }
            else if (menu.buttonPressed == 1)
            {
                NavierStokes cfd;
                cfd.calc();
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

        // RunSimulation simulatie;
        // simulatie.run();
    }

    return 0;
}