#include "simulation.h"
#include "ui/menu.h"
#include "navierstokes.h"
#include "liftFileCode/createFile.h"
#include "ui/loadingScreen.h"

#include <raylib.h>
#include <iostream>
#include <thread>
#include <chrono>

void loadingScreen2(float fontSize, float screenWidth, float screenHeight, bool *repeat) {
    const char *oneLetter = "A";
    const char *text = "Loading";

    float lenghtOfOneLetter = MeasureText(oneLetter, fontSize);
    
    float startX = (screenWidth / 2) - (MeasureText(text, fontSize) / 2);
    float startY = (screenHeight/2 + 50);
    float endX = ((screenWidth / 2) + (MeasureText(text, fontSize) / 2));
    float endY = (screenHeight/2 - 50);
            std::cout << *repeat << " kaas " << std::endl;
    
    int i = 1;
    while (*repeat) { // tried with a case but cases are stupid and it didn't work maybe my own mistake
        std::cout << "Ja" << std::endl;
        int i = i * -1;
        
        if (i == 1) {
            std::cout << "ja" << std::endl;
            BeginDrawing();
                ClearBackground(WHITE);
                DrawText(text, (screenWidth / 2) - (MeasureText(text, fontSize) / 2), screenHeight/2 - 50, fontSize, RED);
                DrawLine(startX - (lenghtOfOneLetter / 2), startY + (lenghtOfOneLetter / 10), endX + lenghtOfOneLetter, startY + (lenghtOfOneLetter / 10), RED);
            EndDrawing();
        } else {
            std::cout << "nee" << std::endl;
            BeginDrawing();
                ClearBackground(WHITE);
                DrawText(text, (screenWidth / 2) - (MeasureText(text, fontSize) / 2), screenHeight/2 - 50, fontSize, BLACK);
            EndDrawing();
            break;
        }
        // sleep_for(nanoseconds(10));
        // sleep(1);
    }
}

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
                std::thread t2(loadingScreen2, 100, screenWidth, screenHeight, &simulatie.loading);
                std::thread t(&RunSimulation::Start, &simulatie, screenHeight, screenWidth);
                // loadingScreen2(100, screenWidth, screenHeight, &simulatie.loading);
                t.join();
                t2.join();
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