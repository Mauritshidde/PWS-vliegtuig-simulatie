#pragma once
#include <vector>
#include <raylib.h>
class Menu
{
private:
    const char *oneLetter = "A";
    const char *simuText = "Simulation";
    const char *cdfText = "CFD";
    const char *addtext = "Add model";
    const char *settingsText = "Settings";
    const char *quitText = "Quit";

    float fontSize;
    int amount;
    float val;

    bool isButtonPressed;
    int selected;
    
    float lenghtOfOneLetter;

    std::vector<const char *> menuTexts;
    std::vector<float> startX, startY, endX, endY;

public:
    int buttonPressed;
    bool running;
    bool startScreen;
    void Draw(int screenWidth, int screenHeight);
    void Update(int screenWidth, int screenHeight);

    Menu(int screenWidth=1920, int screenHeight=1080);
    ~Menu();
};

