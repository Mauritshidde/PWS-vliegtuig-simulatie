#pragma once
#include <vector>
#include <raylib.h>

class CfdMenu
{
private:
    const char *oneLetter = "A";
    // const char *simuText = "nx";
    // const char *simuText = "ny";
    // const char *simuText = "nz";
    // const char *cdfText = "dx";
    // const char *addtext = "dy";
    // const char *settingsText = "dz";
    // const char *quitText = "Start";

    float fontSize;
    int amount;
    float val;

    bool isButtonPressed;
    int selected;

    float lenghtOfOneLetter;

    std::vector<const char *> cfdMenuTexts, cfdMenuTextsC;
    std::vector<int> cfdMenuTextsType;
    std::vector<float> startX, startY, endX, endY;

public:
    int buttonPressed;
    bool running;int thisV;
    bool startScreen;
    void Draw(int screenWidth, int screenHeight);
    void Update(int screenWidth, int screenHeight);

    CfdMenu(int screenWidth = 1920, int screenHeight = 1080);
    ~CfdMenu();
};
