#include "loadingScreen.h"
#include <iostream>
void loadingScreen(float fontSize, float screenWidth, float screenHeight) {
    const char *oneLetter = "A";
    const char *text = "Loading";

    float lenghtOfOneLetter = MeasureText(oneLetter, fontSize);
    
    float startX = (screenWidth / 2) - (MeasureText(text, fontSize) / 2);
    float startY = (screenHeight/2 + 50);
    float endX = ((screenWidth / 2) + (MeasureText(text, fontSize) / 2));
    float endY = (screenHeight/2 - 50);

    BeginDrawing();
        ClearBackground(WHITE);
        DrawText(text, (screenWidth / 2) - (MeasureText(text, fontSize) / 2), screenHeight/2 - 50, fontSize, RED);
        DrawLine(startX - (lenghtOfOneLetter / 2), startY + (lenghtOfOneLetter / 10), endX + lenghtOfOneLetter, startY + (lenghtOfOneLetter / 10), RED);
    EndDrawing();
}

void blinkingLoadingScreen(float fontSize, float screenWidth, float screenHeight, bool *door) {
    const char *oneLetter = "A";
    const char *text = "Loading";

    float lenghtOfOneLetter = MeasureText(oneLetter, fontSize);
    
    float startX = (screenWidth / 2) - (MeasureText(text, fontSize) / 2);
    float startY = (screenHeight/2 + 50);
    float endX = ((screenWidth / 2) + (MeasureText(text, fontSize) / 2));
    float endY = (screenHeight/2 - 50);
    int i=1;
    while (*door) {
        std::cout << "nee"  << *door << std::endl;
        i = i *-1;
        if (i == 1) {
            BeginDrawing();
                ClearBackground(WHITE);
                DrawText(text, (screenWidth / 2) - (MeasureText(text, fontSize) / 2), screenHeight/2 - 50, fontSize, RED);
                DrawLine(startX - (lenghtOfOneLetter / 2), startY + (lenghtOfOneLetter / 10), endX + lenghtOfOneLetter, startY + (lenghtOfOneLetter / 10), RED);
            EndDrawing();
        } else {
            BeginDrawing();
                ClearBackground(WHITE);
                DrawText(text, (screenWidth / 2) - (MeasureText(text, fontSize) / 2), screenHeight/2 - 50, fontSize, BLACK);
                DrawLine(startX - (lenghtOfOneLetter / 2), startY + (lenghtOfOneLetter / 10), endX + lenghtOfOneLetter, startY + (lenghtOfOneLetter / 10), RED);
            EndDrawing();
        }
    }
}