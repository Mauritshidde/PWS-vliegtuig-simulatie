#include "loadingScreen.h"

void loadingScreen(float fontSize, float screenWidth, float screenHeight, bool *repeat) {
    const char *oneLetter = "A";
    const char *text = "Loading";

    float lenghtOfOneLetter = MeasureText(oneLetter, fontSize);
    
    float startX = (screenWidth / 2) - (MeasureText(text, fontSize) / 2);
    float startY = (screenHeight/2 + 50);
    float endX = ((screenWidth / 2) + (MeasureText(text, fontSize) / 2));
    float endY = (screenHeight/2 - 50);

    int i = 1;
    while (*repeat) {
        int i = i * -1;
        
        switch (i)
        {
        case 1:
            BeginDrawing();
                ClearBackground(WHITE);
                DrawText(text, (screenWidth / 2) - (MeasureText(text, fontSize) / 2), screenHeight/2 - 50, fontSize, RED);
                DrawLine(startX - (lenghtOfOneLetter / 2), startY + (lenghtOfOneLetter / 10), endX + lenghtOfOneLetter, startY + (lenghtOfOneLetter / 10), RED);
            EndDrawing();
            break;
        case -1:
            BeginDrawing();
                ClearBackground(WHITE);
                DrawText(text, (screenWidth / 2) - (MeasureText(text, fontSize) / 2), screenHeight/2 - 50, fontSize, RED);
            EndDrawing();
            break;
        }

        
    }
}