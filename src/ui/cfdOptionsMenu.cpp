#include "cfdOptionsMenu.h"

// #define RAYGUI_IMPLEMENTATION
// #include "../include/modules/raygui.h"

CfdMenu::CfdMenu(int screenWidth, int screenHeight)
{
      running = true;
      fontSize = 30;
      startScreen = true;
      selected = 0;
      isButtonPressed = false;
      buttonPressed = 0;

      lenghtOfOneLetter = MeasureText(oneLetter, fontSize);
      cfdMenuTexts = {"nx", "ny", "nz", "dx", "dy", "dz", "steps", "stepsizeYaw", "stepsizePitch"};
      cfdMenuTextsC = {"", "", "", "", "", "", "", "", ""};
      cfdMenuTextsType = {1, 1, 1, 1, 1, 1, 2, 2, 2}; // 2 = slider 1 = textbox
      thisV = 0;
      amount = cfdMenuTexts.size();
      val = screenHeight / (amount + 1);

      for (int i = 0; i < amount; i++)
      {
            startX.push_back((screenWidth / 2) - (MeasureText(cfdMenuTexts.at(i), fontSize) / 2));
            startY.push_back(val * (i + 1) + 50);
            endX.push_back((screenWidth / 2) + (MeasureText(cfdMenuTexts.at(i), fontSize) / 2));
            endY.push_back(val * (i + 1) - 50);
      }
}

CfdMenu::~CfdMenu()
{
}

void CfdMenu::Draw(int screenWidth, int screenHeight)
{
      BeginDrawing();
      ClearBackground(WHITE);

        for (int i = 0; i < cfdMenuTexts.size(); i++)
        {
                if (cfdMenuTextsType.at(i) == 1) {
                  //   DrawText(cfdMenuTexts.at(i), (screenWidth / 2) - (MeasureText(cfdMenuTexts.at(i), fontSize) / 2), val * (i + 1) - 50, fontSize, RED);
                  // Vector2 point = {(screenWidth / 2) - (MeasureText(cfdMenuTextsC.at(i), fontSize) / 2), val * (i + 1) - 50};
                  // Vector2 size = {100.0f,100.0f};
                  Rectangle box = {(screenWidth / 2) - (MeasureText(cfdMenuTextsC.at(i), fontSize) / 2), val * (i + 1) - 50, 100.0f, 100.0f};
                  // GuiTextBox(, cfdMenuTextsC.at(i), fontSize);
                  // GuiValueBox(box, cfdMenuTexts.at(i), &thisV, 1, 100, true);
                } else if (cfdMenuTextsType.at(i) == 2) {
                    DrawText(cfdMenuTexts.at(i), (screenWidth / 2) - (MeasureText(cfdMenuTexts.at(i), fontSize) / 2), val * (i + 1) - 50, fontSize, RED);
                  //   GuiSlider();
                }
                // if (selected == i)
                // {
                //     DrawText(cfdMenuTexts.at(i), (screenWidth / 2) - (MeasureText(cfdMenuTexts.at(i), fontSize) / 2), val * (i + 1) - 50, fontSize, RED);
                //     // DrawLine((screenWidth/2)-(MeasureText(cfdMenuTexts.at(i), 100)/2) - (lenghtOfOneLetter/2), val*(i+1)-50-(screenWidth/20), (screenWidth/2)-(MeasureText(cfdMenuTexts.at(i), 100)/2) + lenghtOfOneLetter, val*(i+1)-50-(screenWidth/20), RED);
                //     DrawLine(startX.at(i) - (lenghtOfOneLetter / 2), startY.at(i) + (lenghtOfOneLetter / 10), endX.at(i) + lenghtOfOneLetter, startY.at(i) + (lenghtOfOneLetter / 10), RED);
                // }
                // else
                // {
                //     DrawText(cfdMenuTexts.at(i), (screenWidth / 2) - (MeasureText(cfdMenuTexts.at(i), fontSize) / 2), val * (i + 1) - 50, fontSize, BLACK);
                // }
        }
      EndDrawing();
}

// void CfdMenu::Update(int screenWidth, int screenHeight)
// {
//       if (IsKeyPressed(KEY_UP))
//       {
//             selected--;
//             if (selected < 0)
//             {
//                   selected = 0;
//             }
//       }
//       else if (IsKeyPressed(KEY_DOWN))
//       {
//             selected++;
//             if (selected > cfdMenuTexts.size())
//             {
//                   selected = amount - 1;
//             }
//       }

//       if (1 == 2)
//       { // screenheight or screenwidth or fontSize changed
//             val = screenHeight / (amount + 1);
//             lenghtOfOneLetter = MeasureText(oneLetter, fontSize);

//             for (int i = 0; i < amount; i++)
//             {
//                   startX.push_back((screenWidth / 2) - (MeasureText(cfdMenuTexts.at(i), fontSize) / 2));
//                   startY.push_back(val * (i + 1) + 50);
//                   endX.push_back((screenWidth / 2) + (MeasureText(cfdMenuTexts.at(i), fontSize) / 2));
//                   endY.push_back(val * (i + 1) - 50);
//             }
//       }

//       if (IsMouseButtonPressed(0))
//       {
//             for (int i = 0; i < amount; i++)
//             {
//                   if (CheckCollisionPointRec(GetMousePosition(), {startX.at(i), endX.at(i), startY.at(i), endY.at(i)}))
//                   {
//                         buttonPressed = i;
//                         isButtonPressed = true;
//                         break;
//                   }
//             }
//       }
//       else if (IsKeyPressed(KEY_ENTER))
//       {
//             buttonPressed = selected;
//             isButtonPressed = true;
//       }

//       if (isButtonPressed)
//       {
//             startScreen = false;
//       }
// }