// class simulationGUI
// {
// private:
//     int screenWidth;
//     int screenHeight;

// public:
//     simulationGUI(int setScreenWidth, int setScreenHeight);
//     ~simulationGUI();

//     void drawGUI();
// };

// simulationGUI::simulationGUI(int setScreenWidth, int setScreenHeight)
// {
//     screenWidth = setScreenWidth;
//     screenHeight = setScreenHeight;
// }

// simulationGUI::~simulationGUI()
// {
// }

// void simulationGUI::drawGUI()
// {
//     GuiPanel((Rectangle){screenWidth - (screenWidth / 8), 0, screenWidth, screenHeight}, NULL);
// }
#pragma once
class Slider
{
private:
public:
    float sliderValue;
    float maxValue;
    float minValue;

    Rectangle sliderRectangle;

    Slider(float setMinValue = 0, float setMaxvalue = 100, int startX = 0, int startY = 0, int width = 10, int height = 10);
    ~Slider();

    void DrawSlider();
};

Slider::Slider(float setMinValue, float setMaxvalue, int startX, int startY, int width, int height)
{
    sliderValue = setMinValue;
    minValue = setMinValue;
    maxValue = setMaxvalue;

    sliderRectangle = (Rectangle){startX, startY, width, height};
}

Slider::~Slider()
{
}

void Slider::DrawSlider()
{
    GuiSlider(sliderRectangle, TextFormat("%0.f", minValue), TextFormat("%0.f", maxValue), &sliderValue, 0, maxValue);
}

class Button
{
private:
    Rectangle buttonRectangle;

public:
    Button(int startX = 0, int startY = 0, int width = 10, int height = 10);
    ~Button();

    void Update();
    void DrawButton();
    bool state;
};

Button::Button(int startX, int startY, int width, int height)
{
    buttonRectangle = (Rectangle){startX, startY, width, height};
    state = false;
}

Button::~Button()
{
}

void Button::DrawButton()
{
    if (GuiButton(buttonRectangle, "TEST"))
    {
        if (state)
        {
            state = false;
        }
        else
        {
            state = true;
        }
    }
}
