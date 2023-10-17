class simulationGUI
{
private:
    int screenWidth;
    int screenHeight;

public:
    simulationGUI(int setScreenWidth, int setScreenHeight);
    ~simulationGUI();

    void drawGUI();
};

simulationGUI::simulationGUI(int setScreenWidth, int setScreenHeight)
{
    screenWidth = setScreenWidth;
    screenHeight = setScreenHeight;
}

simulationGUI::~simulationGUI()
{
}

void simulationGUI::drawGUI()
{
    GuiPanel((Rectangle){screenWidth - (screenWidth / 8), 0, screenWidth, screenHeight}, NULL);
}

class Slider
{
private:
public:
    float sliderValue;
    float maxValue;
    float minValue;
    int startX, startY, width, height;

    Rectangle sliderRectangle;

    Slider(float setMinValue, float setMaxvalue, int setStartX, int setStartY, int setWidth, int setHeight);
    ~Slider();

    void DrawSlider();
};

Slider::Slider(float setMinValue = 0, float setMaxvalue = 100, int setStartX = 0, int setStartY = 0, int setWidth = 10, int setHeight= 10)
{
    sliderValue = setMinValue;
    minValue = setMinValue;
    maxValue = setMaxvalue;
    startX = setStartX;
    startY = setStartY;
    width = setWidth;
    height = setHeight;

    sliderRectangle = (Rectangle){startX, startY, width, height};
}

Slider::~Slider()
{
}

void Slider::DrawSlider()
{
    GuiSlider(sliderRectangle, TextFormat("%0.f", minValue), TextFormat("%0.f", maxValue), &sliderValue, 0, maxValue);
}