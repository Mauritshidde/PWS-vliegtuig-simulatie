class Menu
{
private:
    const char *oneLetter = "A";
    const char *simuText = "Simulation";
    const char *cdfText = "CDf";
    const char *addtext = "Add module";
    const char *settingsText = "Settings";
    const char *quitText = "Quit";


    bool startScreen = true;
    int selected = 0;
    
    float lenghtOfOneLetter;
    std::vector<const char *> menuTexts;
    std::vector<float> startX, startY, endX, endY;

public:
    void Draw(int screenWidth, int screenHeight);
    void Update();

    Menu(/* args */);
    ~Menu();
};

Menu::Menu(/* args */)
{
    lenghtOfOneLetter = MeasureText(oneLetter, 100);
    menuTexts = {simuText, cdfText, addtext, settingsText, quitText};
}

Menu::~Menu()
{
}

void Menu::Draw(int screenWidth, int screenHeight) {
    BeginDrawing();
        ClearBackground(WHITE);

        int amount = menuTexts.size();
        float val = screenHeight/(amount+1);
        for (int i=0; i < menuTexts.size(); i++) {
            if (selected == i) {
                DrawText(menuTexts.at(i), (screenWidth/2)-(MeasureText(menuTexts.at(i), 100)/2), val*(i+1)-50, 100, RED);
                // DrawLine((screenWidth/2)-(MeasureText(menuTexts.at(i), 100)/2) - (lenghtOfOneLetter/2), val*(i+1)-50-(screenWidth/20), (screenWidth/2)-(MeasureText(menuTexts.at(i), 100)/2) + lenghtOfOneLetter, val*(i+1)-50-(screenWidth/20), RED);
            } else {
                DrawText(menuTexts.at(i), (screenWidth/2)-(MeasureText(menuTexts.at(i), 100)/2), val*(i+1)-50, 100, BLACK);
            }
        }
    EndDrawing();
}

void Menu::Update() {
    if (IsKeyPressed(KEY_UP)) {
        selected--;
        if (selected < 0) {
            selected = 0;
        }
    } else if (IsKeyPressed(KEY_DOWN)) {
        selected++;
        if (selected > menuTexts.size()) {
            selected = menuTexts.size()-1;
        }
    }
}