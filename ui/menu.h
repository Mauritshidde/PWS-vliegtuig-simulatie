class Menu
{
private:
    const char *oneLetter = "A";
    const char *simuText = "Simulation";
    const char *cdfText = "CDf";
    const char *addtext = "Add module";
    const char *settingsText = "Settings";
    const char *quitText = "Quit";

    float fontSize;
    int amount;
    float val;

    bool isButtonPressed;
    int buttonPressed;
    bool startScreen;
    int selected;
    
    float lenghtOfOneLetter;

    std::vector<const char *> menuTexts;
    std::vector<float> startX, startY, endX, endY;

public:
    void Draw(int screenWidth, int screenHeight);
    void Update(int screenWidth, int screenHeight);

    Menu(int screenWidth, int screenHeight);
    ~Menu();
};

Menu::Menu(int screenWidth, int screenHeight)
{
    fontSize = 100;
    startScreen = true;
    selected = 0;
    isButtonPressed = false;
    buttonPressed = 0;

    lenghtOfOneLetter = MeasureText(oneLetter, fontSize);
    menuTexts = {simuText, cdfText, addtext, settingsText, quitText};

    amount = menuTexts.size();
    val = screenHeight/(amount+1);
    
        
    for (int i=0; i < amount; i++) {
        startX.push_back((screenWidth/2)-(MeasureText(menuTexts.at(i), fontSize)/2));
        startY.push_back(val*(i+1)+50);
        endX.push_back((screenWidth/2)+(MeasureText(menuTexts.at(i), fontSize)/2));
        endY.push_back(val*(i+1)-50);
    }
}

Menu::~Menu()
{
}

void Menu::Draw(int screenWidth, int screenHeight) {
    BeginDrawing();
        ClearBackground(WHITE);

        for (int i=0; i < menuTexts.size(); i++) {
            if (selected == i) {
                DrawText(menuTexts.at(i), (screenWidth/2)-(MeasureText(menuTexts.at(i), fontSize)/2), val*(i+1)-50, fontSize, RED);
                // DrawLine((screenWidth/2)-(MeasureText(menuTexts.at(i), 100)/2) - (lenghtOfOneLetter/2), val*(i+1)-50-(screenWidth/20), (screenWidth/2)-(MeasureText(menuTexts.at(i), 100)/2) + lenghtOfOneLetter, val*(i+1)-50-(screenWidth/20), RED);
            } else {
                DrawText(menuTexts.at(i), (screenWidth/2)-(MeasureText(menuTexts.at(i), fontSize)/2), val*(i+1)-50, fontSize, BLACK);
            }
        }
    EndDrawing();
}

void Menu::Update(int screenWidth, int screenHeight) {
    if (IsKeyPressed(KEY_UP)) {
        selected--;
        if (selected < 0) {
            selected = 0;
        }
    } else if (IsKeyPressed(KEY_DOWN)) {
        selected++;
        if (selected > menuTexts.size()) {
            selected = amount-1;
        }
    }

    if (1 == 2) { //screenheight or screenwidth or fontSize changed
        val = screenHeight/(amount+1);
        lenghtOfOneLetter = MeasureText(oneLetter, fontSize);

        for (int i=0; i < amount; i++) {
            startX.push_back((screenWidth/2)-(MeasureText(menuTexts.at(i), fontSize)/2));
            startY.push_back(val*(i+1)+50);
            endX.push_back((screenWidth/2)+(MeasureText(menuTexts.at(i), fontSize)/2));
            endY.push_back(val*(i+1)-50);
        }   
    } 


    if (IsMouseButtonPressed(0)) {
        for (int i=0; i < amount; i++) {
            if (CheckCollisionPointRec(GetMousePosition(), {startX.at(i), endX.at(i), startY.at(i), endY.at(i)})) {
                buttonPressed = i;
                isButtonPressed = true;
                break;
            }
        }
    }

    if (isButtonPressed) {
        if (buttonPressed == 0) {
            // start simulation;
        } else if (buttonPressed == 1) {
            // start cdf program for calculating cl and cd;
        } else if (buttonPressed == 2) {
            // load application for loading a model // selecting parts from the model;
        } else if (buttonPressed == 3) {
            // open settings
        } else if (buttonPressed == 4) {
            // quit the application
        }
    }
}