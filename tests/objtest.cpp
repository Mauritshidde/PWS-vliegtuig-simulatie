#include <raylib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#define RAYGUI_IMPLEMENTATION
#include "modules/raygui.h"
#include "gui_layout_name.h"
#include "gui/simulationGui.h"
#include <json/json.h>

// tijdelijke plaats voor variablen die bij een andere class horen
float maxAirspeed; // defined by mach number has to be lower than 1; speed is given in m/s
float engineTrust = 0.0f;
float maxEngineTrust = 1000.0f;
bool Button002Pressed = false;

std::vector<Vector3> v, vn;
std::vector<Vector2> vt;
std::vector<std::vector<std::string>> f;
std::vector<std::vector<Vector3>> data;

void read() {
    double a, b, c;
    std::string as, bs, cs;
    std::string v2;
    std::ifstream fin("tinker.txt");

    while(fin >> v2){
        if (v2 == "f"){
            fin >> as >> bs >> cs;
            f.push_back({as, bs, cs});
        } else if (v2 == "v") {
            fin >> a >> b >> c;
            v.push_back((Vector3){a, b, c});
        } else if (v2 == "vt") {
            fin >> a >> b;
            vt.push_back((Vector2){a, b});
        } else if (v2 == "vn") {
            fin >> a >> b >> c;
            vn.push_back((Vector3){a, b, c});
        }
        // std::cout << a << " " << b << " " <<  c << " " << std::endl;
    }
    
    for (int i=0; i < f.size(); i++) {
        std::vector<Vector3> data2;
        for (int k=0; k < 3; k++) {
            char z = 'e';
            int index = 0;
            int j = 0;
            std::vector<int> indices;
            // std::cout << f.at(i).at(k).size()-1 << " size" << std::endl;
            while (z != '/' && j < f.at(i).at(k).size()) {
                // std::cout << " fuck " << std::endl;
                z = f.at(i).at(k).at(j);
                if (z != '/') {
                    // std::cout << z << " z" << std::endl;
                    // std::cout << "ja " << z-'0' << std::endl;
                    indices.push_back(z - '0');
                }
                j++;

            }
            if (indices.size() >= 1) {
                int k2 = 0;
                // std::cout << "nee" << std::endl;
                for (int k=indices.size()-1; k >= 0; k--) {
                    index += indices.at(k) * pow(10, k2);
                    // std::cout << v.size() << " help me " << std::endl;
                    // std::cout << indices.at(k) * pow(10, k2) << " pow" << std::endl;
                    k2++;
                }
            }
            index--;
            // std::cout << index << std::endl;
            data2.push_back({v.at(index).x, v.at(index).y, v.at(index).z});
        }
        data.push_back(data2);
    }
}

void read2() {
    double a, b, c;
    std::string as, bs, cs;
    std::string v2;
    std::ifstream fin("tinker.txt");

    while(fin >> v2){
        if (v2 == "f"){
            fin >> as >> bs >> cs;
            f.push_back({as, bs, cs});
        } else if (v2 == "v") {
            fin >> a >> b >> c;
            v.push_back((Vector3){a, b, c});
        } else if (v2 == "vt") {
            fin >> a >> b;
            vt.push_back((Vector2){a, b});
        } else if (v2 == "vn") {
            fin >> a >> b >> c;
            vn.push_back((Vector3){a, b, c});
        }
    }
    
    for (int i=0; i < f.size(); i++) {
        std::vector<Vector3> data2;
        for (int k=0; k < 3; k++) {
            char z = 'e';
            int index = 0;
            int j = 0;
            std::vector<int> indices;
            while (z != '/' && j < f.at(i).at(k).size()) {
                z = f.at(i).at(k).at(j);
                if (z != '/') {
                    indices.push_back(z - '0');
                }
                j++;
            }
            if (indices.size() >= 1) {
                int k2 = 0;
                for (int k=indices.size()-1; k >= 0; k--) {
                    index += indices.at(k) * pow(10, k2);
                    k2++;
                }
            }
            index--;
            data2.push_back({v.at(index).x, v.at(index).y, v.at(index).z});
        }
        data.push_back(data2);
    }
}

class RunSimulation
{
private:
    Model airplane;
    Model skybox;
    Slider testtest;
    Button testtest2;

    Vector2 previousMousePosition;

    Texture2D airplaneTexture;
    Texture2D skyboxTexture;
    Vector2 cameraYZPos;
    Vector3 cameraPos;
    Vector2 cameraXYPos;
    Camera mainCamera;

    int renderWidth;
    int renderHeight;

    float angleYAxis = 0;
    float angleXZAxis = 0;
    float cameraCircleRadius = 120;

public:
    float test;
    RunSimulation();
    ~RunSimulation();

    bool notOnGUI(Vector2 mousePosition);
    void moveCamera(float deltaTime);
    void Start(int screenHeight, int screenWidth);
    void Update(float deltaTime);
    void Render();
    void run();
};

RunSimulation::RunSimulation()
{
}

RunSimulation::~RunSimulation()
{
}

float zCirclePosCam(float x, float radius)
{
    float z;

    z = sqrt(pow(radius, 2) - pow(x, 2));

    return z;
}

void RunSimulation::Start(int screenHeight, int screenWidth)
{
    read();
    renderWidth = GetRenderWidth();
    renderHeight = GetRenderHeight();
    test = 3;
    airplane = LoadModel("tinker.obj");
    airplaneTexture = LoadTexture("Untitled1485_20230104061358.png");
    // airplane.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = airplaneTexture;

    cameraPos = {0.0f, 0.0f, -120.0f};
    cameraXYPos = {cameraPos.x, cameraPos.y};
    mainCamera = {0};

    mainCamera.position = cameraPos;                  // Camera position perspective
    mainCamera.target = (Vector3){0.0f, 20.0f, 0.0f}; // Camera looking at point
    mainCamera.up = (Vector3){0.0f, 10.0f, 0.0f};     // Camera up vector (rotation towards target)
    mainCamera.fovy = 30.0f;                          // Camera field-of-view Y
    mainCamera.projection = CAMERA_PERSPECTIVE;

    skybox = LoadModel("skybox.obj");
    skyboxTexture = LoadTexture("skyboxtexture.png");
    skybox.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = skyboxTexture;
    testtest = Slider(0, maxEngineTrust, renderWidth - (renderWidth / 8) + (renderWidth / 38), renderHeight / 3.6, (renderWidth / 8) - (renderWidth / 38) * 2, renderHeight / 54);
    // testtest2 = Button(renderWidth - (renderWidth / 8) + (renderWidth / 38), renderHeight / 3.6, (renderWidth / 8) - (renderWidth / 38) * 2, renderHeight / 54);
}

bool RunSimulation::notOnGUI(Vector2 mousePosition)
{
    if (mousePosition.x >= renderWidth - (renderWidth / 8))
    {
        return false;
    }
    else
    {
        return true;
    }
}

void RunSimulation::moveCamera(float deltaTime)
{
    Vector2 currentMousePos = GetMousePosition();

    if (notOnGUI(currentMousePos))
    {
        if (IsMouseButtonDown(0))
        {
            angleYAxis += ((currentMousePos.x - previousMousePosition.x)) * deltaTime;
            angleXZAxis += ((currentMousePos.y - previousMousePosition.y)) * deltaTime;
        }

        if (IsKeyDown(KEY_RIGHT))
        {
            angleYAxis += deltaTime;
        }
        if (IsKeyDown(KEY_LEFT))
        {
            angleYAxis -= deltaTime;
        }
        if (IsKeyDown(KEY_UP))
        {
            angleXZAxis += deltaTime;
        }
        if (IsKeyDown(KEY_DOWN))
        {
            angleXZAxis -= deltaTime;
        }

        if (GetMouseWheelMove() > 0)
        {
            cameraCircleRadius += 1000 * deltaTime;
        }
        else if (GetMouseWheelMove() < 0)
        {
            cameraCircleRadius -= 1000 * deltaTime;
        }

        float x, y, z;
        x = cameraCircleRadius * sin(angleYAxis) * cos(angleXZAxis);
        y = cameraCircleRadius * sin(angleYAxis) * sin(angleXZAxis);
        z = cameraCircleRadius * cos(angleYAxis);

        mainCamera.position.x = x;
        mainCamera.position.y = y;
        mainCamera.position.z = z;

        previousMousePosition = currentMousePos;
    }
}

void RunSimulation::Update(float deltaTime)
{
    // first value updates over time
    // after that value updates by gui or key inputs
    moveCamera(deltaTime);
}
GuiLayoutNameState state = {0};

void RunSimulation::Render()
{
    Rectangle rec = {20, 40, 200, 150};
    Rectangle panelContentRec = {0, 0, 340, 340};
    Rectangle panelView = {0};
    Vector2 panelScroll = {99, -20};
    Rectangle sliderRec = {renderWidth - 240, 40, 200, 150};
    BeginDrawing();
    ClearBackground(BLACK);

    // DrawFPS(500, 500);

    BeginMode3D(mainCamera);
    if (data.size() >= 1) {
        for (int i=0; i < data.size(); i++) {
            // std::cout << data.at(i).at(0).y << std::endl;
            DrawTriangle3D(data.at(i).at(0), data.at(i).at(1), data.at(i).at(2), RED);
            DrawTriangle3D( {10, 10, 10}, {20, 20, 20}, {40, 30, -10}, RED);

            DrawLine3D(data.at(i).at(0), data.at(i).at(1), RED);
        }
    }
    // DrawModel(skybox, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f, skybox.materials->maps->color);
    // DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f}, (Vector3){180.0f, 0.0f, .0f}, 270.0f, (Vector3){0.4f, 0.4f, 0.4f}, GRAY);
    // DrawLine3D((Vector3){0.0f, 0.0f, 0.0f}, (Vector3){0.0f, 100.0f, 0.0f}, RED);
    DrawGrid(10, 10.0f);
    EndMode3D();

        BeginMode3D(mainCamera);
            
            // DrawModel(skybox, (Vector3){0.0f,0.0f,0.0f}, 1.0f, skybox.materials->maps->color);
            // DrawModelEx(airplane, (Vector3){0.0f, 0.0f, 0.0f }, (Vector3){180.0f, 0.0f, .0f }, 270.0f, (Vector3){0.4f,0.4f,0.4f}, GRAY);
            DrawLine3D((Vector3){0.0f, 0.0f, 0.0f }, (Vector3){0.0f, 100.0f, 0.0f }, RED);  
            DrawGrid(10, 10.0f);
        EndMode3D();
        GuiPanel((Rectangle){renderWidth - (renderWidth / 8), 0, (renderWidth / 8), renderHeight}, NULL);
        // GuiSlider((Rectangle){renderWidth - (renderWidth / 8) + (renderWidth/38), renderHeight/3.6, (renderWidth / 8) - (renderWidth/38)*2, renderHeight/54}, NULL, NULL, &engineTrust, 0, maxEngineTrust);
        testtest.DrawSlider();
        testtest2.DrawButton();
        // Button002Pressed = GuiButton((Rectangle){ 824, 288, 120, 24 }, "SAMPLE TEXT");
        // GuiLayoutName();
        // GuiGroupBox((Rectangle){ 66, 24, 276, 312 }, "STANDARD");
        // GuiSlider((Rectangle){ 96, 48, 216, 16 }, TextFormat("%0.f", value), NULL, &value, 0.0f, 1000.0f);
    EndDrawing();
}

void RunSimulation::run()
{   
    InitWindow(0, 0, "airplane simulation");
    ToggleFullscreen();
    const int screenWidth = GetScreenWidth();
    const int screenHeight = GetScreenHeight();

    SetTargetFPS(60);
    Start(screenWidth, screenHeight);

    while (!WindowShouldClose())
    {
        float deltaTime = GetFrameTime();
        Update(deltaTime);
        Render();
    }
}

int main()
{
    RunSimulation simulatie;
    simulatie.run();

    return 0;
}