#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>
#include <raylib.h>

void createLiftFile(int steps); // steps has to be greater than 2 but realisticly has to be higher than 200 and lower than 1000 for perfomance reasons
void fileWithoutYaw(float precisionFactor, bool simpleStepType);
void fileWithoutPitch(float precisionFactor, bool simpleStepType);
Vector2 getClCdWithYaw(float yawAngle, nlohmann::json* constFile);
Vector2 getClCdWithPitch(float pitchAngle, nlohmann::json* constFile);
Vector2 getClCdWithPitchAndYaw(float pitchAngle, float yawAngle, nlohmann::json* constFile);
Vector2 getConstFromLiftFile(float pitchAngle, float yawAngle, nlohmann::json* constFile);