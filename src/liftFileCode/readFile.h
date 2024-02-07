#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>
#include <raylib.h>

class LiftFileReader
{
private:
    std::string fileName;
    nlohmann::json liftData;
    nlohmann::json liftWithYawData;
    nlohmann::json liftWithPitchData;

public:
    Vector2 getClCdWithYaw(float yawAngle);
    Vector2 getClCdWithPitch(float pitchAngle);
    Vector2 getClCdWithPitchAndYaw(float pitchAngle, float yawAngle);
    Vector2 getConstFromLiftFile(float pitchAngle, float yawAngle, bool withYaw, bool withPitch);

    LiftFileReader(std::string name = "Boeing737");
    ~LiftFileReader();
};