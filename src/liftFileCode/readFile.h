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
    Vector2 getClCdWithYaw(double yawAngle);
    Vector2 getClCdWithPitch(double pitchAngle);
    Vector2 getClCdWithPitchAndYaw(double pitchAngle, double yawAngle);
    Vector2 getConstFromLiftFile(double pitchAngle, double yawAngle, bool withYaw, bool withPitch);

    LiftFileReader(std::string name = "Boeing737");
    ~LiftFileReader();
};