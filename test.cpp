#include <iostream>
#include <vector>
#include "liftfile.h"

void createLiftFile(int steps) { // steps has to be greater than 2 but realisticly has to be higher than 200
    nlohmann::json data;
    float stepSize = 360.0f/(steps-1);

    for (int i=0; i < steps; i++) {
        for (int j=0; j < steps; j++) {
            data["pitch"][std::to_string(i)]["yaw"][std::to_string(j)]["cl"] = 0.445;
            data["pitch"][std::to_string(i)]["yaw"][std::to_string(j)]["cd"] = 0.017;
        }
    }

    data["steps"] = steps; 
    data["stepSize"] = stepSize;

    std::ofstream liftfile;
    liftfile.open ("test.json");
    liftfile << data;
    liftfile.close();
}

Vector2 getClCdWithYaw(float pitchAngle, float yawAngle, nlohmann::json* constFile) {
    return {0,0};
}

Vector2 getClCdWithPitch(float pitchAngle, float yawAngle, nlohmann::json* constFile) {
    return {0,0};
}

Vector2 getClCdWithPitchAndYaw(float pitchAngle, float yawAngle, nlohmann::json* constFile) {
    float cl, cd;

    float stepSize = (*constFile)["stepSize"].get<float>();
    float translatedPitch = pitchAngle/stepSize;
    int pitchIndex = (int) translatedPitch;

    float translatedYaw = yawAngle/stepSize;
    int yawIndex = (int) translatedPitch;

    float Cl1 = (*constFile)["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex)]["cl"].get<double>() * (1 - (translatedPitch - pitchIndex)) * (1 - (translatedYaw - yawIndex));
    float Cl2 = (*constFile)["pitch"][std::to_string(pitchIndex + 1)]["yaw"][std::to_string(yawIndex)]["cl"].get<float>() * (translatedPitch - pitchIndex) * (1 - (translatedYaw - yawIndex));

    float Cl3 = (*constFile)["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex + 1)]["cl"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (translatedYaw - yawIndex);
    float Cl4 = (*constFile)["pitch"][std::to_string(pitchIndex + 1)]["yaw"][std::to_string(yawIndex + 1)]["cl"].get<float>() * (translatedPitch - pitchIndex) * (translatedYaw - yawIndex);

    cl = Cl1 + Cl2 + Cl3 + Cl4;

    float Cd1 = (*constFile)["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (1 - (translatedYaw - yawIndex));
    float Cd2 = (*constFile)["pitch"][std::to_string(pitchIndex + 1)]["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (translatedPitch - pitchIndex) * (1 - (translatedYaw - yawIndex));

    float Cd3 = (*constFile)["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex + 1)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (translatedYaw - yawIndex);
    float Cd4 = (*constFile)["pitch"][std::to_string(pitchIndex + 1)]["yaw"][std::to_string(yawIndex + 1)]["cd"].get<float>() * (translatedPitch - pitchIndex) * (translatedYaw - yawIndex);

    cd = Cd1 + Cd2 + Cd3 + Cd4;
    
    std::cout << translatedPitch << " " << cl << std::endl;

    return {cl, cd};
}

Vector2 getConstFromLiftFile(float pitchAngle, float yawAngle, bool withYaw, bool withPitch, nlohmann::json* constFile) {
    float cl, cd;
    
    if (withYaw || withPitch) {
        Vector2 vals = getClCdWithPitchAndYaw(pitchAngle, yawAngle, constFile);
        cl = vals.x;
        cd = vals.y;

        return {cl, cd};
    } else if (withYaw) {
        Vector2 vals = getClCdWithYaw(pitchAngle, yawAngle, constFile);
        cl = vals.x;
        cd = vals.y;

        return {cl, cd};
    } else {
        Vector2 vals = getClCdWithPitch(pitchAngle, yawAngle, constFile);
        cl = vals.x;
        cd = vals.y;

        return {cl, cd};
    }

}

int main() {
        std::ifstream f("test.json");
        nlohmann::json constFile = nlohmann::json::parse(f);
        f.close();
    float pitchAngle = 0;
    float yawAngle = 0;
    int steps = 100;

    float stepSize = 360.0f/(steps-1);
    float precisionFactor = 0.01;  // i.e. round to nearest one-hundreth
    // createLiftFile(steps, precisionFactor, false);
    Vector2 test = getConstFromLiftFile(3, 4, true, true, &constFile);
    return 0;
}

// int main() {
    // float x = 3.97;
    // float z = 2.4;

    // float r = x/z;

    // int h = (int) r;
    // std::cout << h << " " << r << std::endl;
    
    // return 0;
// }