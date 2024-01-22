#include "createFile.h"

void fileWithoutYaw(float stepSize) {
    nlohmann::json data;

    int index = 0;
    for (float i=0; i <= 360; i+= stepSize) {
        data["pitch"][std::to_string(index)]["cl"] = 0.445; // calc lift and cl here
        data["pitch"][std::to_string(index)]["cd"] = 0.017;
        index++;
    }

    int multiplier = 100;
    int stepSizeData = stepSize * multiplier;
    data["multiplier"] = multiplier;
    data["stepSize"] = stepSizeData;

    std::ofstream liftfile;
    liftfile.open ("planes/liftfiles/Boeing737Pitch.json");
    liftfile << data;
    liftfile.close();
}

void fileWithoutPitch(float stepSize) {
    nlohmann::json data;

    int index = 0;
    for (float i=0; i <= 360; i+= stepSize) {
        data["yaw"][std::to_string(index)]["cl"] = 0.445;
        data["yaw"][std::to_string(index)]["cd"] = 0.017;
        index++;
    }

    int multiplier = 100;
    int stepSizeData = stepSize * multiplier;
    data["multiplier"] = multiplier;
    data["stepSize"] = stepSizeData;

    std::ofstream liftfile;
    liftfile.open ("planes/liftfiles/Boeing737Yaw.json");
    liftfile << data;
    liftfile.close();
}

void fileWithBoth(int steps) {
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
    liftfile.open ("planes/liftfiles/Boeing737.json");
    liftfile << data;
    liftfile.close();
}


void createLiftFiles(int steps, float stepSize) { // steps has to be greater than 2 but realisticly has to be higher than 200
    fileWithoutYaw(stepSize); // for the best results x times stepsize should be 360
    fileWithoutPitch(stepSize);
    fileWithBoth(steps);
}