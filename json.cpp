#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>
#include <raylib.h>

// using json = nlohmann::json;


void createLiftFile(int steps, float precisionFactor, bool simpleStepType) { // steps has to be greater than 2 but realisticly has to be higher than 200
    nlohmann::json data;

    if (simpleStepType) {
        float stepSize = 1 * precisionFactor * 10;

        for (int i=0; i < 360; i+=stepSize) {
            for (int j=0; j < 360; j+=stepSize) {
                data["pitch"][std::to_string(i)]["yaw"][std::to_string(j)]["cl"] = 3;
                data["pitch"][std::to_string(i)]["yaw"][std::to_string(j)]["cd"] = 5;
            }
        }
        std::cout << "ja" << std::endl;
        data["simpleSteps"] = simpleStepType;
        data["precisionFactor"] = precisionFactor; // the factor by which the values are multiplied to decide how much decimal points to keep
        data["steps"] = 0; 
        data["spaceBetween"] = round(stepSize / precisionFactor);
        
    } else {
        float pitchAngle = 0;
        float yawAngle = 0;
        float stepSize = 360.0f/(steps-1);

        for (int i=0; i < steps; i++) {
            int roundedPitchAngle = round(pitchAngle / precisionFactor); // round to closest decimal point with an specific precision so it is as accurate as possible
            std::cout << roundedPitchAngle << " ";
            for (int j=0; j < steps; j++) {
                int roundedYawAngle = round(yawAngle / precisionFactor);
                data["pitch"][std::to_string(roundedPitchAngle)]["yaw"][std::to_string(roundedYawAngle)]["cl"] = 3;
                data["pitch"][std::to_string(roundedPitchAngle)]["yaw"][std::to_string(roundedYawAngle)]["cd"] = 5;
                yawAngle += stepSize;
            }
            pitchAngle += stepSize;
            yawAngle = 0;
        }
        std::cout << std::endl;

        data["simpleSteps"] = simpleStepType;
        data["precisionFactor"] = precisionFactor; // the factor by which the values are multiplied to decide how much decimal points to keep
        data["steps"] = steps; 
        data["spaceBetween"] = round(stepSize / precisionFactor);
    }

    std::ofstream liftfile;
    liftfile.open ("test.json");
    liftfile << data;
    liftfile.close();
}

void getConstFromLiftFile(float pitchAngle, float yawAngle) {
    float cl, cd;
    
    std::ifstream f("test.json");
    nlohmann::json data = nlohmann::json::parse(f);
    f.close();

    float precisionFactor = data["precisionFactor"];
    int stepSize = data["spaceBetween"];
    
    if (data["simpleSteps"]) {
        // pitchAngle 
    } else {
        int steps = data["steps"];
        std::cout << precisionFactor << " " << steps << std::endl;
        int roundedPitchAngle = round(pitchAngle / precisionFactor);


        // modules of pitchangle and yawangle by stepsize
        std::cout << stepSize << " " << (stepSize * 2) % stepSize << std::endl;
    }

    // float cl = 
    // return {cl, cd};
}

int main() {
    // std::ifstream f("2.json");
    // json data = json::parse(f);
    // json data;
    // json ex3 = {
    //     {"happy", true},
    //     {"pi", 3.141},
    // };

    float pitchAngle = 0;
    float yawAngle = 0;
    int steps = 100;

    float stepSize = 360.0f/(steps-1);
    // float x = 360.0f/(steps-1);
    float precisionFactor = 0.01;  // i.e. round to nearest one-hundreth
    createLiftFile(steps, precisionFactor, true);
    getConstFromLiftFile(3, 4);
    // float value = (int)(stepSize / precisionFactor) * precisionFactor;
    // std::cout << value << " " << x<< std::endl;

    // // auto j3 = json::parse(R"({"happy": true, "pi": 3.141})");
    // for (int i=0; i < steps; i++) {
    //     std::cout << stepSize * i << " ";
    //     for (int j=0; j < steps; j++) {
    //         data["pitch"][std::to_string(pitchAngle)]["yaw"][std::to_string(yawAngle)]["cl"] = 3;
    //         data["pitch"][std::to_string(pitchAngle)]["yaw"][std::to_string(yawAngle)]["cd"] = 5;
    //         yawAngle += stepSize;
    //     }
    //     pitchAngle += stepSize;
    //     yawAngle = 0;
    // }
    // std::cout << std::endl;
    // data["precisionFactor"] = precisionFactor; // amount of decimal places for pitch angle and yaw angle
    // // data["pitch"][std::to_string(1)]["yaw"][std::to_string(5)]["cl"] = 3;
    // // data["pitch"][std::to_string(1)]["yaw"][std::to_string(5)]["cd"] = 5;


    // std::ofstream testfile;
    // testfile.open ("test.json");
    // testfile << data;
    // testfile.close();


    return 0;
}
