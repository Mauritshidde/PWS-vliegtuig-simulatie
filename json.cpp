#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>

// using json = nlohmann::json;


void createLiftFile(int steps, float precisionFactor) { // steps has to be greater than 10 but realisticly has to be higher than 200
    nlohmann::json data;

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

    data["precisionFactor"] = precisionFactor; // the factor by which the values are multiplied to decide how much decimal points to keep

    std::ofstream liftfile;
    liftfile.open ("test.json");
    liftfile << data;
    liftfile.close();
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
    createLiftFile(steps, precisionFactor);
    
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
