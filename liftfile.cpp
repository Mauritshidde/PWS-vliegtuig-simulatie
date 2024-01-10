#include "liftfile.h"

void fileWithoutYaw(float stepSize) {
    nlohmann::json data;

    int index = 0;
    for (float i=0; i <= 360; i+= stepSize) {
        data["pitch"][std::to_string(index)]["cl"] = 0.445; // calc lift and cl here
        data["pitch"][std::to_string(index)]["cd"] = 0.017;
        index++;
    }

    data["stepSize"] = stepSize;

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

    data["stepSize"] = stepSize;

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

Vector2 getClCdWithYaw(float yawAngle, nlohmann::json constFile) {
    float cl, cd;

    float stepSize = constFile["stepSize"].get<float>();

    float translatedYaw = yawAngle/stepSize;
    int yawIndex = (int) translatedYaw;

    int indexAdition = 1;
    if (yawAngle + stepSize > 360) {
        indexAdition = -yawIndex;
    }

    float Cl1 = constFile["yaw"][std::to_string(yawIndex)]["cl"].get<float>() * (1 - (translatedYaw - yawIndex));
    float Cl2 = constFile["yaw"][std::to_string(yawIndex + indexAdition)]["cl"].get<float>() * (translatedYaw - yawIndex);

    cl = Cl1 + Cl2;

    float Cd1 = constFile["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (1 - (translatedYaw - yawIndex));
    float Cd2 = constFile["yaw"][std::to_string(yawIndex + indexAdition)]["cd"].get<float>() * (translatedYaw - yawIndex);

    cd = Cd1 + Cd2;
    
    return {cl, cd};
}

Vector2 getClCdWithPitch(float pitchAngle, nlohmann::json constFile) {
    float cl, cd;

    float stepSize = constFile["stepSize"].get<float>();

    float translatedPitch = pitchAngle/stepSize;
    int pitchIndex = (int) translatedPitch;

    int indexAdition = 1;
    if (pitchIndex + stepSize > 360) {
        indexAdition = -pitchIndex;
    }

    float Cl1 = constFile["pitch"][std::to_string(pitchIndex)]["cl"].get<float>() * (1 - (translatedPitch - pitchIndex));
    float Cl2 = constFile["pitch"][std::to_string(pitchIndex + indexAdition)]["cl"].get<float>() * (translatedPitch - pitchIndex);

    cl = Cl1 + Cl2;

    float Cd1 = constFile["pitch"][std::to_string(pitchIndex)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex));
    float Cd2 = constFile["pitch"][std::to_string(pitchIndex + indexAdition)]["cd"].get<float>() * (translatedPitch - pitchIndex);

    cd = Cd1 + Cd2;
    
    return {cl, cd};
}

Vector2 getClCdWithPitchAndYaw(float pitchAngle, float yawAngle, nlohmann::json constFile) {
    float cl, cd;
    float stepSize = constFile["stepSize"].get<float>();
    float translatedPitch = pitchAngle/stepSize;
    int pitchIndex = (int) translatedPitch;

    float translatedYaw = yawAngle/stepSize;
    int yawIndex = (int) translatedYaw;

    int indexAditionYaw = 1;
    if (yawIndex + stepSize > 360) {
        indexAditionYaw = -yawIndex;
    }

    int indexAditionPitch = 1;
    if (pitchIndex + stepSize > 360) {
        indexAditionPitch = -pitchIndex;
    }

    float Cl1 = constFile["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex)]["cl"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (1 - (translatedYaw - yawIndex));
    float Cl2 = constFile["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex)]["cl"].get<float>() * (translatedPitch - pitchIndex) * (1 - (translatedYaw - yawIndex));

    float Cl3 = constFile["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cl"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (translatedYaw - yawIndex);
    float Cl4 = constFile["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cl"].get<float>() * (translatedPitch - pitchIndex) * (translatedYaw - yawIndex);

    cl = Cl1 + Cl2 + Cl3 + Cl4;

    float Cd1 = constFile["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (1 - (translatedYaw - yawIndex));
    float Cd2 = constFile["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (translatedPitch - pitchIndex) * (1 - (translatedYaw - yawIndex));

    float Cd3 = constFile["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (translatedYaw - yawIndex);
    float Cd4 = constFile["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cd"].get<float>() * (translatedPitch - pitchIndex) * (translatedYaw - yawIndex);

    cd = Cd1 + Cd2 + Cd3 + Cd4;
    
    return {cl, cd};
}

Vector2 getConstFromLiftFile(float pitchAngle, float yawAngle, bool withYaw, bool withPitch, std::string fileName) {
    float cl, cd;
    
    if (withYaw || withPitch) {
        std::ifstream f("planes/liftfiles/"+fileName+".json");
        nlohmann::json constFile = nlohmann::json::parse(f);
        f.close();
        Vector2 vals = getClCdWithPitchAndYaw(pitchAngle, yawAngle, constFile);
        cl = vals.x;
        cd = vals.y;

        return {cl, cd};
    } else if (withYaw) {
        std::ifstream f("planes/liftfiles/"+fileName+"Yaw.json");
        nlohmann::json constFile = nlohmann::json::parse(f);
        f.close();
        Vector2 vals = getClCdWithYaw(yawAngle, constFile);
        cl = vals.x;
        cd = vals.y;

        return {cl, cd};
    } else {
        std::ifstream f("planes/liftfiles/"+fileName+"Pitch.json");
        nlohmann::json constFile = nlohmann::json::parse(f);
        f.close();
        Vector2 vals = getClCdWithPitch(pitchAngle, constFile);
        cl = vals.x;
        cd = vals.y;

        return {cl, cd};
    }

    return {0, 0};
}