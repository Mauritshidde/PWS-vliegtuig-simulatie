#include "readFile.h"

Vector2 LiftFileReader::getClCdWithYaw(float yawAngle) {
    float cl, cd;

    float stepSize = liftWithYawData["stepSize"].get<float>();

    float translatedYaw = yawAngle/stepSize;
    int yawIndex = (int) translatedYaw;

    int indexAdition = 1;
    if (yawAngle + stepSize > 360) {
        indexAdition = -yawIndex;
    }

    float Cl1 = liftWithYawData["yaw"][std::to_string(yawIndex)]["cl"].get<float>() * (1 - (translatedYaw - yawIndex));
    float Cl2 = liftWithYawData["yaw"][std::to_string(yawIndex + indexAdition)]["cl"].get<float>() * (translatedYaw - yawIndex);

    cl = Cl1 + Cl2;

    float Cd1 = liftWithYawData["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (1 - (translatedYaw - yawIndex));
    float Cd2 = liftWithYawData["yaw"][std::to_string(yawIndex + indexAdition)]["cd"].get<float>() * (translatedYaw - yawIndex);

    cd = Cd1 + Cd2;
    
    return {cl, cd};
}

Vector2 LiftFileReader::getClCdWithPitch(float pitchAngle) {
    float cl, cd;

    float stepSize = liftWithPitchData["stepSize"].get<float>();

    float translatedPitch = pitchAngle/stepSize;
    int pitchIndex = (int) translatedPitch;

    int indexAdition = 1;
    if (pitchAngle + stepSize > 360) {
        indexAdition = -pitchIndex;
    }

    float Cl1 = liftWithPitchData["pitch"][std::to_string(pitchIndex)]["cl"].get<float>() * (1 - (translatedPitch - pitchIndex));
    float Cl2 = liftWithPitchData["pitch"][std::to_string(pitchIndex + indexAdition)]["cl"].get<float>() * (translatedPitch - pitchIndex);

    cl = Cl1 + Cl2;

    float Cd1 = liftWithPitchData["pitch"][std::to_string(pitchIndex)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex));
    float Cd2 = liftWithPitchData["pitch"][std::to_string(pitchIndex + indexAdition)]["cd"].get<float>() * (translatedPitch - pitchIndex);

    cd = Cd1 + Cd2;
    
    return {cl, cd};
}

Vector2 LiftFileReader::getClCdWithPitchAndYaw(float pitchAngle, float yawAngle) {
    float cl, cd;
    float stepSize = liftData["stepSize"].get<float>();
    float translatedPitch = pitchAngle/stepSize;
    int pitchIndex = (int) translatedPitch;

    float translatedYaw = yawAngle/stepSize;
    int yawIndex = (int) translatedYaw;

    int indexAditionYaw = 1;
    if (yawAngle + stepSize > 360) {
        indexAditionYaw = -yawIndex;
    }

    int indexAditionPitch = 1;
    if (pitchAngle + stepSize > 360) {
        indexAditionPitch = -pitchIndex;
    }

    float Cl1 = liftData["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex)]["cl"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (1 - (translatedYaw - yawIndex));
    float Cl2 = liftData["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex)]["cl"].get<float>() * (translatedPitch - pitchIndex) * (1 - (translatedYaw - yawIndex));

    float Cl3 = liftData["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cl"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (translatedYaw - yawIndex);
    float Cl4 = liftData["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cl"].get<float>() * (translatedPitch - pitchIndex) * (translatedYaw - yawIndex);

    cl = Cl1 + Cl2 + Cl3 + Cl4;

    float Cd1 = liftData["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (1 - (translatedYaw - yawIndex));
    float Cd2 = liftData["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex)]["cd"].get<float>() * (translatedPitch - pitchIndex) * (1 - (translatedYaw - yawIndex));

    float Cd3 = liftData["pitch"][std::to_string(pitchIndex)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cd"].get<float>() * (1 - (translatedPitch - pitchIndex)) * (translatedYaw - yawIndex);
    float Cd4 = liftData["pitch"][std::to_string(pitchIndex + indexAditionPitch)]["yaw"][std::to_string(yawIndex + indexAditionYaw)]["cd"].get<float>() * (translatedPitch - pitchIndex) * (translatedYaw - yawIndex);

    cd = Cd1 + Cd2 + Cd3 + Cd4;
    
    return {cl, cd};
}

Vector2 LiftFileReader::getConstFromLiftFile(float pitchAngle, float yawAngle, bool withYaw, bool withPitch) {
    if (withYaw || withPitch) {
        Vector2 vals = getClCdWithPitchAndYaw(pitchAngle, yawAngle);

        return vals;
    } else if (withYaw) {
        Vector2 vals = getClCdWithYaw(yawAngle);

        return vals;
    } else {
        Vector2 vals = getClCdWithPitch(pitchAngle);

        return vals;
    }

    return {0, 0}; // to make the compiler happy
}


LiftFileReader::LiftFileReader(std::string fileName)
{
    std::ifstream f("planes/liftfiles/"+fileName+".json");
    liftData = nlohmann::json::parse(f);
    f.close();
    
    std::ifstream g("planes/liftfiles/"+fileName+"Yaw.json");
    liftWithYawData = nlohmann::json::parse(g);
    g.close();

    std::ifstream h("planes/liftfiles/"+fileName+"Pitch.json");
    liftWithPitchData = nlohmann::json::parse(h);
    h.close();
}

LiftFileReader::~LiftFileReader()
{
}
