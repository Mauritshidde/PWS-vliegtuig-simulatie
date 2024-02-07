#include "createFile.h"

void fileWithoutYaw(std::vector<Vector2> *consts)
{
    nlohmann::json data;
    double stepSize = 360 / (consts->size() - 1);

    int index = 0;
    for (double i = 0; i <= 360; i += stepSize)
    {
        data["pitch"][std::to_string(index)]["cl"] = consts->at(index).x;
        data["pitch"][std::to_string(index)]["cd"] = consts->at(index).y;
        index++;
    }

    data["stepSize"] = stepSize;

    std::ofstream liftfile;
    liftfile.open("planes/liftfiles/Boeing737Pitch.json");
    liftfile << data;
    liftfile.close();
}

void fileWithoutPitch(std::vector<Vector2> *consts)
{
    nlohmann::json data;
    double stepSize = 360 / (consts->size() - 1);

    int index = 0;
    for (double i = 0; i <= 360; i += stepSize)
    {
        data["yaw"][std::to_string(index)]["cl"] = consts->at(index).x;
        data["yaw"][std::to_string(index)]["cd"] = consts->at(index).y;
        index++;
    }

    data["stepSize"] = stepSize;

    std::ofstream liftfile;
    liftfile.open("planes/liftfiles/Boeing737Yaw.json");
    liftfile << data;
    liftfile.close();
}

void fileWithBoth(std::vector<std::vector<Vector2>> *consts)
{
    nlohmann::json data;
    int steps = consts->size();
    float stepSize = 360.0f / (steps - 1);

    for (int i = 0; i < steps; i++)
    {
        for (int j = 0; j < steps; j++)
        {
            data["pitch"][std::to_string(i)]["yaw"][std::to_string(j)]["cl"] = consts->at(i).at(j).x;
            data["pitch"][std::to_string(i)]["yaw"][std::to_string(j)]["cd"] = consts->at(i).at(j).y;
        }
    }

    data["steps"] = steps;
    data["stepSize"] = stepSize;

    std::ofstream liftfile;
    liftfile.open("planes/liftfiles/Boeing737.json");
    liftfile << data;
    liftfile.close();
}

void createLiftFiles(std::vector<std::vector<Vector2>> *constsBoth, std::vector<Vector2> *constsPitch, std::vector<Vector2> *constsYaw)
{                                // steps has to be greater than 2 but realisticly has to be higher than 200
    fileWithoutYaw(constsPitch); // for the best results x times stepsize should be 360
    //std::cout << "done with Pitch file << " << std::endl;
    fileWithoutPitch(constsYaw);
    //std::cout << "done with Yaw file<< " << std::endl;
    fileWithBoth(constsBoth);
    //std::cout << "done with Yaw and Pitch file<< " << std::endl;
}