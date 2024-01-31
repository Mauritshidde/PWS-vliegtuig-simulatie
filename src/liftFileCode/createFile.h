#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <iostream>
#include <raylib.h>

void fileWithoutYaw(std::vector<Vector2> *consts);
void fileWithoutPitch(std::vector<Vector2> *consts);
void fileWithBoth(std::vector<std::vector<Vector2>> *consts);
void createLiftFiles(std::vector<std::vector<Vector2>> *constsBoth, std::vector<Vector2> *constsPitch, std::vector<Vector2> *constsYaw);