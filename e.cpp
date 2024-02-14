#define WITHOUT_NUMPY
#include "src/Physics/matplotlibcpp.h"

namespace mat = matplotlibcpp;

std::vector<float> linspace(int startX, int endX, int steps)
{
    float stepSize = (endX - startX) / (steps - 1);
    std::vector<float> coords;
    for (int i = 0; i < steps; i++)
    {
        coords.push_back(startX + (stepSize * i));
    }
    return coords;
}

int main() {
    
    
    return 0;
}