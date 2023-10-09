#include <raylib.h>
#include <iostream>

Vector3 crossProduct(Vector3 vec1, Vector3 vec2)
{
    Vector3 result;
    result.x = vec1.y * vec2.z - vec1.z * vec2.y;
    result.y = vec1.z * vec2.x - vec1.x * vec2.z;
    result.z = vec1.x * vec2.y - vec1.y * vec2.x;
    return result;
}

int main()
{
    Vector3 locatie = {3, 2, 0};
    Vector3 Force = {10, 21, 1};
    Vector3 result = crossProduct(locatie, Force);

    std::cout << result.x << " " << result.y << " " << result.z << std::endl;

    return 0;
}