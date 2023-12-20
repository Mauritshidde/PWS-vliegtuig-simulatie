#include "ModelLoader.h"

FluidDynamicsModel::FluidDynamicsModel()
{
}

FluidDynamicsModel::~FluidDynamicsModel()
{
}

void FluidDynamicsModel::readObjectModelFile(std::string *as, std::string *bs, std::string *cs, float *a, float *b, float *c)
{
    std::ifstream fin("models/object/txt/tinker2.txt");
    std::string v2;

    while (fin >> v2)
    {
        if (v2 == "f")
        {
            fin >> *as >> *bs >> *cs;
            f.push_back({*as, *bs, *cs});
        }
        else if (v2 == "v")
        {
            fin >> *a >> *b >> *c;
            v.push_back((Vector3){*a, *b, *c});
        }
        else if (v2 == "vt")
        {
            fin >> *a >> *b;
            vt.push_back((Vector2){*a, *b});
        }
        else if (v2 == "vn")
        {
            fin >> *a >> *b >> *c;
            vn.push_back((Vector3){*a, *b, *c});
        }
    }
}

void FluidDynamicsModel::loadObjectModel()
{
    float a, b, c;
    std::string aString, bString, cString;

    readObjectModelFile(&aString, &bString, &cString, &a, &b, &c);

    for (int i = 0; i < f.size(); i++)
    {
        std::vector<Vector3> additionToOTP;
        for (int k = 0; k < 3; k++)
        {
            char currentChar = 'e';
            int index = 0;
            int j = 0;
            std::vector<int> indexNumbers; // the sum of indexNumbers * pow(10, i) in reverse order of a for loop added together -1 is the index of the point
            while (currentChar != '/' && j < f.at(i).at(k).size())
            {
                currentChar = f.at(i).at(k).at(j);
                if (currentChar != '/')
                {
                    indexNumbers.push_back(currentChar - '0');
                }
                j++;
            }
            if (indexNumbers.size() >= 1)
            {
                int kHelper = 0;
                for (int k = indexNumbers.size() - 1; k >= 0; k--)
                {
                    index += indexNumbers.at(k) * pow(10, kHelper);
                    kHelper++;
                }
            }
            index--; // subtract one from index, because object files start index is at 1 and cpp index starts at 0
            additionToOTP.push_back({v.at(index).x, v.at(index).y, v.at(index).z});
        }
        objectTrianglePoints.push_back(additionToOTP);
    }
}

void FluidDynamicsModel::drawModel()
{
    if (objectTrianglePoints.size() >= 1)
    {
        for (int i = 0; i < objectTrianglePoints.size(); i++)
        {
            DrawTriangle3D(objectTrianglePoints.at(i).at(0), objectTrianglePoints.at(i).at(1), objectTrianglePoints.at(i).at(2), RED);
            DrawLine3D(objectTrianglePoints.at(i).at(0), objectTrianglePoints.at(i).at(1), BLACK);
        }
    }
}

void modelToCubes() {
    // start with a cube in the middle and go in every direction until collision with the object repeat for every cube;
}