#include "matplotlibcpp.h"
namespace mat = matplotlibcpp;
int main()
{
    // u and v are respectively the x and y components of the arrows we're plotting
    std::vector<int> x, y, u, v;
    for (int i = -5; i <= 5; i++) {
        for (int j = -5; j <= 5; j++) {
            x.push_back(i);
            u.push_back(-i);
            y.push_back(j);
            v.push_back(-j);
        }
    }
    // std::cout << "{";
    // for (int i = 0; i < x.size(); i++)
    // {
    //     std::cout << x.at(i) << ", ";
    // }
    // std::cout << "} \n";

    // std::cout << "{";
    // for (int i = 0; i < y.size(); i++)
    // {
    //     std::cout << y.at(i) << ", ";
    // }
    // std::cout << "} \n";

    std::cout << "{";
    for (int i = 0; i < u.size(); i++)
    {
        std::cout << u.at(i) << ", ";
    }
    std::cout << "} \n";

    std::cout << "{";
    for (int i = 0; i < v.size(); i++)
    {
        std::cout << v.at(i) << ", ";
    }
    std::cout << "} \n";

    mat::quiver(x, y, u, v);
    mat::show();
}