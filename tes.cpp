#include <vector>
#include <iostream>

void addForces(std::vector<double> inputVector)
{
    double sum = 0;
    for (int i = 0; i < inputVector.size(); i++)
    {
        sum += inputVector.at(i);
    }
    std::cout << sum << std::endl;
}

int main() {
    // int N = 10000000;
    // std::vector<double> test;
    // for (int i=0; i < N; i++) {
    //     test.push_back(1);
    // }
    // std::cout << "ja" << std::endl;
    // addForces(test);

    std::vector<std::vector<double>> test;
    for (int i=0; i < 10; i++) {
        std::vector<double> helper;
        for (int j=0; j < 10; j++) {
            helper.push_back(0);
        }
        test.push_back(helper);
    }

    // std::cout << test.at(1).at(3) << std::endl;
    // double *t = &test.at(1).at(3);
    // *t = 4;
    // std::cout << test.at(1).at(3) << std::endl;
    double *t = test.at(0).data();
    t[0] = 5;
    std::cout << t[0] << " " <<  test[0][0] << std::endl;
    // std::cout << test.data() << std::endl;

    return 0;
}