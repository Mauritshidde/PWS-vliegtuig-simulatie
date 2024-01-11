// C++ program to demonstrate
// multithreading using three
// different callables.
#include <iostream>
#include <thread>
using namespace std;

class TestThis
{
private:
    /* data */
public:
    int var;
    TestThis(/* args */);
    ~TestThis();
    void Start(int i);
};

TestThis::TestThis(/* args */)
{
}

void TestThis::Start( int i) {
    std::cout << "test test" << i << std::endl;
    var = i;
}

TestThis::~TestThis()
{
}


void func() {
    for (int i=0; i < 1000000; i++) {
        std::cout << "1 ";
    }
    std::cout << std::endl;
}

void func2(int i) {
    for (int i=0; i < 1000000; i++) {
        std::cout << "2 ";
    }
    std::cout << std::endl;
}

int main()
{
    // Start thread t1
    // std::thread t1(func);
    TestThis test;
    std::thread t2(&TestThis::Start, &test, 3);
    // func();
    // Wait for t1 to finish

    // t1.join();
    t2.join();

    std::cout << test.var << std::endl;
    
 
    // t1 has finished do other stuff
}
