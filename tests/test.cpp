#include <iostream>
#include <fstream>
#include <json/json.h>

using namespace std;

int main() {
    ifstream ifs("CloverAlpha.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj); // reader can also read strings
    cout << "Book: " << obj["-8.750"]["Cl"] << endl;
    cout << "Year: " << obj["-8.750"]["Cd"] << endl;
    
    return 0;
}