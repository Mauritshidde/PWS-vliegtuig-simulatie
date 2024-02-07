#include "file.h"

void func() {
    std::string myText;

    // Read from the text file
    std::ifstream MyReadFile("planes/liftfiles/BoeingVLMlift.txt");

    // Use a while loop together with the getline() function to read the file line by line
    std::vector<std::string> vals;

    while (getline (MyReadFile, myText)) {
        vals.push_back(myText);
    }

    std::vector<std::vector<double>> lift;

    for (int i=0; i < vals.size(); i++) {
        std::vector<std::string> result;
        std::stringstream s_stream(vals.at(i)); //create string stream from the string
        while(s_stream.good()) {
            std::string substr;
            getline(s_stream, substr, ','); //get first string delimited by comma
            result.push_back(substr);
        }
        for(int i = 0; i<result.size(); i++) {    //print all splitted strings
            std::cout << std::stod(result.at(i)) << std::endl;
        }

        // lift.push_back(result);
    }


    // Close the file
    MyReadFile.close(); 
}