
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>

using namespace std;

int main() {
    ostringstream stream;
    for (int i=0; i<5;i++) {
        stream << i;
    }
    string output = stream.str();
    cout << output << endl;
}