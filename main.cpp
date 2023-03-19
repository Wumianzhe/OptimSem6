#include <iostream>
#include "potential.h"

using namespace std;

int main(int argc, char* argv[]) {
    transpTask task = readTransport("task.csv");
    cout << task.C;
    task.buildInit();
    task.isOptimal();
    cout << task.X;
    for (auto [i,j] : task.buildCycle()) {
        cout << "{" << i << "," << j <<"}\n";
    }
    return 0;
}
