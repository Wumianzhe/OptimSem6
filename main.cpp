#include <iostream>
#include "potential.h"

using namespace std;

int main(int argc, char* argv[]) {
    transpTask task = readTransport("task.csv");
    cout << task.C;
    task.buildInit();
    cout << task.X;
    return 0;
}
