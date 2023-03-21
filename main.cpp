#include <iostream>
#include "linmethods.h"

using namespace std;

double f(double x) {
    return x / (1 + x*x);
}

int main(int argc, char* argv[]) {
    Ratio<double> rat;
    array<LinMethod<double>*,3> methods;
    methods[0] = new Ratio<double>;
    methods[1] = new Dichotomy<double>;
    methods[2] = new TestPoints<double>;
    for (auto method: methods) {
        method->write({1e-3,1e-9}, f, -1.5, 1.5, 2, -1);
        delete method;
    }
    return 0;
}
