#include <iostream>
#include "linmethods.h"

using namespace std;

double f(double x) {
    return x / (1 + x*x);
}

int main(int argc, char* argv[]) {
    Ratio<double> rat;
    double x = rat.solve(f, -1.5, 1.5, 1e-6);
    cout << x << endl;
    Dichotomy<double> dich;
    x = dich.solve(f, -1.5, 1.5, 1e-6);
    cout << x << endl;
    TestPoints<double> test;
    x = test.solve(f, -1.5, 1.5, 1e-6);
    cout << x << endl;
    return 0;
}
