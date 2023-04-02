#include "gradient.h"
#include "linmethods.h"
#include "vector.h"
#include <iostream>

using namespace std;

double f(vector_t x) {
    return x[0]*x[0] + 4*x[1]*x[1] + sin(6*x[0] + 7 * x[1]) + 3*x[0] + 2*x[1];
}

vector_t grad(vector_t x) {
    vector_t res(2);
    res[0] = 2 * x[0] + 6 * cos(6 * x[0] + 7 * x[1]) + 3;
    res[1] = 8 * x[1] + 7 * cos(6 * x[0] + 7 * x[1]) + 2;
    return res;
}

int main(int argc, char* argv[]) {
    array<LinMethod<vector_t>*,3> methods;
    methods[0] = new Ratio<vector_t>;
    methods[1] = new Dichotomy<vector_t>;
    // Julia gradient descent
    vector_t sol = {-1.9043886748, -0.3679467219};
    for (auto method: methods) {
        Grad1 g(*method);
        g.write({1e-3,1e-9},f,grad,{0,0},sol);
        delete method;
    }
    return 0;
}
