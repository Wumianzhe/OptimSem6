#include "gradient.h"
#include "linmethods.h"
#include "vector.h"
#include <iostream>

using namespace std;

double f(vector_t x) { return x[0] * x[0] + 4 * x[1] * x[1] + sin(6 * x[0] + 7 * x[1]) + 3 * x[0] + 2 * x[1]; }

vector_t grad(vector_t x) {
    vector_t res(2);
    res[0] = 2 * x[0] + 6 * cos(6 * x[0] + 7 * x[1]) + 3;
    res[1] = 8 * x[1] + 7 * cos(6 * x[0] + 7 * x[1]) + 2;
    return res;
}

Matrix Hessian(vector_t x) {
    Matrix res(2, 2);
    res(0, 0) = -36 * sin(6 * x[0] + 7 * x[1]) + 2;
    res(1, 0) = res(0, 1) = -42 * sin(6 * x[0] + 7 * x[1]);
    res(1, 1) = -49 * sin(6 * x[0] + 7 * x[1]) + 8;
    return res;
}

int main(int argc, char* argv[]) {
    array<LinMethod*, 2> methods;
    methods[0] = new Ratio;
    methods[1] = new Dichotomy;
    // Julia gradient descent solution (global min)
    vector_t sol = {-1.9043886748, -0.3679467219};
    Grad1 gr(f, grad, sol);
    gr.write({1e-1, 1e-3}, *methods[0], {0, 0});
    gr.solve(f, grad, {0, 0}, 0.1, *methods[0], "d");
    delete methods[0];
    // local minimum near 0,0
    vector_t sol2 = {-4.050064e-1, 6.937312e-2};
    Grad1 gd(f, grad, sol2);
    gd.write({1e-1, 1e-3}, *methods[1], {0, 0});
    delete methods[1];

    Grad2 g2(f, grad, Hessian, sol2);
    g2.write({1e-1, 1e-3}, {0, -0});
    cout << (1 - 0.1) / 58.5 << "     "
         << "(1-eps)/R" << endl;
    g2.solve(f, grad, {0, -0}, 1e-7, "d");
    return 0;
}
