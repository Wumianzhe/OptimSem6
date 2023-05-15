#include "gradient.h"
#include "linmethods.h"
#include "vector.h"
#include <iostream>

using namespace std;

double f(vector_t x) { return x[0] * x[0] + 4 * x[1] * x[1] + sin(6 * x[0] + 3 * x[1]) + 3 * x[0] + 2 * x[1]; }

vector_t grad(vector_t x) {
    vector_t res(2);
    res[0] = 2 * x[0] + 6 * cos(6 * x[0] + 3 * x[1]) + 3;
    res[1] = 8 * x[1] + 3 * cos(6 * x[0] + 3 * x[1]) + 2;
    return res;
}

int main(int argc, char* argv[]) {
    array<LinMethod*, 2> methods;
    methods[0] = new Ratio;
    methods[1] = new UniSearch;
    // Julia gradient descent solution (global min)
    vector_t sol = {-1.217384, -0.21467307551};
    Grad1 gr(f, grad, sol);
    gr.write({1e-2, 1e-8}, *methods[0], {-0, -0});
    delete methods[0];
    Grad1 gu(f, grad, sol);
    gu.write({1e-2, 1e-8}, *methods[1], {-1, -0});
    delete methods[1];

    return 0;
}
