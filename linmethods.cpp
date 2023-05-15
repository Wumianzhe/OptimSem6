#include "linmethods.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

double Ratio::solve(function<double(double)> f, double a, const double p_0, const double eps) const {
    const double alpha = (3 - sqrt(5)) / 2;
    double p = p_0;
    double lambda = a + p * alpha;
    double mu = a + p * (1 - alpha);
    array<double, 2> vals;
    vals[0] = f(lambda);
    vals[1] = f(mu);
    while (abs(lambda - mu) > eps * (1 - 2 * alpha)) {
        if (vals[0] > vals[1]) {
            p *= 1 - alpha;
            a = lambda;
            double n_mu = a + p * (1 - alpha);
            lambda = mu;
            mu = n_mu;
            vals[0] = vals[1];
            vals[1] = f(mu);
        } else {
            p *= 1 - alpha;
            double n_l = a + p * alpha;
            mu = lambda;
            lambda = n_l;
            vals[1] = vals[0];
            vals[0] = f(lambda);
        }
    }
    return lambda;
}

double Dichotomy::solve(function<double(double)> f, double a, const double p_0, const double eps) const {
    double p = p_0;
    double delta = p * 1e-3;
    double x1 = a + p / 2 - delta;
    double x2 = a + p / 2 + delta;
    array<double, 2> vals;
    vals[0] = f(x1);
    vals[1] = f(x2);
    while (abs(a - x2) > eps) {
        if (vals[1] > vals[0]) {
            p = p / 2 + delta;
            delta = p * 1e-3;
            x1 = a + p / 2 - delta;
            x2 = a + p / 2 + delta;
            vals[0] = f(x1);
            vals[1] = f(x2);
        } else {
            p = p / 2 + delta;
            a = x1;
            delta = p * 1e-3;
            x1 = a + p / 2 - delta;
            x2 = a + p / 2 + delta;
            vals[0] = f(x1);
            vals[1] = f(x2);
        }
    }
    return (x1 + x2) / 2;
}

double TestPoints::solve(function<double(double)> f, double a, const double p_0, const double eps) const {
    double p = p_0;
    array<double, 3> x;
    array<double, 3> vals;
    for (int i = 0; i < 3; i++) {
        x[i] = a + p / 4 * (i + 1);
        vals[i] = f(x[i]);
    }
    while (abs(x[0] - x[2]) > eps) {
        if (vals[0] < vals[1]) {
            p = p / 2;
            x[1] = x[0];
            vals[1] = vals[0];
            x[0] = x[1] - p / 4;
            x[2] = x[1] + p / 4;
            vals[0] = f(x[0]);
            vals[2] = f(x[2]);
        } else {
            if (vals[1] < vals[2]) {
                p = p / 2;
                x[0] = x[1] - p / 4;
                x[2] = x[1] + p / 4;
                vals[0] = f(x[0]);
                vals[2] = f(x[2]);
            } else {
                p = p / 2;
                x[1] = x[2];
                vals[1] = vals[2];
                x[0] = x[1] - p / 4;
                x[2] = x[1] + p / 4;
                vals[0] = f(x[0]);
                vals[2] = f(x[2]);
            }
        }
    }
    return x[1];
}

const int n=4;
double UniSearch::solve(std::function<double (double)> f, double a, const double p_0, const double eps) const {
    double p = p_0;
    array<double, n> x;
    array<double, n> vals;
    int iMin = 0;
    do {
        for (int i = 0; i < n; i++) {
            x[i] = a + p / n * (i + 1);
            vals[i] = f(x[i]);
        }
        iMin = 0;
        for (int i=1; i < n; i++) {
            if (vals[i] < vals[iMin]) {
                iMin = i;
            }
        }
        a = (iMin>0)?x[iMin-1]:a;
        p = p * 2 / (n+1);
    } while (p > eps);
    return x[iMin];
}
