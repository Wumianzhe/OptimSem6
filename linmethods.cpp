#include "linmethods.h"
#include <cmath>

using namespace std;

template <typename T> T Ratio<T>::solve(function<double(T)> f, T a, const T p_0, const double eps) {
    const double alpha = (3 - sqrt(5)) / 2;
    T p = p_0;
    T lambda = a + p * alpha;
    T mu = a + p * (1 - alpha);
    array<double, 2> vals;
    vals[0] = f(lambda);
    vals[1] = f(mu);
    while (abs(vals[1] - vals[0]) > eps) {
        if (vals[0] > vals[1]) {
            p *= 1 - alpha;
            a = lambda;
            T n_mu = a + p * (1 - alpha);
            lambda = mu;
            mu = n_mu;
            vals[0] = vals[1];
            vals[1] = f(mu);
        } else {
            p *= 1 - alpha;
            T n_l = a + p * alpha;
            mu = lambda;
            lambda = n_l;
            vals[1] = vals[0];
            vals[0] = f(lambda);
        }
    }
    return lambda;
}

template <typename T> T Dichotomy<T>::solve(function<double(T)> f, T a, const T p_0, const double eps) {
    T p = p_0;
    T delta = p * 1e-2;
    T x1 = a + p / 2 - delta;
    T x2 = a + p / 2 + delta;
    array<double, 2> vals;
    vals[0] = f(x1);
    vals[1] = f(x2);
    while (abs(vals[0] - vals[1]) > eps) {
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
    return x1;
}

template <typename T> T TestPoints<T>::solve(function<double(T)> f, T a, const T p_0, const double eps) {
    T p = p_0;
    array<double, 3> x;
    array<double, 3> vals;
    for (int i = 0; i < 3; i++) {
        x[i] = a + p / 4 * (i + 1);
        vals[i] = f(x[i]);
    }
    while (true) {
        if (abs(vals[0] - vals[1]) < eps) {
            return x[0];
        }
        if (vals[0] < vals[1]) {
            p = p / 2;
            x[1] = x[0];
            vals[1] = vals[0];
            x[0] = x[1] - p / 4;
            x[2] = x[1] + p / 4;
            vals[0] = f(x[0]);
            vals[2] = f(x[2]);
        } else {
            if (abs(vals[1] - vals[2]) < eps) {
                return x[1];
            }
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
}

template class Ratio<double>;
template class Dichotomy<double>;
template class TestPoints<double>;
