#include "linmethods.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include "vector.h"

using namespace std;

double abs(const vector_t& v) {
    return norm(v);
}

template <typename T>
void LinMethod<T>::write(std::pair<double, double> eps_bounds, std::function<double(T)> f, T a, T p, double fact, T x) {
    array<T, 3> approx;
    array<double, 3> epses;
    array<int, 3> counts;
    epses[0] = eps_bounds.first;
    epses[2] = eps_bounds.second;
    // sort epsilons
    if (epses[0] < epses[2]) {
        double tmp = epses[0];
        epses[0] = epses[2];
        epses[2] = tmp;
    }
    epses[1] = sqrt(epses[0] * epses[2]);

    int counter = 0;
    auto f_count = [&](T x) {
        counter++;
        return f(x);
    };

    approx[0] = solve(f_count, a, p, epses[0]);
    counts[0] = counter;

    counter = 0;
    approx[1] = solve(f_count, a, p, epses[1]);
    counts[1] = counter;

    counter = 0;
    approx[2] = solve(f_count, a, p, epses[2]);
    counts[2] = counter;

    cout << m_name << endl;
    cout << std::scientific;
    for (int i = 0; i < 3; i++) {
        cout << "eps: " << epses[i] << " res: " << approx[i] << " calls: " << counts[i] << endl;
        cout << "     delta: " << abs(approx[i] - x) << " fdelta: " << abs(f(approx[i]) - f(x));
        cout << endl;
    }
    cout << endl;

    ofstream file("res/" + m_name);
    file << "eps,delta,fdelta,iters" << endl;
    file << std::scientific;
    for (double eps = epses[0]; eps > epses[2]; eps /= fact) {
        counter = 0;
        T approx = solve(f_count, a, p, eps);
        file << eps << "," << max(1e-16, abs(approx - x)) << "," << max(1e-16, abs(f(approx) - f(x))) << "," << counter
             << endl;
    }
}

template <typename T> T Ratio<T>::solve(function<double(T)> f, T a, const T& p_0, const double eps) {
    const double alpha = (3 - sqrt(5)) / 2;
    T p = p_0;
    T lambda = a + p * alpha;
    T mu = a + p * (1 - alpha);
    array<double, 2> vals;
    vals[0] = f(lambda);
    vals[1] = f(mu);
    while (abs(lambda - mu) > eps) {
        if (vals[0] > vals[1]) {
            p = p * (1 - alpha);
            a = lambda;
            T n_mu = a + p * (1 - alpha);
            lambda = mu;
            mu = n_mu;
            vals[0] = vals[1];
            vals[1] = f(mu);
        } else {
            p = p * (1 - alpha);
            T n_l = a + p * alpha;
            mu = lambda;
            lambda = n_l;
            vals[1] = vals[0];
            vals[0] = f(lambda);
        }
    }
    return lambda;
}

template <typename T> T Dichotomy<T>::solve(function<double(T)> f, T a, const T& p_0, const double eps) {
    T p = p_0;
    T delta = p * 1e-3;
    T x1 = a + p / 2 - delta;
    T x2 = a + p / 2 + delta;
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
    return x1;
}

template <typename T> T TestPoints<T>::solve(function<double(T)> f, T a, const T& p_0, const double eps) {
    T p = p_0;
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

template class Ratio<vector_t>;
template class Dichotomy<vector_t>;
