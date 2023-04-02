#include "gradient.h"

using namespace std;

vector_t Grad1::solve(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, const vector_t& x_0, const double eps, std::string of) {
    vector_t xk = x_0;
    ofstream ofstr;
    if (of != "") {
        ofstr.open(of);
    }
    vector_t g_x = grad(xk)* -1;
    while (norm(g_x) > eps) {
        if (ofstr.is_open()) {
            ofstr << xk.T();
        }
        const double epsInner = eps/25;
        xk = linMin.solve(f, xk, g_x, epsInner);
        g_x = grad(xk)* -1;
    }
    return xk;
}

void Grad1::write(std::pair<double, double> eps_bounds, function<double (vector_t)> f, function<vector_t (vector_t)> grad, const vector_t& x_0, const vector_t& x) {
    vector<vector_t> approx(3,vector_t(2));
    array<double, 3> epses;
    array<pair<int,int>, 3> counts;
    epses[0] = eps_bounds.first;
    epses[2] = eps_bounds.second;
    // sort epsilons
    if (epses[0] < epses[2]) {
        double tmp = epses[0];
        epses[0] = epses[2];
        epses[2] = tmp;
    }
    epses[1] = sqrt(epses[0] * epses[2]);
    int fcounter = 0;
    auto f_count = [&](vector_t x) {
        fcounter++;
        return f(x);
    };
    int gcounter = 0;
    auto g_count = [&](vector_t x) {
        gcounter++;
        return grad(x);
    };
    approx[0] = solve(f_count, g_count, x_0, epses[0]);
    counts[0] = {fcounter,gcounter};

    fcounter = gcounter = 0;
    approx[1] = solve(f_count, g_count, x_0, epses[1]);
    counts[1] = {fcounter, gcounter};

    fcounter = gcounter = 0;
    approx[2] = solve(f_count, g_count, x_0, epses[2]);
    counts[2] = {fcounter, gcounter};

    cout << "Gradient descent using " << linMin.m_name << endl;
    cout << std::scientific;
    for (int i = 0; i < 3; i++) {
        cout << "eps: " << epses[i] << " res: " << approx[i].T() << " function calls: " << counts[i].first << " gradient calls: " << counts[i].second << endl;
        cout << "     delta: " << norm(approx[i] - x) << " fdelta: " << abs(f(approx[i]) - f(x)) << " gradient norm: " << norm(grad(approx[i]));
        cout << endl;
    }
    cout << endl;

    solve(f,grad,x_0,epses[1],"res/" + linMin.m_name);
}
