#include "gradient.h"

using namespace std;
Matrix LinEqsolve(Matrix A, Matrix b);

vector_t Grad1::solve(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, const vector_t& x_0,
                      const double eps, const LinMethod& method, std::string of) {
    vector_t xk = x_0;
    vector_t xp = x_0;
    ofstream ofstr;
    if (of != "") {
        ofstr.open(of);
        cout << "Вывод промежуточных значений для eps = " << eps << endl;
        // cout << "||∇ f(xk)||^2, 4R{f(xk)-f(xk+1)}, q" << endl;
    }
    vector_t g_x = grad(xk) * -1;
    int iterations = 0;
    while (norm(g_x) > eps) {
        if (iterations++ > 1000) {
            break;
        }
        if (ofstr.is_open()) {
            ofstr << xk.T();
            // cout << norm(g_x)*norm(g_x)  << "   " << 4 * 58.5 * (m_f(xp) - m_f(xk)) << "   " << norm(xk -
            // m_sol)/norm(xp - m_sol) << endl;
            cout << norm(xk - m_sol) << "      ||x_k-x*||" << endl;
        }
        const double epsInner = 0.01;
        auto fa = [&](double a) { return f(xk + g_x * a); };
        double alpha = method.solve(fa, 0, 1, epsInner);
        xp = xk;
        xk = xk + g_x * alpha;
        g_x = grad(xk) * -1;
    }
    return xk;
}

void Grad1::write(std::pair<double, double> eps_bounds, const LinMethod& method, const vector_t& x_0) {
    vector<vector_t> approx(3, vector_t(2));
    array<double, 3> epses;
    array<pair<int, int>, 3> counts;
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
        return m_f(x);
    };
    int gcounter = 0;
    auto g_count = [&](vector_t x) {
        gcounter++;
        return m_grad(x);
    };

    for (int i = 0; i < 3; i++) {
        fcounter = gcounter = 0;
        approx[i] = solve(f_count, g_count, x_0, epses[i], method);
        counts[i] = {fcounter, gcounter};
    }

    cout << "Gradient descent using " << method.m_name << endl;
    cout << std::scientific;
    for (int i = 0; i < 3; i++) {
        cout << "eps: " << epses[i] << " res: " << approx[i].T() << " function calls: " << counts[i].first << endl;
        cout << "     delta: " << norm(approx[i] - m_sol) << " fdelta: " << abs(m_f(approx[i]) - m_f(m_sol))
             << " gradient norm: " << norm(m_grad(approx[i]));
        cout << endl;
    }
    cout << endl;

    solve(m_f, m_grad, x_0, epses[1], method, "res/" + method.m_name);
}

vector_t Grad2::solve(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, const vector_t& x_0,
                      const double eps, std::string of) {
    vector_t xk = x_0;
    vector_t xp = x_0;
    ofstream ofstr;
    if (of != "") {
        ofstr.open(of);
        cout << "Вывод промежуточных значений для eps = " << eps << endl;
    }
    vector_t g_x = grad(xk);
    int iterations = 0;
    while (norm(g_x) > eps) {
        if (iterations++ > 1000) {
            break;
        }
        if (ofstr.is_open()) {
            ofstr << xk.T();
            // cout << norm(g_x) << " " << xk.T();
            cout << norm(xk - m_sol) / (norm(xp - m_sol) * norm(xp - m_sol)) << "   " << 58.5 / 37.5 << " R/m" << endl;
        }
        vector_t p = LinEqsolve(m_Hessian(xk), g_x * -1);
        double alpha = 1;
        double e = 1.0 / 3;
        double l = 2.0 / 3;
        while (f(xk + p * alpha) - f(xk) > alpha * e * dot(g_x, p)) {
            alpha *= l;
        }
        xp = xk;
        xk = xk + p * alpha;
        g_x = grad(xk);
    }
    cout << norm(g_x) << " " << xk.T();
    return xk;
}

void Grad2::write(std::pair<double, double> eps_bounds, const vector_t& x_0) {
    vector<vector_t> approx(3, vector_t(2));
    array<double, 3> epses;
    array<pair<int, int>, 3> counts;
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
        return m_f(x);
    };
    int gcounter = 0;
    auto g_count = [&](vector_t x) {
        gcounter++;
        return m_grad(x);
    };

    for (int i = 0; i < 3; i++) {
        fcounter = gcounter = 0;
        approx[i] = solve(f_count, g_count, x_0, epses[i]);
        counts[i] = {fcounter, gcounter};
    }

    cout << "Newton method" << endl;
    cout << std::scientific;
    for (int i = 0; i < 3; i++) {
        cout << "eps: " << epses[i] << " res: " << approx[i].T() << " function calls: " << counts[i].first << endl;
        cout << "     delta: " << norm(approx[i] - m_sol) << " fdelta: " << abs(m_f(approx[i]) - m_f(m_sol))
             << " gradient norm: " << norm(m_grad(approx[i]));
        cout << endl;
    }
    cout << endl;

    solve(m_f, m_grad, x_0, epses[1], "res/second");
}

Matrix LinEqsolve(Matrix A, Matrix b) {
    int m = A.rows;
    auto [Qt, R] = Matrix::QtRdecomp(A);
    b = Qt * b;
    // gaussBack
    vector_t res(m);
    for (int i = m - 1; i >= 0; i--) {
        double s = 0;
        for (int k = i + 1; k < m; k++) {
            s += R(i, k) * res[k];
        }
        double x = (b(i, 0) - s) / R(i, i);
        res[i] = x;
    }
    return res;
}
