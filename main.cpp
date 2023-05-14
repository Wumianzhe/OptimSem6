#include <iostream>
#include <functional>
#include "simplex.h"
#include <sstream>

using namespace std;

const double eps = 1e-6;

std::tuple<task_t, set<int>, set<int>> read(string filename);
vector<string> split(string line, char sep);
vector_t restore(vector_t& xCan, int size, set<int>& unbound);
vector_t naiveSolver(const vector_t& x0, const task_t& task, set<int>unbound, set<int> noneq);
vector_t fastSolver(const task_t& task, set<int> unbound, set<int> noneq);

double f(vector_t x) { return x[0] * x[0] + 4 * x[1] * x[1] + sin(6 * x[0] + 7 * x[1]) + 3 * x[0] + 2 * x[1]; }
double phi(vector_t x) { return x[0] * x[0] + 4 * x[1] * x[1] + sin(6 * x[0] + 7 * x[1]) + 3 * x[0] + 2 * x[1] -x[2]; }
double c(vector_t x) { return (x[0] + 2)*(x[0] + 2) + 3 * (x[1] + 0.5) * (x[1] + 0.5) - 0.5; }

vector_t gradf(vector_t x) {
    vector_t res(3);
    res[0] = 2 * x[0] + 6 * cos(6 * x[0] + 7 * x[1]) + 3;
    res[1] = 8 * x[1] + 7 * cos(6 * x[0] + 7 * x[1]) + 2;
    res[2] = -1;
    return res;
}

vector_t gradc(vector_t x) {
    vector_t res(3);
    res[0] = 2* (x[0] + 2);
    res[1] = 6* (x[1] + 0.5);
    res[2] = 0;
    return res;
}

int main(int argc, char* argv[]) {
    auto [task, unbound, noneq] = read("task.csv");
    task_t primCan = genCanon(task,unbound,noneq,extrem::min);
    vector_t sol = {-1.9043886748, -0.3679467219};
    // cout << primCan.A;

    auto init = initBasic(primCan);
    auto x0 = simplex(primCan,init,{});
    auto x0r = restore(x0, task.C.size(), unbound);

    vector_t xf = fastSolver(task, unbound, noneq);
    cout << xf.T();
    cout << norm(sol-xf) << endl;

    vector_t x = naiveSolver(x0r, task, unbound, noneq);
    cout << norm(sol-x) << endl;
    return 0;
}

std::tuple<task_t, set<int>, set<int>> read(string filename) {
    fstream in(filename);
    int height, width;
    string line;
    // sizes
    getline(in, line);
    sscanf(line.c_str(), "%d,%d\n", &height, &width);
    vector_t C(width);
    Matrix A(width, height);
    vector_t b(height);
    // target function
    getline(in, line);
    auto funCoefs = split(line, ',');
    transform(funCoefs.begin(), funCoefs.end(), C.begin(), [](string& str) { return stod(str); });
    // reading coefs
    set<int> noneqIndices;
    for (int i = 0; i < height; i++) {
        getline(in, line);
        auto words = split(line, ',');
        transform(words.begin(), words.begin() + width, A.begin() + i * width, [](string& str) { return stod(str); });
        b[i] = stod(words[width + 1]);
        if (words[width] != "EQ") {
            noneqIndices.insert(i);
            if (words[width] == "LT") {
                transform(A.begin() + i * width, A.begin() + (i + 1) * width, A.begin() + i * width,
                          [](double d) { return -1 * d; });
                b[i] *= -1;
            }
        }
    }
    // unbounded variables
    getline(in, line);
    auto words = split(line, ',');
    set<int> unboundIndices;
    transform(words.begin(), words.end(), inserter(unboundIndices, unboundIndices.end()),
              [](string& str) { return stoi(str) - 1; });
    // I almost forgot i store in column-major, so consecutive areas are columns, not rows
    return {{C, A.T(), b}, unboundIndices, noneqIndices};
}

vector<string> split(string line, char sep) {
    stringstream iss(line);
    vector<string> ret;
    for (string word; getline(iss, word, sep);) {
        ret.push_back(word);
    }
    return ret;
}

vector_t restore(vector_t& xCan, int size, set<int>& unbound) {
    vector_t res(size);
    copy(xCan.begin(), xCan.begin() + size, res.begin());
    for (int i : unbound) {
        int ord = distance(unbound.begin(), find(unbound.begin(), unbound.end(), i));
        if (xCan[size + ord] > 0) {
            res[i] -= xCan[size + ord];
        }
    }
    return res;
}

vector_t naiveSolver(const vector_t& x0, const task_t& task, set<int>unbound, set<int> noneq) {
    // 8 seconds of real time on PC
    // I see why alternative is needed
    task_t taskInt = task;
    vector_t xk = x0;
    vector_t xp(3);
    int iter = 0;
    do {
        vector_t ak = (phi(xk) < c(xk))?gradc(xk):gradf(xk);
        double bkn = -((phi(xk) < c(xk))?c(xk):phi(xk)) + dot(ak,xk);
        Matrix Ak(taskInt.A.rows+1,taskInt.A.cols);
        vector_t bk(taskInt.b.rows+1);
        for (int i=0; i < taskInt.A.rows; i++) {
            for (int j=0; j < taskInt.A.cols; j++) {
                Ak(i,j) = taskInt.A(i,j);
            }
        }
        for (int j=0; j < taskInt.A.cols; j++) {
            Ak(taskInt.A.rows,j) = -ak[j];
        }
        for (int j=0; j < taskInt.b.rows; j++) {
            bk[j] = taskInt.b[j];
        }
        bk[taskInt.b.rows] = -bkn;
        taskInt.A = Ak;
        taskInt.b = bk;
        xp = xk;
        noneq.insert(taskInt.A.rows-1);

        task_t can = genCanon(taskInt,unbound,noneq,extrem::min);
        auto initk = initBasic(can);
        auto xkraw = simplex(can,initk,{});
        xk = restore(xkraw,taskInt.C.size(), unbound);
        cout << "Итерация: " <<++iter << ", значение: ";
        cout << xk.T();
    } while (norm(xk - xp) > eps);
    vector_t sol(2);
    sol[0] = xk[0];
    sol[1] = xk[1];
    return sol;
}

vector_t fastSolver(const task_t& task, set<int> unbound, set<int> noneq) {
    int iter = 0;

    auto [dual, dUnbound, dNoneq] = genDual(task, unbound, noneq);
    dual.C = dual.C *-1;
    auto dInit = initBasic(dual);
    auto dualSol = simplex(dual,dInit,{});
    set<int> dualNoneq;
    for (int i=0; i < dualSol.size(); i++) {
        if (dualSol[i] !=0) {
            dualNoneq.insert(i);
        }
    }
    vector_t xk = solve(dual.A[dualNoneq].T(),dual.C.T()[dualNoneq].T()*-1);
    vector_t xp(3);
    do {
        vector_t ak = ((phi(xk) < c(xk))?gradc(xk):gradf(xk))*-1;
        double bk = ((phi(xk) < c(xk))?c(xk):phi(xk)) + dot(ak,xk);

        Matrix Ak = dual.A;
        Ak.resize(dual.A.rows,dual.A.cols+1); // appending to the right is easy
        vector_t Ck = dual.C;
        Ck.resize(dual.C.rows+1);
        vector_t init = dualSol;
        init.resize(dual.C.rows+1);

        copy(ak.begin(),ak.end(),Ak.begin() + dual.A.size());
        Ck[dual.C.rows] = -bk;
        init[dual.C.rows] = 0;
        dualSol = simplex({Ck,Ak,dual.b},init,{});
        set<int> dualNoneq;
        for (int i = 0; i < dualSol.size(); i++) {
            if (dualSol[i] != 0) {
                dualNoneq.insert(i);
            }
        }

        dual.A = Ak;
        dual.C = Ck;
        xp = xk;
        xk = solve(Ak[dualNoneq].T(),Ck.T()[dualNoneq].T()*-1);
        cout << "Итерация: " <<++iter << ", значение: ";
        cout << xk.T();
    } while (norm(xk - xp) > eps);

    vector_t sol(2);
    sol[0] = xk[0];
    sol[1] = xk[1];
    return sol;
}
