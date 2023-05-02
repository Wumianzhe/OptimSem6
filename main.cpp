#include <iostream>
#include <functional>
#include "simplex.h"
#include <sstream>

using namespace std;

std::tuple<task_t, set<int>, set<int>> read(string filename);
vector<string> split(string line, char sep);
vector_t restore(vector_t& xCan, int size, set<int>& unbound);

double f(vector_t x) { return x[0] * x[0] + 4 * x[1] * x[1] + sin(6 * x[0] + 7 * x[1]) + 3 * x[0] + 2 * x[1]; }

vector_t grad(vector_t x) {
    vector_t res(3);
    res[0] = 2 * x[0] + 6 * cos(6 * x[0] + 7 * x[1]) + 3;
    res[1] = 8 * x[1] + 7 * cos(6 * x[0] + 7 * x[1]) + 2;
    res[2] = 0;
    return res;
}

int main(int argc, char* argv[]) {
    auto [task, unbound, noneq] = read("task.csv");
    task_t primCan = genCanon(task,unbound,noneq,extrem::min);

    auto init = initBasic(primCan);
    auto x0 = simplex(primCan,init,{});
    auto x0r = restore(x0, task.C.size(), unbound);
    // auto [dual, dUnbound, dNoneq] = genDual(task, unbound, noneq);
    // auto dualCan = genCanon(dual,dUnbound,dNoneq,extrem::max);
    // dualCan.C = dualCan.C *-1;
    // cout << dualCan.C.T() << endl << dualCan.A << endl;
    // auto dInit = initBasic(dualCan);
    // auto dualCanSol = simplex(dualCan,dInit,{});
    // cout << dualCanSol.T();
    // auto dualSol = restore(dualCanSol, dual.C.size(), dUnbound);

    // Matrix res19 = task.C.T() - dualSol.T() * task.A;
    // cout << res19;
    // cout << res19 * x0r;

    vector_t xk = x0r;
    vector_t xp(3);
    do {
        vector_t ak = grad(xk);
        double bkn = -f(xk) + dot(ak,xk);
        Matrix Ak(task.A.rows+1,task.A.cols);
        vector_t bk(task.b.rows+1);
        for (int i=0; i < task.A.rows; i++) {
            for (int j=0; j < task.A.cols; j++) {
                Ak(i,j) = task.A(i,j);
            }
        }
        for (int j=0; j < task.A.cols; j++) {
            Ak(task.A.rows,j) = -ak[j];
        }
        for (int j=0; j < task.b.rows; j++) {
            bk[j] = -task.b[j];
        }
        bk[task.b.rows] = bkn;
        task.A = Ak;
        task.b = bk;
        xp = xk;
        noneq.insert(task.A.rows-1);

        task_t can = genCanon(task,unbound,noneq,extrem::min);
        auto initk = initBasic(can);
        auto xkraw = simplex(can,initk,{});
        xk = restore(xkraw,can.C.size(), unbound);
    } while (norm(xk - xp) > 1e-3);
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
