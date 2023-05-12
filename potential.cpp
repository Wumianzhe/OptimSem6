#include "potential.h"
#include <functional>
#include <iostream>
#include <sstream>
#include <stack>

enum dir { hor, vert };

using coord_t = std::pair<int, int>;
using namespace std;

vector<string> split(string line, char sep);
bool buildCycleInt(const Matrix& X, vector<coord_t>& cycle, coord_t fin, vector<coord_t>& cand, dir d);
transpTask readTransport(string filename) {
    fstream in(filename);
    int height, width;
    string line;
    // sizes
    getline(in, line);
    sscanf(line.c_str(), "%d,%d\n", &height, &width);
    Matrix C(width, height);
    vector_t a(height);
    vector_t b(width);
    vector_t p(width); // shorthand for penalties

    // most coefficients
    for (int i = 0; i < height; i++) {
        getline(in, line);
        auto words = split(line, ',');
        transform(words.begin(), words.begin() + width, C.begin() + i * width, [](string& str) { return stod(str); });
        a[i] = stod(words[width]);
    }

    // required amounts
    getline(in, line);
    vector<string> words = split(line, ',');
    transform(words.begin(), words.end(), b.begin(), [](string& str) { return stod(str); });

    // penalties
    getline(in, line);
    words = split(line, ',');
    if (words.empty()) {
        p = vector_t(width);
    } else {
        transform(words.begin(), words.end(), p.begin(), [](string& str) { return stod(str); });
    }

    // open to closed conversion
    int sum = 0;
    // providers
    for (int i = 0; i < height; i++) {
        sum += a[i];
    }
    // requirements
    for (int j = 0; j < width; j++) {
        sum -= b[j];
    }
    // more required than provided
    if (sum < 0) {
        C.resize(width, height + 1);
        a.resize(height + 1);
        for (int i = 0; i < width; i++) {
            // penalties
            C(i, height) = p[i];
        }
        // provide all missing
        a[height] = -sum;
    }

    C = C.T();

    // more provided than required
    if (sum > 0) {
        C.resize(height, width + 1);
        b.resize(width + 1);
        // require all surplus
        b[width] = sum;
        for (int i = 0; i < height; i++) {
            C(i, width) = 0;
        }
    }
    return {C, a, b, sum};
}

vector<string> split(string line, char sep) {
    stringstream iss(line);
    vector<string> ret;
    for (string word; getline(iss, word, sep);) {
        ret.push_back(word);
    }
    return ret;
}

// NW corner
void transpTask::buildInit() {
    int m = a.size();
    int n = b.size();
    auto tmpa = a;
    auto tmpb = b;
    int i = 0, j = 0;
    for (int k = 0; k < m + n - 1; k++) {
        int x_k = min(tmpa[i], tmpb[j]);
        tmpa[i] -= x_k;
        tmpb[j] -= x_k;
        X(i, j) = x_k;
        if (tmpa[i] == 0) {
            i++;
            continue;
        }
        j++;
        continue;
    }
}

void transpTask::solvePot() {
    while (!isOptimal()) {
        auto cycle = buildCycle();
        coord_t coords = cycle[1];
        int theta = X(coords.first, coords.second);
        // k = 1 is already "done"
        for (int k = 3; k < cycle.size(); k += 2) {
            if (X(cycle[k]) < theta) {
                coords = cycle[k];
                theta = X(cycle[k]);
            }
        }
        X(cycle[0]) = theta;
        int sign = -1;
        for (int k = 1; k < cycle.size(); k++) {
            X(cycle[k]) += theta * sign;
            sign *= -1;
        }
        X(coords) = -1;
    }
}

std::vector<std::pair<int, int>> transpTask::buildCycle() {
    int m = a.size();
    int n = b.size();
    vector<coord_t> cycle;
    auto [i, j] = dMaxCell.second;
    vector<coord_t> initCand;
    for (int i = 0; i < m; i++) {
        if (X(i, j) >= 0) {
            initCand.push_back({i, j});
        }
    }
    buildCycleInt(X, cycle, dMaxCell.second, initCand, dir::hor);
    cycle.insert(cycle.begin(), dMaxCell.second);
    dMaxCell = {0, {0, 0}};
    return cycle;
}

bool buildCycleInt(const Matrix& X, vector<coord_t>& cycle, coord_t fin, vector<coord_t>& cand, dir d) {
    int m = X.rows;
    int n = X.cols;
    vector<coord_t> nextCand;
    for (auto coords : cand) {
        auto [i, j] = coords;
        if (d == dir::vert) {
            // search in direction for end or for other candidates
            if (fin.second == j) {
                cycle.push_back(coords);
                return true;
            }
            for (int i = 0; i < m; i++) {
                if (coords != pair{i, j} && X(i, j) >= 0) {
                    nextCand.push_back({i, j});
                }
            }
        } else {
            if (fin.first == i) {
                cycle.push_back(coords);
                return true;
            }
            for (int j = 0; j < n; j++) {
                if (coords != pair{i, j} && X(i, j) >= 0) {
                    nextCand.push_back({i, j});
                }
            }
        }
        if (buildCycleInt(X, cycle, fin, nextCand, (d == dir::vert) ? dir::hor : dir::vert)) {
            cycle.push_back(coords);
            return true;
        }
    }
    return false;
}

bool transpTask::isOptimal() {
    int m = a.size();
    int n = b.size();
    vector_t v(n);
    vector_t u(m);

    stack<pair<coord_t, dir>> evals;
    for (int j = 0; j < n; j++) {
        if (X(0, j) >= 0) {
            evals.push({{0, j}, dir::vert});
        }
    }
    // пока остались нерешённые ограничения
    while (!evals.empty()) {
        auto [coords, dir] = evals.top();
        auto [i, j] = coords;
        evals.pop();
        if (dir == dir::vert) {
            v[j] = C(i, j) + u[i];
            for (int i = 0; i < m; i++) {
                if (coords != pair{i, j} && X(i, j) >= 0) {
                    evals.push({{i, j}, dir::hor});
                }
            }
        } else {
            u[i] = v[j] - C(i, j);
            for (int j = 0; j < n; j++) {
                if (coords != pair{i, j} && X(i, j) >= 0) {
                    evals.push({{i, j}, dir::vert});
                }
            }
        }
    }

    bool optim = true;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int delta = C(i, j) - v[j] + u[i];
            if (optim && delta < 0) {
                optim = false;
                if (abs(delta) >= dMaxCell.first) {
                    dMaxCell.first = abs(delta);
                    dMaxCell.second = {i, j};
                }
            }
        }
    }
    return optim;
}

int transpTask::results(const Matrix& Xr) {
    int m = a.size();
    int n = b.size();
    int cost = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cost += C(i, j) * ((Xr(i, j) < 0) ? 0 : Xr(i, j));
        }
    }
    return cost;
}

void transpTask::prettyPrint() {
    int m = a.size();
    int n = b.size();
    if (diff > 0) {
        n--;
    } else if (diff < 0) {
        m--;
    }
    Matrix Cf(m, n + 1);
    Matrix Xf(m, n);
    vector_t af(m);
    vector_t bf(n);
    copy(a.begin(), a.begin() + m, af.begin());
    copy(b.begin(), b.begin() + n, bf.begin());
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Cf(i, j) = C(i, j);
            Xf(i, j) = (X(i, j) < 0) ? 0 : X(i, j);
        }
    }
    copy(af.begin(), af.end(), Cf.begin() + n * m);
    cout << "Условия задачи:" << endl;
    cout << Cf << bf.T();
    if (diff < 0) {
        bool p = false;
        for (int j = 0; j < n; j++) {
            if (C(m, j) > 0) {
                p = true;
                break;
            }
        }
        if (p) {
            cout << "Штрафы: ";
            for (int j = 0; j < n - 1; j++) {
                cout << C(m, j) << ' ';
            }
            cout << C(m, n - 1) << endl;
        }
    }
    cout << "Значение функции цели: " << results(X) << endl;
    cout << "Матрица перевозок:" << endl;
    cout << Xf;
    if (diff < 0) {
        cout << "Недопоставки: ";
        for (int j = 0; j < n - 1; j++) {
            cout << ((X(m, j) < 0) ? 0 : X(m, j)) << ' ';
        }
        cout << X(m, n - 1) << endl;
    }
}
