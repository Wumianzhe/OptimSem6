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
// phi uses x[2], which is additional variable reflecting value of f
double phi(vector_t x) { return x[0] * x[0] + 4 * x[1] * x[1] + sin(6 * x[0] + 7 * x[1]) + 3 * x[0] + 2 * x[1] -x[2]; }
double c(vector_t x) { return (x[0] + 2)*(x[0] + 2) + 3 * (x[1] + 0.5) * (x[1] + 0.5) - 0.5; }

vector_t gradphi(vector_t x) {
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
    vector_t sol = {-1.9043886748, -0.3679467219};
    // cout << primCan.A;

    vector_t xf = fastSolver(task, unbound, noneq);
    cout << xf.T();
    cout << norm(sol-xf) << endl;

    task_t primCan = genCanon(task,unbound,noneq,extrem::min);
    auto init = initBasic(primCan);
    auto x0 = simplex(primCan,init,{});
    auto x0r = restore(x0, task.C.size(), unbound);
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
    // 8 seconds of real time on PC for 1e-3 eps, 3 minutes for 1e-6
    // I see why alternative is needed
    task_t taskInt = task;
    vector_t xk = x0;
    vector_t xp(3);
    int iter = 0;
    do {
        // тернарный оператор для того, чтобы выбрать между двумя нелинейными ограничениями
        // как субградиент используется градиент большей из функций
        vector_t ak = (phi(xk) < c(xk))?gradc(xk):gradphi(xk);
        // то же самое для значения
        double bkn = -((phi(xk) < c(xk))?c(xk):phi(xk)) + dot(ak,xk);

        // копирование и добавление строчки
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

        //сохранение в форму условий задачи ЛП
        taskInt.A = Ak;
        taskInt.b = bk;
        noneq.insert(taskInt.A.rows-1);

        // следующая итерация
        xp = xk;
        task_t can = genCanon(taskInt,unbound,noneq,extrem::min);
        auto initk = initBasic(can);
        auto xkraw = simplex(can,initk,{});
        xk = restore(xkraw,taskInt.C.size(), unbound);

        cout << "Итерация: " <<++iter << ", значение: ";
        cout << xk.T();
    } while (norm(xk - xp) > eps);
    // возврат значения в виде двух аргументов
    vector_t sol(2);
    sol[0] = xk[0];
    sol[1] = xk[1];
    return sol;
}

vector_t fastSolver(const task_t& task, set<int> unbound, set<int> noneq) {
    int iter = 0;

    // генерация двойственной задачи
    auto [dual, dUnbound, dNoneq] = genDual(task, unbound, noneq);
    // симплекс только минимизирует
    dual.C = dual.C *-1;
    // начальное приближение и решение симплексом
    auto dInit = initBasic(dual);
    auto dualSol = simplex(dual,dInit,{});

    // восстановление решения по двойственной
    set<int> dualNoneq;
    for (int i=0; i < dualSol.size(); i++) {
        if (dualSol[i] !=0) {
            dualNoneq.insert(i);
        }
    }
    vector_t xk = solve(dual.A[dualNoneq].T(),dual.C.T()[dualNoneq].T()*-1);
    vector_t xp(3);
    do {
        // то же, что и раньше, но другие знаки.
        vector_t ak = ((phi(xk) < c(xk))?gradc(xk):gradphi(xk))*-1;
        double bk = (phi(xk) < c(xk))?c(xk):phi(xk) + dot(ak,xk);

        Matrix Ak = dual.A;
        Ak.resize(dual.A.rows,dual.A.cols+1); // appending to the right is easy
        vector_t Ck = dual.C;
        Ck.resize(dual.C.rows+1);
        vector_t init = dualSol;
        init.resize(dual.C.rows+1);

        // дополнение значениями
        copy(ak.begin(),ak.end(),Ak.begin() + dual.A.size());
        Ck[dual.C.rows] = -bk;
        init[dual.C.rows] = 0;

        // решение двойственной и восстановление
        dualSol = simplex({Ck,Ak,dual.b},init,{});
        set<int> dualNoneq;
        for (int i = 0; i < dualSol.size(); i++) {
            if (dualSol[i] != 0) {
                dualNoneq.insert(i);
            }
        }

        // сохранение и следующая итерация
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
