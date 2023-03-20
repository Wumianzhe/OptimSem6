#ifndef LINMETHODS_H_
#define LINMETHODS_H_
#include <functional>
#include <string>

template <typename T>
class LinMethod {
    public:
    // returns l: a + l*p minimizes f
    virtual double solve(std::function<double(T)> f, T a, T p, double eps);
    void write(std::pair<double,double> eps_bounds, std::function<double(T)> f, T a, T p, double fact);
    protected:
    std::string m_name;
};

template <typename T>
class Ratio : LinMethod<T> {
    double solve(std::function<double(T)> f, T a, T p, double eps) override;
};

template <typename T>
class Dichotomy : LinMethod<T> {
    double solve(std::function<double(T)> f, T a, T p, double eps) override;
};

template <typename T>
class TestPoints : LinMethod<T> {
    double solve(std::function<double(T)> f, T a, T p, double eps) override;
};
#endif // LINMETHODS_H_
