#ifndef LINMETHODS_H_
#define LINMETHODS_H_
#include <functional>
#include <string>

template <typename T> class LinMethod {
  public:
    // returns l: a + l*p minimizes f
    virtual T solve(std::function<double(T)> f, T a, const T p, const double eps) = 0;
    LinMethod() = delete;
    virtual ~LinMethod(){};
    void write(std::pair<double, double> eps_bounds, std::function<double(T)> f, T a, T p, double fact);

  protected:
    LinMethod(std::string name) : m_name(name) {}
    std::string m_name;
};

template <typename T> class Ratio : public LinMethod<T> {
  public:
    T solve(std::function<double(T)> f, T a, const T p, const double eps) override;
    Ratio() : LinMethod<T>("ratio") {}
};

template <typename T> class Dichotomy : public LinMethod<T> {
  public:
    T solve(std::function<double(T)> f, T a, const T p, const double eps) override;
    Dichotomy() : LinMethod<T>("dichotomy") {}
};

template <typename T> class TestPoints : public LinMethod<T> {
  public:
    T solve(std::function<double(T)> f, T a, const T p, const double eps) override;
    TestPoints() : LinMethod<T>("testPoints") {}
};

#endif // LINMETHODS_H_
