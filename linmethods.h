#ifndef LINMETHODS_H_
#define LINMETHODS_H_
#include <functional>
#include <string>

class LinMethod {
  public:
    // returns l: a + l*p minimizes f
    virtual double solve(std::function<double(double)> f, double a, const  double p, const double eps) const = 0;
    LinMethod() = delete;
    virtual ~LinMethod(){};
    LinMethod(std::string name) : m_name(name) {}
    std::string m_name;
};

 class Ratio : public LinMethod {
  public:
    double solve(std::function<double(double)> f, double a, const  double p, const double eps) const override;
    Ratio() : LinMethod("ratio") {}
};

 class Dichotomy : public LinMethod {
  public:
    double solve(std::function<double(double)> f, double a, const  double p, const double eps) const override;
    Dichotomy() : LinMethod("dichotomy") {}
};

 class TestPoints : public LinMethod {
  public:
    double solve(std::function<double(double)> f, double a, const  double p, const double eps) const override;
    TestPoints() : LinMethod("testPoints") {}
};

#endif // LINMETHODS_H_
