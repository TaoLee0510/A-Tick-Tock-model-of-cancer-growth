#ifndef recovery_parameters_hpp
#define recovery_parameters_hpp

#include <vector>

class RecoveryParameters
{
public:
    explicit RecoveryParameters(int parameter_count = 39);

    void resize(int parameter_count);
    void operator=(double value);
    double &operator()(int index);
    double operator()(int index) const;

private:
    int parameter_count_ = 0;
    std::vector<double> values_;
};

#endif
