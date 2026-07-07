//
//  recovery_parameters.hpp
//  ATCG
//

#ifndef recovery_parameters_hpp
#define recovery_parameters_hpp

#include <algorithm>
#include <vector>

class RecoveryParameters
{
public:
    explicit RecoveryParameters(int parameter_count = 39)
    {
        resize(parameter_count);
    }

    void resize(int parameter_count)
    {
        parameter_count_ = parameter_count;
        values_.assign(parameter_count_ + 1, 0.0);
    }

    void operator=(double value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }

    double &operator()(int index)
    {
        return values_[index];
    }

    double operator()(int index) const
    {
        return values_[index];
    }

private:
    int parameter_count_ = 0;
    std::vector<double> values_;
};

#endif /* recovery_parameters_hpp */
