#include "recovery_parameters.hpp"

#include <algorithm>

RecoveryParameters::RecoveryParameters(int parameter_count)
{
    resize(parameter_count);
}

void RecoveryParameters::resize(int parameter_count)
{
    parameter_count_ = parameter_count;
    values_.assign(parameter_count_ + 1, 0.0);
}

void RecoveryParameters::operator=(double value)
{
    std::fill(values_.begin(), values_.end(), value);
}

double &RecoveryParameters::operator()(int index)
{
    return values_[index];
}

double RecoveryParameters::operator()(int index) const
{
    return values_[index];
}
