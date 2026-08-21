#pragma once

#include <vector>
#include "sampled_positions.hpp"





// Equal-mass spherical sampler
SampledPositions sample_equal_mass(
    const std::vector<std::vector<double>>& gas_density,
    const std::vector<double>& Zcenter,
    const std::vector<double>& Rcenter,
    double target_mass,
    int NDisk
);


SampledPositions sample_cartesian(
    const std::vector<double>& LBox,
    const std::vector<double>& dx,
    const SampledPositions* initial = nullptr
);