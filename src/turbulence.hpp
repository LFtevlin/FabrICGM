#ifndef TURBULENCE_HPP
#define TURBULENCE_HPP

#include <array>
#include <vector>
#include "sampled_positions.hpp"



void save_turbulence_field(
    const SampledPositions& field,
    const std::string& filename
);

std::vector<std::array<double, 3>> create_turbulence(
    double L,
    int N,
    double k_inj,
    double exponent = -5.0/3.0,
    bool dimensions = false,
    bool divergence_free = false
);

double interpolate_scalar(
    double x,
    double y,
    double z,
    const SampledPositions& source,
    const std::vector<double>& field
);

void interpolate_turbulence(
    const SampledPositions& source,
    SampledPositions& target
);

double calculate_sigma(
    const std::vector<double>& values
);

double outer_smoothing_factor(
    double x,
    double y,
    double z,
    double R_outer,
    double delta_outer
);

#endif