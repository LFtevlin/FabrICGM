#ifndef SAMPLED_POSITIONS_HPP
#define SAMPLED_POSITIONS_HPP

#include <vector>
#include <array>

struct SampledPositions
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    std::vector<double> rho;
    std::vector<double> P;
    std::vector<double> T;
    std::vector<double> V;
    std::vector<double> Z;

    std::vector<std::array<double,3>> v;
    std::vector<std::array<double,3>> B;


};

#endif