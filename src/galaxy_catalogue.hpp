#pragma once

#include "parameters.hpp"
#include <vector>

class GalaxyCatalogue
{
private:
    std::vector<double> M200_bins;
    std::vector<double> redshift_bins;
    std::vector<std::vector<double>> Rgas_table;
    std::vector<std::vector<double>> Rstar_table;
    std::vector<std::vector<double>> Hgas_table;
    std::vector<std::vector<double>> Hstar_table;
    std::vector<std::vector<double>> concentration_table;
    std::vector<std::vector<double>> Mgas_table;
    std::vector<std::vector<double>> Mstar_table;
    std::vector<std::vector<double>> Mdot_table;
    std::vector<std::vector<double>> Zgas_table;
    int find_mass_index(double M200);
    int find_redshift_index(double z);

public:
    GalaxyCatalogue() = default;
    void ensure_catalogue_exists();
    void load_catalogue();
    void initialize();
    GalaxyParameters lookup(double M200, double z, const GalaxyVariance& variance);
};