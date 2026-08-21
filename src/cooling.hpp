#ifndef COOLING_HPP
#define COOLING_HPP

#include <vector>
#include <string>
#include <H5Cpp.h>
#include <tuple>
#include <utility>
#include <functional>
#include "constants.hpp"



// #pragma once

// class CoolingTable
// {
// public:

//     CoolingTable(
//         const std::string& filename
//     );

//     double Lambda(
//         double rho,
//         double T,
//         double Z,
//         double z
//     );

// private:

//     int find_nearest(
//         const std::vector<double>& array,
//         double value
//     );

//     std::vector<double> nH_tab;
//     // std::vector<double> T_tab;
//     std::vector<double> Z_tab;
//     std::vector<double> U_tab;
//     std::vector<double> z_tab;

//     std::vector<double> cooling;
//     std::vector<double> heating;

//     std::vector<hsize_t> cooling_shape;
//     std::vector<hsize_t> heating_shape;
// };

// #endif

#pragma once



class CoolingTable
{

public:

    CoolingTable(const std::string& filename);


    double Lambda(double rho,double T,double Z,double z);


private:

    int find_nearest(const std::vector<double>& array,double value);


    std::vector<double> nH_tab;
    std::vector<double> Z_tab;
    std::vector<double> z_tab;
    std::vector<double> T_tab;


    std::vector<double> cooling;
    std::vector<double> heating;


    std::vector<hsize_t> cooling_shape;
    std::vector<hsize_t> heating_shape;


    int cool_TotalPrim = -1;
    int cool_TotalMetal = -1;

    int heat_TotalPrim = -1;
    int heat_TotalMetal = -1;

    std::vector<std::string> cool_ids;
    std::vector<std::string> heat_ids;

    int n_species_cool;
    int n_species_heat;

};


#endif