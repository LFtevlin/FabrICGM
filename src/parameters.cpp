#include "parameters.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <random>
#include "galaxy_catalogue.hpp"
#include "functions.hpp"
#include "constants.hpp"

void Parameters::finalize()
{
    if(galaxy.Mgas <= 0) throw std::runtime_error("Mgas not set");
    if(galaxy.Rgas <= 0) throw std::runtime_error("Rgas not set");
    if(galaxy.Mstellar <= 0) throw std::runtime_error("Mstellar not set");
    if(galaxy.Rstar <= 0) throw std::runtime_error("Rstar not set");
    if(galaxy.Hstar <= 0) throw std::runtime_error("Hstar not set");
    if(galaxy.Hgas <= 0) throw std::runtime_error("Hgas not set");
    if(galaxy.Z0 <= 0) throw std::runtime_error("Z0 not set");
    if(galaxy.concentration <= 0) throw std::runtime_error("C not set");
    if(galaxy.Mdot <= 0) throw std::runtime_error("Mdot not set");
    if(galaxy.Rsonic <= 0) throw std::runtime_error("Rsonic not set");

    if(galaxy.fgas < 0)
    {
        if(M200 > 0 && galaxy.Mgas > 0 && galaxy.Mstellar > 0)
        {
            galaxy.fgas = (galaxy.Mgas + 2 * galaxy.Mstellar) / M200;
        }
        else
        {
            throw std::runtime_error("Cannot compute fgas: missing M200, Mgas, or Mstellar");
        }
    }
}

Parameters read_parameters(const std::string& filename)
{
    Parameters params;
    std::ifstream file(filename);

    if(!file)
        throw std::runtime_error("Could not open parameter file: " + filename);

    std::string line;
    bool found_Rcmax = false;

    while(getline(file, line))
    {
        if(line.empty() || line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string key;
        std::string value;

        getline(ss, key, '=');

        if(!(ss >> value))
            throw std::runtime_error("Could not read value for parameter: " + key);

        if(key == "z")
            params.z = std::stod(value);
        else if(key == "M200")
            params.M200 = std::stod(value) * Msun_to_g;
        else if(key == "runtime")
            params.runtime = std::stod(value);
        else if(key == "tolerance")
            params.tolerance = std::stod(value);
        else if(key == "boxsize")
            params.boxsize = std::stod(value) * kpc_to_cm;
        else if(key == "Z_slope")
            params.Z_slope = std::stod(value);
        else if(key == "Z_profile")
            params.Z_profile = value;
        else if(key == "name")
            params.name = value;
        else if(key == "create_random_galaxy_catalogue")
            params.create_random_galaxy_catalogue = (value == "1" || value == "true");
        else if(key == "magnetic")
            params.magnetic = (value == "1" || value == "true");
        else if(key == "validate_galaxies")
            params.validate_galaxies = (value == "1" || value == "true");
        else if(key == "allow_variance")
            params.allow_variance = (value == "1" || value == "true");

        if(params.allow_variance)
        {
            if(key == "variance_Rgas")
                params.galaxy_variance.Rgas = std::stod(value);
            else if(key == "variance_Rstar")
                params.galaxy_variance.Rstar = std::stod(value);
            else if(key == "variance_Hgas")
                params.galaxy_variance.Hgas = std::stod(value);
            else if(key == "variance_Hstar")
                params.galaxy_variance.Hstar = std::stod(value);
            else if(key == "variance_Mgas")
                params.galaxy_variance.Mgas = std::stod(value);
            else if(key == "variance_Mstar")
                params.galaxy_variance.Mstellar = std::stod(value);
            else if(key == "variance_Z0")
                params.galaxy_variance.Z0 = std::stod(value);
            else if(key == "variance_concentration")
                params.galaxy_variance.concentration = std::stod(value);
            else if(key == "variance_Mdot")
                params.galaxy_variance.Mdot = std::stod(value);
            else if(key == "variance_Rsonic")
                params.galaxy_variance.Rsonic = std::stod(value);
        }

        if(!params.create_random_galaxy_catalogue)
        {
            if(key == "Mstar")
                params.galaxy.Mstellar = random_variation(std::stod(value) * Msun_to_g, params.galaxy_variance.Mstellar);
            else if(key == "Mgas")
                params.galaxy.Mgas = random_variation(std::stod(value) * Msun_to_g, params.galaxy_variance.Mgas);
            else if(key == "Rgas")
                params.galaxy.Rgas = random_variation(std::stod(value) * kpc_to_cm, params.galaxy_variance.Rgas);
            else if(key == "Hgas")
                params.galaxy.Hgas = random_variation(std::stod(value) * kpc_to_cm, params.galaxy_variance.Hgas);
            else if(key == "Rstar")
                params.galaxy.Rstar = random_variation(std::stod(value) * kpc_to_cm, params.galaxy_variance.Rstar);
            else if(key == "Hstar")
                params.galaxy.Hstar = random_variation(std::stod(value) * kpc_to_cm, params.galaxy_variance.Hstar);
            else if(key == "Z0")
                params.galaxy.Z0 = random_variation(std::stod(value), params.galaxy_variance.Z0);
            else if(key == "concentration")
                params.galaxy.concentration = random_variation(std::stod(value), params.galaxy_variance.concentration);
            else if(key == "Mdot")
                params.galaxy.Mdot = random_variation(std::stod(value) * Msun_per_year_to_g_per_s, params.galaxy_variance.Mdot);
            else if(key == "Rsonic")
            {
                double lower = std::stod(value) * kpc_to_cm;
                double upper = params.galaxy_variance.Rsonic * params.galaxy.Rgas;

                if(upper <= lower)
                {
                    params.galaxy.Rsonic = lower;
                }
                else
                {
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_real_distribution<double> dist(lower, upper);
                    params.galaxy.Rsonic = dist(gen);
                }
            }
            else if(key == "fgas")
                params.galaxy.fgas = std::stod(value);
        }

        if(key == "Rcmax")
        {
            params.Rcmax = std::stod(value);
            found_Rcmax = true;
        }

        if(!found_Rcmax)
            params.Rcmax = params.galaxy.Rsonic;

        if(key == "sampling")
            params.sampling = value;
        else if(key == "NBox")
            params.NBox = std::stoi(value);
        else if(key.rfind("LBox", 0) == 0)
        {
            if(params.sampling != "equalmass")
            {
                int index = std::stoi(key.substr(4));

                if(index > 0 && index <= params.NBox)
                {
                    if(params.LBox.size() < static_cast<size_t>(params.NBox))
                        params.LBox.resize(params.NBox);

                    params.LBox[index - 1] = std::stod(value) * kpc_to_cm;
                }
            }
        }
        else if(key.rfind("dx", 0) == 0)
        {
            if(params.sampling != "equalmass")
            {
                int index = std::stoi(key.substr(2));

                if(index > 0 && index <= params.NBox)
                {
                    if(params.dx.size() < static_cast<size_t>(params.NBox))
                        params.dx.resize(params.NBox);

                    params.dx[index - 1] = std::stod(value) * kpc_to_cm;
                }
            }
        }
        else if(key == "turbulence")
            params.turbulence = (value == "true" || value == "1");
        else if(key == "turbulence_scaling")
            params.turbulence_scaling = (value == "true" || value == "1");
        else if(key == "L_inj")
            params.L_inj = std::stod(value) * kpc_to_cm;
        else if(key == "turbulence_rho_percentage")
            params.turbulence_rho_percentage = std::stod(value);
        else if(key == "turbulence_B_percentage")
            params.turbulence_B_percentage = std::stod(value);
        else if(key == "turbulence_v_percentage")
            params.turbulence_v_percentage = std::stod(value);
        else if(key == "Tfloor")
            params.Tfloor = std::stod(value);
        else if(key == "Nneg")
            params.Nneg = std::stod(value);
        else if(key == "M_cutoff")
            params.m_cutoff = std::stod(value);

        if(key == "mtarget")
        {
            if(params.sampling == "equalmass" || params.sampling == "both")
                params.mtarget = std::stod(value) * Msun_to_g;
        }
    }

    if(params.M200 < 0)
        throw std::runtime_error("Missing parameter: M200");

    if(params.z < 0)
        throw std::runtime_error("Missing parameter: z");

    return params;
}