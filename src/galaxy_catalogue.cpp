#include "galaxy_catalogue.hpp"
#include "constants.hpp"
#include "functions.hpp"
#include <cmath>
#include <filesystem>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <H5Cpp.h>
#include <random>

namespace fs = std::filesystem;

void GalaxyCatalogue::ensure_catalogue_exists()
{
    std::string filename = "./catalogue/galaxy_catalogue.hdf5";
    if(fs::exists(filename))
    {
        std::cout << "Found galaxy catalogue: " << filename << "\n";
        return;
    }
    std::cout << "Galaxy catalogue not found. Creating...\n";
    int ret = std::system("python ./catalogue/create_catalogue.py");
    if(ret != 0)
        throw std::runtime_error("Failed to create galaxy catalogue");
    if(!fs::exists(filename))
        throw std::runtime_error("Catalogue creation finished but file missing");
    std::cout << "Created galaxy catalogue successfully\n";
}

void GalaxyCatalogue::load_catalogue()
{
    H5::H5File file("./catalogue/galaxy_catalogue.hdf5", H5F_ACC_RDONLY);

    auto load_1d = [&](std::string name)
    {
        H5::DataSet dataset = file.openDataSet(name);
        H5::DataSpace space = dataset.getSpace();
        hsize_t dims[1];
        space.getSimpleExtentDims(dims);
        std::vector<double> data(dims[0]);
        dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
        return data;
    };

    M200_bins = load_1d("M200");
    redshift_bins = load_1d("redshift");

    auto load_2d = [&](std::string name)
    {
        H5::DataSet dataset = file.openDataSet(name);
        H5::DataSpace space = dataset.getSpace();
        hsize_t dims[2];
        space.getSimpleExtentDims(dims);
        std::vector<std::vector<double>> data(dims[0], std::vector<double>(dims[1]));
        std::vector<double> buffer(dims[0] * dims[1]);
        dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
        for(size_t i = 0; i < dims[0]; i++)
            for(size_t j = 0; j < dims[1]; j++)
                data[i][j] = buffer[i * dims[1] + j];
        return data;
    };

    Rgas_table = load_2d("Rgas");
    Rstar_table = load_2d("Rstar");
    Hgas_table = load_2d("Hgas");
    Hstar_table = load_2d("Hstar");
    concentration_table = load_2d("Concentration");
    Mgas_table = load_2d("Mgas");
    Mstar_table = load_2d("Mstar");
    Mdot_table = load_2d("Mdot");
    Zgas_table = load_2d("Zgas");
}

int GalaxyCatalogue::find_mass_index(double M200)
{
    double logM = std::log10(M200 / Msun_to_g);
    int index = 0;
    double min_diff = std::abs(std::log10(M200_bins[0]) - logM);
    for(size_t i = 1; i < M200_bins.size(); i++)
    {
        double diff = std::abs(std::log10(M200_bins[i]) - logM);
        if(diff < min_diff)
        {
            min_diff = diff;
            index = i;
        }
    }
    return index;
}

int GalaxyCatalogue::find_redshift_index(double z)
{
    int index = 0;
    double min_diff = std::abs(redshift_bins[0] - z);
    for(size_t i = 1; i < redshift_bins.size(); i++)
    {
        double diff = std::abs(redshift_bins[i] - z);
        if(diff < min_diff)
        {
            min_diff = diff;
            index = i;
        }
    }
    return index;
}

GalaxyParameters GalaxyCatalogue::lookup(double M200, double z, const GalaxyVariance& variance)
{
    int im = find_mass_index(M200);
    int iz = find_redshift_index(z);
    GalaxyParameters g;

    g.Rgas = random_variation(Rgas_table[iz][im], variance.Rgas) * kpc_to_cm;
    g.Rstar = random_variation(Rstar_table[iz][im], variance.Rstar) * kpc_to_cm;
    g.Hgas = random_variation(Hgas_table[iz][im], variance.Hgas) * kpc_to_cm;
    g.Hstar = random_variation(Hstar_table[iz][im], variance.Hstar) * kpc_to_cm;
    g.concentration = random_variation(concentration_table[iz][im], variance.concentration);
    g.Mgas = random_variation(Mgas_table[iz][im], variance.Mgas) * Msun_to_g;
    g.Mstellar = random_variation(Mstar_table[iz][im], variance.Mstellar) * Msun_to_g;
    g.Mdot = random_variation(Mdot_table[iz][im], variance.Mdot) * Msun_per_year_to_g_per_s;
    g.Z0 = random_variation(Zgas_table[iz][im], variance.Z0);

    double lower = g.Rgas;
    double upper = variance.Rsonic * g.Rgas;

    if(upper > lower)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(lower, upper);
        g.Rsonic = dist(gen);
    }
    else
    {
        g.Rsonic = lower;
    }

    return g;
}

void GalaxyCatalogue::initialize()
{
    ensure_catalogue_exists();
    load_catalogue();
}