#include <string>
#include <vector>

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>

struct GalaxyParameters
{
    double Mgas = -1.;
    double Mstellar = -1.;
    double fgas = -1.;

    double Rgas = -1.;
    double Rstar = -1.;

    double Hgas = -1.;
    double Hstar = -1.;

    double concentration = -1.;
    double Z0 = -1.;

    double Mdot = -1.;
    double Rsonic = -1.;
};

struct GalaxyVariance
{
    double Rgas = 1.;
    double Rstar = 1.;
    double Hgas = 1.;
    double Hstar = 1.;

    double Mgas = 1.;
    double Mstellar = 1.;

    double Z0 = 1.;

    double concentration = 1.;

    double Mdot = 1.0;

    double Rsonic = 0. ;
};


struct Parameters
{

    std::string name;

    bool create_random_galaxy_catalogue = false;
    bool validate_galaxies = false;
    bool allow_variance = false;

    GalaxyParameters galaxy;
    GalaxyVariance galaxy_variance;

    //
    // numerical options
    //
    double boxsize = -1.0;


    //
    // halo
    //
    double M200 = -1.0;
    double concentration = -1.0;
    double z = -1.0;

    double Rcmax = -1.0;
    std::string sampling = "equalmass";

    int NBox = 0;

    std::vector<double> LBox;
    std::vector<double> dx;



    double BulgeFraction = 0.0;
    double BulgeScaleLength = -1.0;


    double Z_slope = 0.0;
    std::string Z_profile = "constant";

    double runtime;
    double tolerance;

    double mtarget = -1.0;

    bool magnetic = false;

    bool turbulence = false;
    bool turbulence_scaling = false;

    double L_inj = 0.0;

    double turbulence_rho_percentage = 0.0;
    double turbulence_B_percentage   = 0.0;
    double turbulence_v_percentage   = 0.0;

    double Tfloor = 1e2;
    double Nneg = 500;
    double m_cutoff=4;



    void finalize();

};


Parameters read_parameters(
    const std::string& filename
);

#endif