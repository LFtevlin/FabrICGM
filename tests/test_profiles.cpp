#include <iostream>
#include <fstream>
#include <cmath>

#include "functions.hpp"


int main()
{

    // ============================
    // Units: cgs
    // ============================

    double Msun = 1.98847e33;
    double kpc  = 3.085677581e21;


    // ============================
    // Parameters
    // ============================

    double z = 0.0;


    // Halo
    double M200 =
        1e12 * Msun;

    double C =
        10.0;

    double GasFraction =
        0.157;


    // Outer component
    double Rvir =
        200.0 * kpc;


    // Disk
    double Mdisk =
        5e10 * Msun;

    double Rdisk =
        3.0 * kpc;

    double Hdisk =
        0.3 * kpc;



    // ============================
    // Create objects
    // ============================


    HernquistHalo halo(
        M200,
        C,
        GasFraction,
        z
    );


    OuterProfile outer(
        Rvir,
        z
    );


    DoubleExponentialDisk disk(
        Mdisk,
        Rdisk,
        Hdisk
    );



    // ============================
    // Radius grid
    // ============================

    int N = 200;


    double r_min =
        0.01 * kpc;

    double r_max =
        1000.0 * kpc;


    std::ofstream file(
        "profiles.txt"
    );


    file
    << "# r[kpc] "
    << "M_halo[Msun] "
    << "rho_halo[g/cm3] "
    << "phi_halo "
    << "M_outer[Msun] "
    << "rho_outer[g/cm3] "
    << "phi_outer "
    << "M_disk[Msun] "
    << "phi_disk "
    << "\n";



    // ============================
    // Evaluate profiles
    // ============================

    for(int i=0; i<N; i++)
    {

        double frac =
            double(i)/(N-1);


        // logarithmic spacing

        double r =
            r_min *
            pow(
                r_max/r_min,
                frac
            );



        double M_halo =
            halo.MassProfile(r);


        double rho_halo =
            halo.DensityProfile(r);


        double phi_halo =
            halo.PotentialProfile(r);



        double M_outer =
            outer.MassProfile(r);


        double rho_outer =
            outer.DensityProfile(r);


        double phi_outer =
            outer.PotentialProfile(r);



        double M_disk =
            disk.MassProfile(r);


        double phi_disk =
            disk.PotentialProfile(r);



        file
        << r/kpc << " "
        << M_halo/Msun << " "
        << rho_halo << " "
        << phi_halo << " "
        << M_outer/Msun << " "
        << rho_outer << " "
        << phi_outer << " "
        << M_disk/Msun << " "
        << phi_disk
        << "\n";

    }


    file.close();


    std::cout
        << "Written profiles.txt\n";


    return 0;

}