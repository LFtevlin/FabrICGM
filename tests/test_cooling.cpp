#include "functions.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


int main()
{

    CoolingTable table(
        "./UVB_dust1_CR1_G1_shield0.hdf5"
    );


    std::ofstream outfile("cooling_test.txt");

    double rho;
    double T;
    double Z;
    double z;


    if(!outfile)
    {
        std::cerr << "Could not open output file\n";
        return 1;
    }


    /*
       Test cases:

       vary:
       - temperature
       - density
       - metallicity
       - redshift

       Output:
       x-axis value
       Lambda
       label parameter
    */





    // -------------------------
    // 1) Temperature scan
    // -------------------------

    rho = 1e-27;     // g/cm3
    Z   = 1.0;       // solar metallicity
    z   = 0.0;


    for(int i=0;i<100;i++)
    {

        double T =
            pow(10.0, 1.5 + i*(7.0-1.5)/99.0);


        double lambda =
            table.Lambda(
                rho,
                T,
                Z,
                z
            );


        outfile
            << "T "
            << T
            << " "
            << lambda
            << "\n";
    }



    // -------------------------
    // 2) Density scan
    // -------------------------

    T = 1e6;
    Z = 1.0;
    z = 0.0;


    for(int i=0;i<100;i++)
    {

        double nH =
            pow(10.0,-6.0 + i*(4.0+6.0)/99.0);


        // inverse of:
        // nH = rho/(mu*m_p)*0.76

        double rho =
            nH*0.597*1.6726219e-24/0.76;


        double lambda =
            table.Lambda(
                rho,
                T,
                Z,
                z
            );


        outfile
            << "nH "
            << nH
            << " "
            << lambda
            << "\n";
    }



    // -------------------------
    // 3) Metallicity scan
    // -------------------------

    rho =
        1e-27;

    T =
        1e6;

    z =
        0.0;


    for(int i=0;i<100;i++)
    {

        double Z =
            pow(10.0,-3.0+i*3.0/99.0);


        double lambda =
            table.Lambda(
                rho,
                T,
                Z,
                z
            );


        outfile
            << "Z "
            << Z
            << " "
            << lambda
            << "\n";

    }



    // -------------------------
    // 4) Redshift scan
    // -------------------------

    rho = 1e-27;
    T = 1e6;
    Z = 1.0;


    for(int i=0;i<100;i++)
    {

        double z =
            i*10.0/99.0;


        double lambda =
            table.Lambda(
                rho,
                T,
                Z,
                z
            );


        outfile
            << "z "
            << z
            << " "
            << lambda
            << "\n";

    }


    outfile.close();


    std::cout
        << "Finished. Written cooling_test.txt\n";


    return 0;
}