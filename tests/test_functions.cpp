#include "functions.hpp"

#include <iostream>
#include <vector>
#include <cmath>


int main()
{

    // =====================================================
    // Load cooling table
    // =====================================================

    CoolingTable cooling(
        "./UVB_dust1_CR1_G1_shield0.hdf5"
    );


    std::cout << "\n=== Basic physics ===\n";


    // cgs inputs
    double M =
        1e12 * 1.98847e33;      // solar masses -> g

    double r =
        100.0 * 3.08567758e21;  // 100 kpc -> cm

    double rho =
        1e-27;                 // g/cm3

    double T =
        1e6;                   // K

    double v =
        -1e7;                  // cm/s

    double Z =
        1.0;

    double z =
        0.0;


    std::cout
        << "vc [cm/s] = "
        << v_c(M,r)
        << "\n";


    std::cout
        << "Tc [K] = "
        << T_c(M,r)
        << "\n";


    std::cout
        << "cs [cm/s] = "
        << c_s(T)
        << "\n";


    std::cout
        << "Mach = "
        << Mach(T,v)
        << "\n";


    std::cout
        << "t_flow [s] = "
        << t_flow(r,v)
        << "\n";


    std::cout
        << "t_ff [s] = "
        << t_ff(r,M)
        << "\n";


    std::cout << "before t_cool\n";

    std::cout << "before Lambda\n";

    double lam =
        cooling.Lambda(
            rho,
            T,
            Z,
            z
        );

    std::cout << "Lambda = "
            << lam
            << "\n";

    std::cout << "after Lambda\n";

    double tc =
        t_cool(
            T,
            rho,
            Z,
            z,
            cooling
        );

    std::cout << "after t_cool\n";

    std::cout
        << "t_cool [s] = "
        << tc
        << "\n";


    double pot =
        -G*M/r;


    std::cout
        << "Bernoulli [cm2/s2] = "
        << Bernoulli(v,T,pot)
        << "\n";


    double P =
        rho/(mu*m_p)*k_B*T;


    std::cout
        << "Entropy = "
        << entropy(P,rho)
        << "\n";


    std::cout
        << "Bound? "
        << check_Bernoulli(v,T,pot)
        << "\n";


    std::cout
        << "Mach OK? "
        << check_Mach(v,T)
        << "\n";



    // =====================================================
    // Self-similar tests
    // =====================================================

    std::cout
        << "\n=== Self similar ===\n";


    double Rvir =
        200.0 * 3.08567758e21;


    double Mdot =
        1.0e25;       // g/s


    double A = 2.0;
    double B = 1.1;
    double m = 0.5;
    double X = 0.76;


    std::cout
        << "vc_ss = "
        << v_c_ss(
            Rvir,
            M,
            r,
            m
        )
        << "\n";


    std::cout
        << "T_ss = "
        << T_ss(
            M,
            r,
            A
        )
        << "\n";


    std::cout
        << "K_ss = "
        << K_ss(
            r,
            B
        )
        << "\n";


    double Lambda =
        cooling.Lambda(
            rho,
            T,
            Z,
            z
        );


    std::cout
        << "Mach_ss = "
        << Mach_ss(
            Rvir,
            M,
            r,
            m,
            A,
            B,
            Mdot,
            Lambda,
            X
        )
        << "\n";


    std::cout
        << "nH_ss = "
        << n_H_ss(
            Rvir,
            M,
            r,
            Mdot,
            B,
            A,
            Lambda,
            m
        )
        << "\n";



    // =====================================================
    // RK test setup
    // =====================================================

    std::cout
        << "\n=== RK4 ===\n";


    int N = 100;


    std::vector<double> radius(N);
    std::vector<double> Mcum(N);
    std::vector<double> Zarr(N);


    for(int i=0;i<N;i++)
    {
        radius[i] =
            (10.0+i)
            *
            3.08567758e21;

        Mcum[i] =
            M;

        Zarr[i] =
            1.0;
    }


    double dr =
        radius[10]*0.01;


    auto result =
        RK4_step(
            radius[10],
            T,
            rho,
            v,
            Mcum,
            radius,
            dr,
            Mdot,
            Zarr,
            z,
            cooling
        );


    std::cout
        << "RK4 step v = "
        << result.first
        << "\n";


    std::cout
        << "RK4 step T = "
        << result.second
        << "\n";



    // =====================================================
    // find_Tinit test
    // =====================================================

    std::cout
        << "\n=== find_Tinit ===\n";


    auto init =
        find_Tinit(
            radius,
            Mcum,
            T,
            rho,
            v,
            Mdot,
            0.5,
            radius[10],
            10,
            1e-2,
            Zarr,
            z,
            cooling
        );


    std::cout
        << "Tinit = "
        << std::get<0>(init)
        << "\n";


    std::cout
        << "vinit = "
        << std::get<1>(init)
        << "\n";


    std::cout
        << "rhoinit = "
        << std::get<2>(init)
        << "\n";


    std::cout
        << "\nFinished successfully\n";


    return 0;
}