#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// All units are cgs

constexpr double kpc_to_cm = 3.0856775814913673e21;
constexpr double Msun_to_g = 1.98847e33;
constexpr double year_to_s = 3.15576e7;

constexpr double Msun_per_year_to_g_per_s = Msun_to_g / year_to_s;

constexpr double GAMMA = 5.0/3.0;
constexpr double mu = 0.597;
constexpr double m_p = 1.67262192369e-24;
constexpr double k_B = 1.380649e-16;
constexpr double G = 6.67430e-8;

constexpr double H0 = 67.0 * 1e5 / 3.085677581e24; 

constexpr double muG = 1e-6; 

#endif