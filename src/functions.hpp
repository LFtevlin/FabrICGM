#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>
#include <string>
#include <H5Cpp.h>
#include <tuple>
#include <utility>
#include <functional>
#include "constants.hpp"
#include "cooling.hpp"

double random_variation(double value, double factor);

double v_c(double M, double r);

double T_c(double M, double r);

double t_flow(double r, double v);

double c_s(double T);

double Mach(double T, double v);

double t_cool(double T, double rho, double Z_CGM, double z, CoolingTable& cooling);

double t_ff(double r, double M);

double Bernoulli(double v, double T, double pot);

double entropy(double P, double rho);

bool check_Bernoulli(double v, double T, double pot);

bool check_Mach(double v, double T);

double dlnrho_dlnr(double T, double rho, double v, double M, double r, double Z_CGM, double z, CoolingTable& cooling);

double dlnT_dlnr(double T, double rho, double v, double M, double r, double Z_CGM, double z, CoolingTable& cooling);

double dlnv_dlnr(double T, double rho, double v, double M, double r, double Z_CGM, double z, CoolingTable& cooling);

double interp_linear(const std::vector<double>& x, const std::vector<double>& y, double x0);

void derivatives(double T, double rho, double v, double M, double r, double Z, double z, CoolingTable& cooling, double& dvdr, double& dTdr);

double compute_dMCGM(double r_curr, double dr, double rho);

std::pair<double,double> RK4_step(double r_curr, double T, double rho, double v, const std::vector<double>& M, const std::vector<double>& r, double dr, double M_dot, const std::vector<double>& Z_CGM, double z, CoolingTable& cooling);

std::pair<double,double> RK4(double r_curr, double T, double rho, double v, const std::vector<double>& M, const std::vector<double>& r, double dr_total, double M_dot, const std::vector<double>& Z_CGM, double z, CoolingTable& cooling, double tol = 1e-1, int max_iter = 30);

// double v_c_ss(double R_vir, double M_vir, double r, double m);

// double T_ss(double M, double r, double A);

// double K_ss(double r, double B);

// double Mach_ss(double R_vir, double M_vir, double r, double m, double A, double B, double M_dot, double Lambdav, double X);

// double n_H_ss(double R_vir, double M_vir, double r, double M_dot, double B, double A, double Lambdav, double m);

// double tc_tf_ss(double B);

// double tc_tff_ss(double R_vir, double M_vir, double r, double m, double A, double B, double M_dot, double Lambdav, double X);

void new_setup(double v_init, double rho_init, double T_new, int index, const std::vector<double>& radius, std::vector<double>& velocities, std::vector<double>& temperatures, std::vector<double>& rho);

double condition(double r, double M, double T, double rho, double v, double M_dot, double z, double Z_CGM, CoolingTable& cooling);

std::tuple<double,double,double> find_Tinit(const std::vector<double>& radius, const std::vector<double>& M_cum, double T_init, double rho_init, double v_init, double M_dot, double Mach_init, double R_sonic, int index, double epsilon, const std::vector<double>& Z_CGM, double z, CoolingTable& cooling, int n_iter = 40);

class HernquistHalo
{
public:
    HernquistHalo(double M200, double C, double GasFraction, double z);
    double DensityProfile(double r) const;
    double PotentialProfile(double r) const;
    double AccelerationProfile(double r) const;
    double MassProfile(double r) const;
    double R200;
    double a;
    double M200;
    double C;
    double RhoCrit;
    double Overdensity;
private:
    double GasFraction;
};

class HernquistBulge
{
public:
    HernquistBulge(double Mass, double a);
    double DensityProfile(double r) const;
    double PotentialProfile(double r) const;
    double AccelerationProfile(double r) const;
    double MassProfile(double r) const;
private:
    double Mass_bulge;
    double a_bulge;
};

class OuterProfile
{
public:
    OuterProfile(double R, double z);
    double DensityProfile(double r) const;
    double MassProfile(double r) const;
    double AccelerationProfile(double r) const;
    double PotentialProfile(double r, double factor2 = 2.0) const;
private:
    double rho_mean;
    double R;
    double factor;
};

class DoubleExponentialDisk
{
public:
    DoubleExponentialDisk(double Mass, double DiskScaleRadius, double DiskScaleHeight);
    double MassProfile(double r) const;
    double AccelerationProfile(double r) const;
    double PotentialProfile(double r) const;
private:
    double Mass_disk;
    double L;
    double H;
};

double Gamma0(double x);

double GammaMinus1(double x);

double spherical_theta(double r, double z);

double interpolate(double x, const std::vector<double>& x_array, const std::vector<double>& y_array);

double powerlaw_transition(double r, double r_min, double r_max, double y_min, double y_max);

void rho_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& rho, const std::vector<double>& Z, std::vector<std::vector<double>>& P_grid, double R_cmax);

void P_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& P, const std::vector<double>& Z, std::vector<std::vector<double>>& P_grid, double R_cmax);

void vr_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& velocity, const std::vector<double>& Z, std::vector<std::vector<double>>& vr_grid, double R_cmax);

void vt_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& velocity, const std::vector<double>& Z, std::vector<std::vector<double>>& vt_grid, double R_cmax);

void vp_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& M, const std::vector<double>& Z, std::vector<std::vector<double>>& vp_grid, double R_cmax);

double interpolate_spherical_grid(double r, double theta, const std::vector<std::vector<double>>& r_grid, const std::vector<std::vector<double>>& theta_grid, const std::vector<std::vector<double>>& field);

double interpolate_cylindrical_grid(double R, double Z, const std::vector<double>& R_grid, const std::vector<double>& Z_grid, const std::vector<std::vector<double>>& field);

double ellipsoid_smoothing_factor(double x, double y, double z, double R_scale, double H_scale, double m_cutoff);

double estimate_max_turbulence_percentage(size_t N_cells, double max_negative_cells);

#endif