#include "functions.hpp"
#include "constants.hpp"

#include <H5Cpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <tuple>




double interpolate(double x, const std::vector<double>& xs, const std::vector<double>& ys)
{
    if(x < xs.front() || x > xs.back()) return 0.0;
    auto it = std::lower_bound(xs.begin(), xs.end(), x);
    int i = std::distance(xs.begin(), it);
    if(i == 0) return ys[0];
    double x1 = xs[i-1];
    double x2 = xs[i];
    double y1 = ys[i-1];
    double y2 = ys[i];
    double w = (x - x1) / (x2 - x1);
    return y1 + w * (y2 - y1);
}

// void cartesian_to_spherical(
//     double x,
//     double y,
//     double z,
//     double& r,
//     double& theta,
//     double& phi
// )
// {
//     r = std::sqrt(
//         x*x + y*y + z*z
//     );

//     if(r > 0)
//     {
//         phi = std::atan2(y, x);
//     }
//     else
//     {
//         theta = 0.0;
//     }
//     phi = std::atan(y/x);
// }

// void cartesian_to_cylindrical(
//     double x,
//     double y,
//     double z,
//     double& R,
//     double& Z
// )
// {
//     R = std::sqrt(
//         x*x + y*y
//     );

//     Z = z;
// }

double interpolate_spherical_grid(double r, double theta, const std::vector<std::vector<double>>& r_grid, const std::vector<std::vector<double>>& theta_grid, const std::vector<std::vector<double>>& field)
{
    int NR = r_grid.size();
    int NZ = r_grid[0].size();
    int i = 0;
    double dr_min = std::abs(r_grid[0][0] - r);
    for(int ii = 1; ii < NR; ii++)
    {
        double dr = std::abs(r_grid[ii][0] - r);
        if(dr < dr_min)
        {
            dr_min = dr;
            i = ii;
        }
    }
    int j = 0;
    double dtheta_min = std::abs(theta_grid[i][0] - theta);
    for(int jj = 1; jj < NZ; jj++)
    {
        double dtheta = std::abs(theta_grid[i][jj] - theta);
        if(dtheta < dtheta_min)
        {
            dtheta_min = dtheta;
            j = jj;
        }
    }
    return field[i][j];
}

double interpolate_cylindrical_grid(double R, double Z, const std::vector<double>& R_grid, const std::vector<double>& Z_grid, const std::vector<std::vector<double>>& field)
{
    int NR = R_grid.size();
    int NZ = Z_grid.size();
    int i = 0;
    double dR_min = std::abs(R_grid[0] - R);
    for(int ii = 1; ii < NR; ii++)
    {
        double dR = std::abs(R_grid[ii] - R);
        if(dR < dR_min)
        {
            dR_min = dR;
            i = ii;
        }
    }
    int j = 0;
    double dZ_min = std::abs(Z_grid[0] - Z);
    for(int jj = 1; jj < NZ; jj++)
    {
        double dZ = std::abs(Z_grid[jj] - Z);
        if(dZ < dZ_min)
        {
            dZ_min = dZ;
            j = jj;
        }
    }
    return field[i][j];
}




double random_variation(
    double value,
    double factor
)
{
    if (factor <= 1.0)
        return value;


    static std::mt19937 rng(
        std::random_device{}()
    );


    std::uniform_real_distribution<double> dist(-1.0, 1.0);


    double x = dist(rng);


    return value * std::pow(factor, x);
}


double v_c(double M, double r)
{
    return sqrt(G * M / r); // cm/s
}

double T_c(double M, double r)
{
    double vc = v_c(M, r);
    return mu * m_p * vc * vc / (GAMMA * k_B); // K
}

double t_flow(double r, double v)
{
    return r / fabs(v); // s
}

double c_s(double T)
{
    return sqrt(GAMMA * k_B * T / (mu * m_p)); // cm/s
}

double Mach(double T, double v)
{
    return fabs(v) / c_s(T);
}

double t_cool(double T, double rho, double Z_CGM, double z, CoolingTable& cooling)
{
    double n = rho / (mu * m_p);
    double P = n * k_B * T;
    double nH = 0.76 * rho / (mu * m_p);
    double Lambda = cooling.Lambda(rho, T, Z_CGM, z);
    return P / ((GAMMA - 1.0) * nH * nH * Lambda);
}

double t_ff(double r, double M)
{
    return sqrt(2.0) * r / v_c(M, r);
}

double Bernoulli(double v, double T, double pot)
{
    double cs = c_s(T);
    return 0.5 * v * v + cs * cs / (GAMMA - 1.0) + pot;
}

double entropy(double P, double rho)
{
    double n = rho / (mu * m_p);
    return P / std::pow(n, GAMMA);
}

bool check_Bernoulli(double v, double T, double pot)
{
    return Bernoulli(v, T, pot) < 0.0;
}

bool check_Mach(double v, double T)
{
    return Mach(T, v) < 1 - 1e-10;
}

double dlnrho_dlnr(double T, double rho, double v, double M, double r, double Z_CGM, double z, CoolingTable& cooling)
{
    return -2.0 - dlnv_dlnr(T, rho, v, M, r, Z_CGM, z, cooling);
}

double dlnT_dlnr(double T, double rho, double v, double M, double r, double Z_CGM, double z, CoolingTable& cooling)
{
    double tcool = t_cool(T, rho, Z_CGM, z, cooling);
    double tflow = t_flow(r, v);
    double mach = fabs(v) / c_s(T);
    double denom = 1.0 - mach * mach;
    return t_flow(r, v) / t_cool(T, rho, Z_CGM, z, cooling) + (GAMMA - 1.0) * dlnrho_dlnr(T, rho, v, M, r, Z_CGM, z, cooling);
}

double dlnv_dlnr(double T, double rho, double v, double M, double r, double Z_CGM, double z, CoolingTable& cooling)
{
    double tf = std::abs(t_flow(r, v));
    double mach = Mach(T, v);
    double cs = c_s(T);
    double vc = v_c(M, r);
    double tcool = t_cool(T, rho, Z_CGM, z, cooling);
    double denom = mach * mach - 1.0;
    double numerator = 2.0 - vc * vc / (cs * cs) - tf / (GAMMA * tcool);
    const double sonic_tol = 1e-8;
    if(std::abs(denom) < sonic_tol)
    {
        if(std::abs(numerator) < sonic_tol)
        {
            return 0.0;
        }
        throw std::runtime_error("Non-critical sonic point.");
    }
    return numerator / denom;
}

double interp_linear(const std::vector<double>& x, const std::vector<double>& y, double x0)
{
    if(x0 <= x.front()) return y.front();
    if(x0 >= x.back()) return y.back();
    for(size_t i = 0; i < x.size() - 1; i++)
    {
        if(x0 >= x[i] && x0 <= x[i + 1])
        {
            double f = (x0 - x[i]) / (x[i + 1] - x[i]);
            return y[i] + f * (y[i + 1] - y[i]);
        }
    }
    return 0.0;
}

void derivatives(double T, double rho, double v, double M, double r, double Z, double z, CoolingTable& cooling, double& dvdr, double& dTdr)
{
    double dlnv = dlnv_dlnr(T, rho, v, M, r, Z, z, cooling);
    double dlnT = dlnT_dlnr(T, rho, v, M, r, Z, z, cooling);
    dvdr = dlnv * v / r;
    dTdr = dlnT * T / r;
}

double compute_dMCGM(double r_curr, double dr, double rho)
{
    double integrand1 = 4.0 * M_PI * r_curr * r_curr * dr * rho;
    return integrand1;
}


std::pair<double,double> RK4_step(double r_curr, double T, double rho, double v, const std::vector<double>& M, const std::vector<double>& r, double dr, double M_dot, const std::vector<double>& Z_CGM, double z, CoolingTable& cooling)
{
    double M_curr = interp_linear(r, M, r_curr);
    double Z_curr = interp_linear(r, Z_CGM, r_curr);
    double k1v, k1T;
    derivatives(T, rho, v, M_curr, r_curr, Z_curr, z, cooling, k1v, k1T);
    double r_mid = r_curr + 0.5 * dr;
    double T_mid = T + 0.5 * k1T * dr;
    if(T_mid <= 0 || !std::isfinite(T_mid))
    {
        throw std::runtime_error("RK4 midpoint temperature invalid");
    }
    double v_mid = v + 0.5 * k1v * dr;
    double rho_mid = M_dot / (4.0 * M_PI * r_mid * r_mid * std::abs(v_mid));
    double M_CGMmid = compute_dMCGM(r_curr, 0.5 * dr, rho_mid);
    double M_mid = interp_linear(r, M, r_mid) + M_CGMmid;
    double Z_mid = interp_linear(r, Z_CGM, r_mid);
    double k2v, k2T;
    derivatives(T_mid, rho_mid, v_mid, M_mid, r_mid, Z_mid, z, cooling, k2v, k2T);
    T_mid = T + 0.5 * k2T * dr;
    if(T_mid <= 0 || !std::isfinite(T_mid))
    {
        throw std::runtime_error("RK4 midpoint temperature invalid");
    }
    v_mid = v + 0.5 * k2v * dr;
    rho_mid = M_dot / (4.0 * M_PI * r_mid * r_mid * std::abs(v_mid));
    double k3v, k3T;
    derivatives(T_mid, rho_mid, v_mid, M_mid, r_mid, Z_mid, z, cooling, k3v, k3T);
    double r_end = r_curr + dr;
    double T_end = T + k3T * dr;
    if(T_end <= 0 || !std::isfinite(T_end))
    {
        throw std::runtime_error("RK4 midpoint temperature invalid");
    }
    double v_end = v + k3v * dr;
    double rho_end = M_dot / (4.0 * M_PI * r_end * r_end * std::abs(v_end));
    double M_CGMend = compute_dMCGM(r_curr, dr, rho_end);
    double M_end = interp_linear(r, M, r_end) + M_CGMend;
    double Z_end = interp_linear(r, Z_CGM, r_end);
    double k4v, k4T;
    derivatives(T_end, rho_end, v_end, M_end, r_end, Z_end, z, cooling, k4v, k4T);
    double vn = v + dr * (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0;
    double Tn = T + dr * (k1T + 2 * k2T + 2 * k3T + k4T) / 6.0;
    return {vn, Tn};
}


std::pair<double,double> RK4(double r_curr, double T, double rho, double v, const std::vector<double>& M, const std::vector<double>& r, double dr_total, double M_dot, const std::vector<double>& Z_CGM, double z, CoolingTable& cooling, double tol, int max_iter)
{
    double dr_done = 0.0;
    double T_curr = T;
    double v_curr = v;
    while(dr_done < dr_total)
    {
        double dr_try = dr_total - dr_done;
        int iter = 0;
        while(true)
        {
            auto result = RK4_step(r_curr, T_curr, rho, v_curr, M, r, dr_try, M_dot, Z_CGM, z, cooling);
            double vn = result.first;
            double Tn = result.second;
            if(!std::isfinite(vn) || !std::isfinite(Tn) || Tn <= 0.0)
            {
                dr_try *= 0.5;
                iter++;
                if(iter > max_iter) throw std::runtime_error("RK4 failed");
                continue;
            }
            double dv_rel = std::abs((vn - v_curr) / v_curr);
            double dT_rel = std::abs((Tn - T_curr) / T_curr);
            if(dv_rel > tol || dT_rel > tol)
            {
                dr_try *= 0.5;
                iter++;
                if(iter > max_iter) throw std::runtime_error("Adaptive RK4 failed");
                continue;
            }
            v_curr = vn;
            T_curr = Tn;
            r_curr += dr_try;
            dr_done += dr_try;
            rho = M_dot / (4.0 * M_PI * r_curr * r_curr * std::abs(v_curr));
            break;
        }
    }
    return {v_curr, T_curr};
}



// double v_c_ss(double R_vir, double M_vir, double r, double m)
// {
//     return v_c(M_vir, R_vir) * std::pow(r / R_vir, m);
// }

// double T_ss(double M, double r, double A)
// {
//     return T_c(M, r) / A;
// }

// double K_ss(double r, double B)
// {
//     return std::pow(r, B);
// }

// double Mach_ss(double R_vir, double M_vir, double r, double m, double A, double B, double M_dot, double Lambdav, double X)
// {
//     double vc = v_c_ss(R_vir, M_vir, r, m);
//     return X * A / (m_p * vc * vc) * std::sqrt(5.0 * M_dot * Lambdav / (18.0 * M_PI * B * r));
// }

// double n_H_ss(double R_vir, double M_vir, double r, double M_dot, double B, double A, double Lambdav, double m)
// {
//     double vc = v_c_ss(R_vir, M_vir, r, m);
//     return std::sqrt(9.0 * B * M_dot / (40.0 * M_PI * A * r * r * r * Lambdav)) * vc;
// }

// double tc_tf_ss(double B)
// {
//     return 1.0 / B;
// }

// double tc_tff_ss(double R_vir, double M_vir, double r, double m, double A, double B, double M_dot, double Lambdav, double X)
// {
//     return std::sqrt(A / 2.0) / B / Mach_ss(R_vir, M_vir, r, m, A, B, M_dot, Lambdav, X);
// }

void new_setup(double v_init, double rho_init, double T_new, int index, const std::vector<double>& radius, std::vector<double>& velocities, std::vector<double>& temperatures, std::vector<double>& rho)
{
    int N = radius.size();
    velocities.assign(N, 0.0);
    temperatures.assign(N, 0.0);
    rho.assign(N, 0.0);
    velocities[index] = v_init;
    temperatures[index] = T_new;
    rho[index] = rho_init;
}




double condition(double r, double M, double T, double rho, double v, double M_dot, double z, double Z_CGM, CoolingTable& cooling)
{
    return 2.0 - std::pow(v_c(M, r), 2) / std::pow(c_s(T), 2) - t_flow(r, v) / (GAMMA * t_cool(T, rho, Z_CGM, z, cooling));
}



std::tuple<double,double,double> find_Tinit(const std::vector<double>& radius, const std::vector<double>& M_cum, double T_init, double rho_init, double v_init, double M_dot, double Mach_init, double R_sonic, int index, double epsilon, const std::vector<double>& Z_CGM, double z, CoolingTable& cooling, int n_iter)
{
    int count = 0;
    auto update_state = [&]()
    {
        double cs_init = c_s(T_init);
        v_init = Mach_init * cs_init;
        rho_init = M_dot / (4.0 * M_PI * R_sonic * R_sonic * v_init);
    };
    double cond = condition(radius[index], M_cum[index], T_init, rho_init, v_init, M_dot, z, Z_CGM[index], cooling);
    if(epsilon > 0)
    {
        while((cond > epsilon) || (cond < 0.0))
        {
            if(cond < 0.0) T_init *= 1.05;
            else if(cond > epsilon) T_init *= 0.95;
            update_state();
            cond = condition(radius[index], M_cum[index], T_init, rho_init, v_init, M_dot, z, Z_CGM[index], cooling);
            count++;
            if(count > n_iter) break;
        }
    }
    if(epsilon < 0)
    {
        while((cond < epsilon) || (cond > 0.0))
        {
            if(cond > 0.0) T_init *= 0.95;
            else if(cond < epsilon) T_init *= 1.05;
            update_state();
            cond = condition(radius[index], M_cum[index], T_init, rho_init, v_init, M_dot, z, Z_CGM[index], cooling);
            count++;
            if(count > n_iter) break;
        }
    }
    return {T_init, v_init, rho_init};
}

HernquistHalo::HernquistHalo(double M200_in, double C_in, double GasFraction_in, double z)
{
    double Omega_m = 0.3;
    double Omega_l = 0.7;
    double a_scale = 1.0 / (1.0 + z);
    double H = H0 * sqrt(Omega_m / pow(a_scale, 3) + Omega_l);
    RhoCrit = 3.0 * H * H / (8.0 * M_PI * G);
    M200 = M200_in;
    C = C_in;
    GasFraction = GasFraction_in;
    Overdensity = (200.0 / 3.0) * pow(C, 3) / (log(1 + C) - C / (1 + C));
    R200 = cbrt(3.0 * M200 / (4.0 * M_PI * 200.0 * RhoCrit));
    a = cbrt(G * M200 / (100.0 * H * H)) / C * sqrt(2.0 * (log(1 + C) - C / (1 + C)));
}

double HernquistHalo::DensityProfile(double r) const
{
    return M200 * (1.0 - GasFraction) / (2.0 * M_PI * pow(a, 3)) / (r / a) / pow(1.0 + r / a, 3);
}

double HernquistHalo::PotentialProfile(double r) const
{
    return -G * M200 * (1.0 - GasFraction) / (a + r);
}

double HernquistHalo::AccelerationProfile(double r) const
{
    return -G * M200 * (1.0 - GasFraction) / pow(a + r, 2);
}

double HernquistHalo::MassProfile(double r) const
{
    return M200 * (1.0 - GasFraction) * pow(r, 2) / pow(r + a, 2);
}

HernquistBulge::HernquistBulge(double Mass, double a)
    : Mass_bulge(Mass), a_bulge(a)
{
}

double HernquistBulge::DensityProfile(double r) const
{
    return Mass_bulge / (2.0 * M_PI * pow(a_bulge, 3)) / (r / a_bulge) / pow(1.0 + r / a_bulge, 3);
}

double HernquistBulge::PotentialProfile(double r) const
{
    return -G * Mass_bulge / (a_bulge + r);
}

double HernquistBulge::AccelerationProfile(double r) const
{
    return -G * Mass_bulge / pow(a_bulge + r, 2);
}

double HernquistBulge::MassProfile(double r) const
{
    return Mass_bulge * pow(r, 2) / pow(r + a_bulge, 2);
}

OuterProfile::OuterProfile(double R_in, double z)
{
    double Omega_m = 0.3;
    double Omega_l = 0.7;
    double a = 1.0 / (1.0 + z);
    double H = H0 * sqrt(Omega_m / pow(a, 3) + Omega_l);
    double rho_cr = 3.0 * H * H / (8.0 * M_PI * G);
    rho_mean = 0.31 * rho_cr;
    R = R_in;
    factor = 5.0;
}

double OuterProfile::DensityProfile(double r) const
{
    return rho_mean * (pow(r / (R * factor), -1.5) + 1.0);
}

double OuterProfile::MassProfile(double r) const
{
    return (4.0 * M_PI * rho_mean * (pow(r, 3) + 2.0 * pow(R, 1.5) * pow(factor, 1.5) * pow(r, 1.5))) / 3.0;
}

double OuterProfile::AccelerationProfile(double r) const
{
    return -G * MassProfile(r) / pow(r, 2);
}

double OuterProfile::PotentialProfile(double r, double factor2) const
{
    double r_lower = R * factor * factor2;
    double lower = -(2.0 * M_PI * G * rho_mean * (pow(r_lower, 2) + 8.0 * pow(R, 1.5) * pow(factor, 1.5) * sqrt(r_lower))) / 3.0;
    double result = -(2.0 * M_PI * G * rho_mean * (pow(r, 2) + 8.0 * pow(R, 1.5) * pow(factor, 1.5) * sqrt(r))) / 3.0 - lower;
    return -result;
}

DoubleExponentialDisk::DoubleExponentialDisk(double Mass, double DiskScaleRadius, double DiskScaleHeight)
{
    Mass_disk = Mass;
    L = DiskScaleRadius;
    H = DiskScaleHeight;
}

double DoubleExponentialDisk::MassProfile(double r) const
{
    double r_transition = H;
    if(r <= r_transition)
    {
        return (Mass_disk / (2.0 * H * pow(L, 2))) * pow(r, 3);
    }
    else
    {
        double x = r / L;
        return Mass_disk * (1.0 - exp(-x) * (1.0 + x));
    }
}

double DoubleExponentialDisk::AccelerationProfile(double r) const
{
    return -G * MassProfile(r) / pow(r, 2);
}

double Gamma0(double x)
{
    // Gamma(0,x) = E1(x)
    // numerical approximation
    int N = 10000;
    double upper = 100.0;
    double dx = (upper - x) / N;
    double sum = 0.0;
    for(int i = 0; i < N; i++)
    {
        double x1 = x + i * dx;
        double x2 = x + (i + 1) * dx;
        double f1 = exp(-x1) / x1;
        double f2 = exp(-x2) / x2;
        sum += 0.5 * (f1 + f2) * dx;
    }
    return sum;
}

double GammaMinus1(double x)
{
 
    return (exp(-x) / x - Gamma0(x));
}

double DoubleExponentialDisk::PotentialProfile(double r) const
{
    double x = r / L;
    double gamma0 = Gamma0(x);
    double gamma_minus1 = GammaMinus1(x);
    return G * Mass_disk * ((gamma0 + gamma_minus1) * r - L) / (L * r);
}


double powerlaw_transition(double r, double r_min, double r_match, double value_min, double value_match)
{
    double alpha = std::log(value_match / value_min) / std::log(r_match / r_min);
    double A = value_min / std::pow(r_min, alpha);
    return A * std::pow(r, alpha);
}

void rho_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& rho, const std::vector<double>& Z, std::vector<std::vector<double>>& rho_grid, double R_cmax)
{
    int Nr = R.size();
    int NZ = Z.size();
    double r_match = 2.0 * R_cmax;
    int i_match = 0;
    int j_match = 0;
    double min_diff = 1e100;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::acos(Z[i] / r);
            double diff = std::abs(r - r_match);
            if(diff < min_diff)
            {
                min_diff = diff;
                i_match = i;
                j_match = j;
            }
        }
    }
    double rho_match_base = std::abs(rho[i_match][j_match]);
    double r_min = R_cmax * 1e-4;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::atan2(R[i], Z[j]);
            if(std::isnan(th))
            {
                rho_grid[i][j] = 0.0;
                continue;
            }
            double rho_r = std::abs(rho[i][j]);
            double angular = 11.0 / 4.0 * std::pow(std::sin(th), 2) - 35.0 / 24.0;
            double solution = rho_r * (1.0 + R_cmax * R_cmax / (r * r) * angular);
            double rho_match = rho_match_base * (1.0 + R_cmax * R_cmax / (r_match * r_match) * angular);
            double rho_min = rho_match / 100;
            if(r < r_match)
            {
                rho_grid[i][j] = powerlaw_transition(r, r_min, r_match, rho_min, rho_match);
            }
            else
            {
                rho_grid[i][j] = solution;
            }
        }
    }
}




void P_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& P, const std::vector<double>& Z, std::vector<std::vector<double>>& P_grid, double R_cmax)
{
    int Nr = R.size();
    int NZ = Z.size();
    double r_match = 2.0 * R_cmax;
    int i_match = 0;
    int j_match = 0;
    double min_diff = 1e100;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::acos(Z[i] / r);
            double diff = std::abs(r - r_match);
            if(diff < min_diff)
            {
                min_diff = diff;
                i_match = i;
                j_match = j;
            }
        }
    }
    double P_match_base = std::abs(P[i_match][j_match]);
    double r_min = R_cmax * 1e-4;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::atan2(R[i], Z[i]);
            if(std::isnan(th))
            {
                P_grid[i][j] = 0.0;
                continue;
            }
            double P_r = P[i][j];
            double angular = 4.0 / 3.0 * std::pow(std::sin(th), 2) - 5.0 / 8.0;
            double solution = P_r * (1.0 + R_cmax * R_cmax / (r * r) * angular);
            double P_match = P_match_base * (1.0 + R_cmax * R_cmax / (r_match * r_match) * angular);
            double P_min = P_match / 10000;
            if(r < r_match)
            {
                P_grid[i][j] = powerlaw_transition(r, r_min, r_match, P_min, P_match);
            }
            else
            {
                P_grid[i][j] = solution;
            }
        }
    }
}



void vr_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& velocity, const std::vector<double>& Z, std::vector<std::vector<double>>& vr_grid, double R_cmax)
{
    int Nr = R.size();
    int NZ = Z.size();
    double r_match = R_cmax;
    int i_match = 0;
    int j_match = 0;
    double min_diff = 1e100;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::acos(Z[i] / r);
            double diff = std::abs(r - r_match);
            if(diff < min_diff)
            {
                min_diff = diff;
                i_match = i;
                j_match = j;
            }
        }
    }
    double v_match_base = std::abs(velocity[i_match][j_match]);
    double r_min = R_cmax * 1e-4;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::atan2(R[i], Z[j]);
            if(std::isnan(th))
            {
                vr_grid[i][j] = 0.0;
                continue;
            }
            double v_r = -std::abs(velocity[i][j]);
            double angular = 23.0 / 12.0 * std::pow(std::sin(th), 2) - 65.0 / 72.0;
            double solution = v_r * (1.0 - R_cmax * R_cmax / (r * r) * angular);
            double v_match = v_match_base * (1.0 - R_cmax * R_cmax / (r_match * r_match) * angular);
            double v_min = v_match / 100.0;
            if(r < r_match)
            {
                vr_grid[i][j] = powerlaw_transition(r, r_min, r_match, v_min, v_match);
            }
            else
            {
                vr_grid[i][j] = solution;
            }
        }
    }
}




void vt_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& velocity, const std::vector<double>& Z, std::vector<std::vector<double>>& vt_grid, double R_cmax)
{
    int Nr = R.size();
    int NZ = Z.size();
    double r_match = R_cmax;
    int i_match = 0;
    int j_match = 0;
    double min_diff = 1e100;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::acos(Z[i] / r);
            double diff = std::abs(r - r_match);
            if(diff < min_diff)
            {
                min_diff = diff;
                i_match = i;
                j_match = j;
            }
        }
    }
    double v_match_base = std::abs(velocity[i_match][j_match]);
    double r_min = R_cmax * 1e-4;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::atan2(R[i], Z[j]);
            if(std::isnan(th))
            {
                vt_grid[i][j] = 0.0;
                continue;
            }
            double v_r = -std::abs(velocity[i][j]);
            double solution = v_r * 5.0 / 18.0 * R_cmax * R_cmax / (r * r) * std::sin(2.0 * th);
            double v_match = v_match_base * 5.0 / 18.0 * std::sin(2.0 * th);
            double v_min = v_match / 100.0;
            if(r < r_match)
            {
                vt_grid[i][j] = powerlaw_transition(r, r_min, r_match, v_min, v_match);
            }
            else
            {
                vt_grid[i][j] = solution;
            }
        }
    }
}



void vp_2D(const std::vector<double>& R, const std::vector<std::vector<double>>& M, const std::vector<double>& Z, std::vector<std::vector<double>>& vp_grid, double R_cmax)
{
    int Nr = R.size();
    int NZ = Z.size();
    double r_match = R_cmax;
    int i_match = 0;
    int j_match = 0;
    double min_diff = 1e100;
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::acos(Z[i] / r);
            double diff = std::abs(r - r_match);
            if(diff < min_diff)
            {
                min_diff = diff;
                i_match = i;
                j_match = j;
            }
        }
    }
    double M_match_base = std::abs(M[i_match][j_match]);
    double r_min = R_cmax * 1e-4;
    double M_match = M[i_match][j_match];
    double vc_match = v_c(M_match, r_match);
    double omega_match = vc_match * R_cmax / (r_match * r_match);
    for(int i = 0; i < Nr; i++)
    {
        for(int j = 0; j < NZ; j++)
        {
            double r = std::sqrt(R[i] * R[i] + Z[j] * Z[j]);
            double th = std::atan2(R[i], Z[j]);
            if(std::isnan(th))
            {
                vp_grid[i][j] = 0.0;
                continue;
            }
            double vc = v_c(M[i][j], r);
            double omega = vc * R_cmax / (r * r);
            double solution = omega * r * std::sin(th);
            double vp_match = omega_match * r_match * std::sin(th);
            double vp_min = vp_match / 100.0;
            if(r < r_match)
            {
                vp_grid[i][j] = powerlaw_transition(r, r_min, r_match, vp_min, vp_match);
            }
            else
            {
                vp_grid[i][j] = solution;
            }
        }
    }
}

double ellipsoid_smoothing_factor(double x, double y, double z, double R_scale, double H_scale, double m_cutoff)
{
    double dx = x;
    double dy = y;
    double dz = z;
    double R = std::sqrt(dx * dx + dy * dy);
    double m = std::sqrt((R / R_scale) * (R / R_scale) + (dz / H_scale) * (dz / H_scale));
    const double m_min = 1.0;
    if(m >= m_cutoff) return 1.0;
    if(m <= m_min) return 0.0;
    const double log_factor_min = -10.0;
    const double exponent = log_factor_min * (1.0 - std::log10(m / m_min) / std::log10(m_cutoff / m_min));
    return std::pow(10.0, exponent);
}

double estimate_max_turbulence_percentage(size_t N_cells, double max_negative_cells)
{
    if(N_cells == 0) return 0.0;
    auto expected_negative = [N_cells](double p)
    {
        const double sigma_ratio = std::sqrt(p);
        const double probability = 0.5 * std::erfc(1.0 / (std::sqrt(2.0) * sigma_ratio));
        return static_cast<double>(N_cells) * probability;
    };
    double p_low = 0.0;
    double p_high = 1.0;
    if(expected_negative(p_high) < max_negative_cells) return p_high;
    for(int iter = 0; iter < 100; iter++)
    {
        const double p_mid = 0.5 * (p_low + p_high);
        if(expected_negative(p_mid) < max_negative_cells)
        {
            p_low = p_mid;
        }
        else
        {
            p_high = p_mid;
        }
    }
    return p_low;
}

// double spherical_theta(
//     double r,
//     double z
// )
// {
//     if(r <= 0.0)
//     {
//         return 0.0;
//     }


//     double cos_theta = z / r;


//     // avoid numerical errors slightly outside [-1,1]
//     cos_theta = std::max(
//         -1.0,
//         std::min(1.0, cos_theta)
//     );


//     return std::acos(cos_theta);
// }