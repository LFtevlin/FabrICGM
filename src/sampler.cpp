#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include "sampler.hpp"
#include "simulation.hpp"

std::vector<double> linspace(double start, double end, int N)
{
    std::vector<double> result(N);
    double step = (end - start) / (N - 1);
    for (int i = 0; i < N; i++)
    {
        result[i] = start + i * step;
    }
    return result;
}

double interpolate_density(
    double Z,
    double R,
    const std::vector<double>& Zcenter,
    const std::vector<double>& Rcenter,
    const std::vector<std::vector<double>>& density
)
{
    int Nz = Zcenter.size();
    int Nr = Rcenter.size();
    auto itZ = std::lower_bound(Zcenter.begin(), Zcenter.end(), Z);
    int i = std::distance(Zcenter.begin(), itZ) - 1;
    i = std::max(0, std::min(i, Nz - 2));
    auto itR = std::lower_bound(Rcenter.begin(), Rcenter.end(), R);
    int j = std::distance(Rcenter.begin(), itR) - 1;
    j = std::max(0, std::min(j, Nr - 2));
    double z1 = Zcenter[i];
    double z2 = Zcenter[i+1];
    double r1 = Rcenter[j];
    double r2 = Rcenter[j+1];
    double fz = (Z-z1)/(z2-z1);
    double fr = (R-r1)/(r2-r1);
    double f11 = density[i][j];
    double f12 = density[i][j+1];
    double f21 = density[i+1][j];
    double f22 = density[i+1][j+1];
    return f11*(1-fz)*(1-fr) + f21*fz*(1-fr) + f12*(1-fz)*fr + f22*fz*fr;
}

double random_uniform(
    double a,
    double b,
    std::mt19937 &rng
)
{
    std::uniform_real_distribution<double> dist(a,b);
    return dist(rng);
}

void sample_one_point(
    double r_min, double r_max, double theta_min, double theta_max, double rho_max, std::mt19937 &rng,
    const std::vector<double>& Zcenter, const std::vector<double>& Rcenter, const std::vector<std::vector<double>>& density,
    std::vector<double>& r_sample, std::vector<double>& theta_sample, std::vector<double>& phi_sample
)
{
    if (rho_max <= 0)
        return;

    bool check = true;
    int attempts = 0;

    while (check)
    {
        attempts++;

        if (attempts > 100000)
        {
            throw std::runtime_error("Too many rejection sampling attempts");
        }

        double mu_rnd = random_uniform(std::cos(theta_min), std::cos(theta_max), rng);
        double theta_rnd = std::acos(mu_rnd);
        double r_rnd = std::cbrt(random_uniform(r_min*r_min*r_min, r_max*r_max*r_max, rng));
        double Z = std::abs(r_rnd * std::cos(theta_rnd));
        double R = r_rnd * std::sin(theta_rnd);
        double rho_val = interpolate_density(Z, R, Zcenter, Rcenter, density);
        double a = random_uniform(0.0, 1.0, rng);

        if (a <= rho_val/rho_max)
        {
            check = false;
            r_sample.push_back(r_rnd);
            theta_sample.push_back(theta_rnd);
            double phi = random_uniform(0.0, 2.0*M_PI, rng);
            phi_sample.push_back(phi);
        }
    }
}



SampledPositions sample_equal_mass(const std::vector<std::vector<double>>& gas_density, const std::vector<double>& Zcenter, const std::vector<double>& Rcenter, double target_mass, int NDisk)
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    if (gas_density.empty() || gas_density[0].empty())
    {
        throw std::runtime_error("gas_density must be a 2D non-empty grid");
    }
    if (Zcenter.empty() || Rcenter.empty())
    {
        throw std::runtime_error("Zcenter and Rcenter must be non-empty");
    }
    std::random_device rd;
    std::mt19937 rng(rd());
    std::vector<double> rbounds(Rcenter.size());
    if (Rcenter.size() != Zcenter.size())
        throw std::runtime_error("Rcenter and Zcenter must have same size");
    for (size_t i = 0; i < Rcenter.size(); i++)
    {
        rbounds[i] = std::sqrt(Rcenter[i]*Rcenter[i] + Zcenter[i]*Zcenter[i]);
    }
    int N_data = 1024;
    double r_min = *std::min_element(rbounds.begin(), rbounds.end());
    double r_max = *std::max_element(rbounds.begin(), rbounds.end());
    std::vector<double> r_data = linspace(r_min, r_max, N_data);
    std::vector<double> dr(N_data-1);
    for (int i = 0; i < N_data-1; i++)
    {
        dr[i] = r_data[i+1] - r_data[i];
    }
    std::vector<double> r_center_data(N_data-1);
    for (int i = 0; i < N_data-1; i++)
    {
        r_center_data[i] = r_data[i] + 0.5 * dr[i];
    }
    std::vector<double> theta_data = linspace(0.0, M_PI, N_data);
    std::vector<double> dtheta(N_data-1);
    for (int i = 0; i < N_data-1; i++)
    {
        dtheta[i] = theta_data[i+1] - theta_data[i];
    }
    std::vector<double> theta_center_data(N_data-1);
    for (int i = 0; i < N_data-1; i++)
    {
        theta_center_data[i] = theta_data[i] + 0.5 * dtheta[i];
    }
    std::vector<double> r_sample;
    std::vector<double> phi_sample;
    std::vector<double> theta_sample;
    r_sample.reserve(NDisk);
    theta_sample.reserve(NDisk);
    phi_sample.reserve(NDisk);
    std::vector<double> M_sample;
    for (int k = 0; k < N_data - 1; k++)
    {
        for (int i = 0; i < N_data - 1; i++)
        {
            double r_min, r_max;
            r_min = r_center_data[k] - dr[k]/2;
            r_max = r_center_data[k] + dr[k]/2;
            double theta_min, theta_max;
            theta_min = theta_center_data[i] - dtheta[i]/2;
            theta_max = theta_center_data[i] + dtheta[i]/2;
            const int N_sub = 5;
            double dr_local = (r_max - r_min) / (N_sub - 1);
            double dtheta_local = (theta_max - theta_min) / (N_sub - 1);
            double M_element_int = 0.0;
            double rho_max = 0.0;
            double R = 0.0;
            double Z = 0.0;
            for (int ii = 0; ii < N_sub - 1; ii++)
            {
                double r_sub_center = r_min + (ii + 0.5) * dr_local;
                for (int jj = 0; jj < N_sub - 1; jj++)
                {
                    double theta_sub_center = theta_min + (jj + 0.5) * dtheta_local;
                    Z = std::abs(r_sub_center * std::cos(theta_sub_center));
                    R = r_sub_center * std::sin(theta_sub_center);
                    double rho0 = interpolate_density(Z, R, Zcenter, Rcenter, gas_density);
                    M_element_int += 2.0 * M_PI * r_sub_center * r_sub_center * dr_local * std::sin(theta_sub_center) * rho0 * dtheta_local;
                    rho_max = std::max(rho_max, rho0);
                }
            }
            double N_element = M_element_int / target_mass;
            int N_sample = (int)(N_element);
            double remain = N_element - N_sample;
            double random = random_uniform(0.0, 1.0, rng);
            if (remain >= random)
            {
                sample_one_point(r_min, r_max, theta_min, theta_max, rho_max, rng, Zcenter, Rcenter, gas_density, r_sample, theta_sample, phi_sample);
            }
            if (N_element >= 1)
            {
                for (int m = 0; m < N_sample; m++)
                {
                    sample_one_point(r_min, r_max, theta_min, theta_max, rho_max, rng, Zcenter, Rcenter, gas_density, r_sample, theta_sample, phi_sample);
                }
            }
        }
    }
    std::vector<double> x_sample(r_sample.size());
    std::vector<double> y_sample(r_sample.size());
    std::vector<double> z_sample(r_sample.size());
    for (int i = 0; i < r_sample.size(); i++)
    {
        double r = r_sample[i];
        double theta = theta_sample[i];
        double phi = phi_sample[i];
        x_sample[i] = r * std::sin(theta) * std::cos(phi);
        y_sample[i] = r * std::sin(theta) * std::sin(phi);
        z_sample[i] = r * std::cos(theta);
    }
    SampledPositions result;
    result.x = x_sample;
    result.y = y_sample;
    result.z = z_sample;
    return result;
}

double jitter(double center, double dx, std::mt19937& rng, double strength=0.5
)
{
    std::uniform_real_distribution<double> dist(-0.5,0.5);

    return center + dist(rng)*dx*strength;
}

SampledPositions sample_cartesian(const std::vector<double>& LBox, const std::vector<double>& dx, const SampledPositions* initial)
{
    SampledPositions result;
    if(initial != nullptr)
    {
        result.x = initial->x;
        result.y = initial->y;
        result.z = initial->z;
    }
    std::mt19937 rng(std::random_device{}());
    for(size_t b = 0; b < LBox.size(); b++)
    {
        double L = LBox[b];
        double cell = dx[b];
        int N = static_cast<int>(L/cell);
        int Ncell = N;
        std::vector<int> counts(Ncell*Ncell*Ncell, 0);
        for(size_t p = 0; p < result.x.size(); p++)
        {
            int ix = std::floor((result.x[p]+L/2.0)/cell);
            int iy = std::floor((result.y[p]+L/2.0)/cell);
            int iz = std::floor((result.z[p]+L/2.0)/cell);
            if(ix >= 0 && ix < Ncell && iy >= 0 && iy < Ncell && iz >= 0 && iz < Ncell)
            {
                int id = ix*Ncell*Ncell + iy*Ncell + iz;
                counts[id]++;
            }
        }
        for(int ix = 0; ix < Ncell; ix++)
        {
            for(int iy = 0; iy < Ncell; iy++)
            {
                for(int iz = 0; iz < Ncell; iz++)
                {
                    int id = ix*Ncell*Ncell + iy*Ncell + iz;
                    if(counts[id] == 0)
                    {
                        double xc = -L/2.0 + (ix+0.5)*cell;
                        double yc = -L/2.0 + (iy+0.5)*cell;
                        double zc = -L/2.0 + (iz+0.5)*cell;
                        result.x.push_back(jitter(xc, cell, rng));
                        result.y.push_back(jitter(yc, cell, rng));
                        result.z.push_back(jitter(zc, cell, rng));
                    }
                }
            }
        }
    }
    return result;
}
