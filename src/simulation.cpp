#include "simulation.hpp"
#include <fftw3.h>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <H5Cpp.h>
#include <iostream>
#include "constants.hpp"
#include "sampler.hpp"
#include "functions.hpp"

std::ofstream logfile;

std::vector<double> Simulation::create_radius_grid(double rmin, double rmax, int N)
{
    std::vector<double> r(N);
    double log_min = std::log10(rmin);
    double log_max = std::log10(rmax);
    for(int i = 0; i < N; i++)
    {
        double x = static_cast<double>(i)/(N-1);
        r[i] = std::pow(10.0, log_min + x*(log_max-log_min));
    }
    return r;
}

void Simulation::update_initial_state()
{
    v_init = -MACH_INIT*c_s(T_init);
    rho_init = params.galaxy.Mdot/(4.0*M_PI*params.galaxy.Rsonic*params.galaxy.Rsonic*(-v_init));
}

void Simulation::print_initial_conditions() const
{
    logfile << "\n===== Initial Conditions =====\n";
    logfile << "R200        = " << halo->R200/kpc_to_cm << " kpc\n";
    logfile << "boxsize     = " << params.boxsize/kpc_to_cm << " kpc\n";
    logfile << "Rmax        = " << Rmax/kpc_to_cm << " kpc\n";
    logfile << "sonic_index = " << sonic_index << "\n";
    logfile << "Rsonic      = " << radius[sonic_index]/kpc_to_cm << " kpc\n";
    logfile << "T_init      = " << T_init << " K\n";
    logfile << "v_init      = " << v_init/1e5 << " km/s\n";
    logfile << "n_init    = " << rho_init/(m_p*0.597) << " 1/cm^3\n";
    logfile << "==============================\n";
}

void Simulation::create_profiles()
{
    halo = std::make_unique<HernquistHalo>(params.M200, params.galaxy.concentration, params.galaxy.fgas, params.z);

    if(params.boxsize < 0)
    {
        params.boxsize = 3.0*halo->R200;
        logfile << "==============================\n";
        logfile << "Warning: boxsize not specified.\n ";
        logfile << "Using boxsize = 3 R200.\n";
        logfile << "==============================\n";
    }

    Rmax = std::sqrt(3.0)*params.boxsize;

    radius = create_radius_grid(0.1*kpc_to_cm, Rmax + 100*kpc_to_cm, 256);

    gas_disk = std::make_unique<DoubleExponentialDisk>(params.galaxy.Mgas, params.galaxy.Rgas, params.galaxy.Hgas);

    stellar_disk = std::make_unique<DoubleExponentialDisk>(params.galaxy.Mstellar, params.galaxy.Rstar, params.galaxy.Hstar);

    bulge = std::make_unique<HernquistBulge>(params.galaxy.Mstellar*params.BulgeFraction, params.BulgeScaleLength*params.galaxy.Rstar);

    outer = std::make_unique<OuterProfile>(halo->R200, params.z);

    size_t N = radius.size();

    Z.resize(N);
    M_total.resize(N);
    potential.resize(N);

    double r0 = 2.0*params.galaxy.Rgas;

    for(size_t i = 0; i < N; i++)
    {
        if(params.Z_profile == "constant")
        {
            Z[i] = params.galaxy.Z0;
        }
        else if(params.Z_profile == "powerlaw")
        {
            Z[i] = params.galaxy.Z0*std::pow(radius[i]/r0, params.Z_slope);
        }
        else
        {
            throw std::runtime_error("Unknown metallicity profile");
        }
    }

    for(size_t i = 0; i < N; i++)
    {
        double r = radius[i];
        M_total[i] = halo->MassProfile(r) + outer->MassProfile(r) + gas_disk->MassProfile(r) + stellar_disk->MassProfile(r) + bulge->MassProfile(r);
        potential[i] = halo->PotentialProfile(r) + outer->PotentialProfile(r) + gas_disk->PotentialProfile(r) + stellar_disk->PotentialProfile(r) + bulge->PotentialProfile(r);
    }

    logfile << "\n===== Galaxy profile parameters =====\n";
    logfile << "M200          = " << params.M200/Msun_to_g << " Msun\n";
    logfile << "concentration = " << params.galaxy.concentration << "\n";
    logfile << "Mgas          = " << params.galaxy.Mgas/Msun_to_g << " Msun\n";
    logfile << "Mstellar      = " << params.galaxy.Mstellar/Msun_to_g << " Msun\n";
    logfile << "fgas          = " << params.galaxy.fgas << "\n";
    logfile << "Rgas          = " << params.galaxy.Rgas/kpc_to_cm << " kpc\n";
    logfile << "Rstar         = " << params.galaxy.Rstar/kpc_to_cm << " kpc\n";
    logfile << "Hgas          = " << params.galaxy.Hgas/kpc_to_cm << " kpc\n";
    logfile << "Hstar         = " << params.galaxy.Hstar/kpc_to_cm << " kpc\n";
    logfile << "Mdot          = " << params.galaxy.Mdot/Msun_per_year_to_g_per_s << " Msun/yr\n";
    logfile << "Z0            = " << params.galaxy.Z0 << " Zsun\n";
    logfile << "R200          = " << halo->R200/kpc_to_cm << " kpc\n";
    logfile << "Rmax          = " << Rmax/kpc_to_cm << " kpc\n";
    logfile << "Rsonic        = " << params.galaxy.Rsonic/kpc_to_cm << " kpc\n";
    logfile << "====================================\n\n";
}

void Simulation::initialize_conditions()
{
    sonic_index = 0;
    double min_dist = 1e30;
    for(size_t i = 0; i < radius.size(); i++)
    {
        double d = std::abs(radius[i] - params.galaxy.Rsonic);
        if(d < min_dist)
        {
            min_dist = d;
            sonic_index = i;
        }
    }
    logfile << "Rsonic found = " << radius[sonic_index]/kpc_to_cm << " kpc\n" << "Rsonic from params = " << params.galaxy.Rsonic/kpc_to_cm << " kpc\n";
    T_init = T_c(M_total[sonic_index], radius[sonic_index]);
    logfile << "Initial temperature guess based on circular temperature T_init = " << T_init << " K\n";
    update_initial_state();
    auto result = find_Tinit(radius, M_total, T_init, rho_init, v_init, params.galaxy.Mdot, MACH_INIT, params.galaxy.Rsonic, sonic_index, epsilon, Z, params.z, cooling, 400);
    T_init = std::get<0>(result);
    v_init = std::get<1>(result);
    rho_init = std::get<2>(result);
}

#include <chrono>
#include <cmath>

void Simulation::write_catalogue_result(bool reached_Rmax)
{
    std::cout << "CATALOGUE_RESULT " << params.M200/Msun_to_g << " " << params.z << " " << reached_Rmax << " " << radius_max/kpc_to_cm << " " << params.galaxy.concentration << " " << params.galaxy.Z0 << " " << params.galaxy.Mgas/Msun_to_g << " " << params.galaxy.Mstellar/Msun_to_g << " " << params.galaxy.Rgas/kpc_to_cm << " " << params.galaxy.Rstar/kpc_to_cm << " " << params.galaxy.Hgas/kpc_to_cm << " " << params.galaxy.Hstar/kpc_to_cm << " " << params.galaxy.Mdot/Msun_per_year_to_g_per_s << " " << params.galaxy.Rsonic/kpc_to_cm << " ";
    if(reached_Rmax)
    {
        std::cout << Mcgm_R200/Msun_to_g;
    }
    else
    {
        std::cout << -1;
    }
    std::cout << std::endl;
}

void Simulation::compute_CGM_mass()
{
    double M_R200 = 0.0;
    double M_Rmax = 0.0;
    int N = radius.size();
    M_CGM.resize(N);
    for(size_t i = sonic_index; i < radius.size() - 1; i++)
    {
        double r1 = radius[i];
        double r2 = radius[i + 1];
        double rho1 = density[i];
        double rho2 = density[i + 1];
        double integrand1 = 4.0*M_PI*r1*r1*rho1;
        double integrand2 = 4.0*M_PI*r2*r2*rho2;
        if(r1 < halo->R200)
        {
            double r2_use = std::min(r2, halo->R200);
            double f = (r2_use - r1)/(r2 - r1);
            double integrand2_use = integrand1 + f*(integrand2 - integrand1);
            M_R200 += 0.5*(integrand1 + integrand2_use)*(r2_use - r1);
        }
        if(r1 >= halo->R200)
        {
            Mcgm_R200 = M_R200;
        }
        if(r1 < Rmax)
        {
            double r2_use = std::min(r2, Rmax);
            double f = (r2_use - r1)/(r2 - r1);
            double integrand2_use = integrand1 + f*(integrand2 - integrand1);
            M_Rmax += 0.5*(integrand1 + integrand2_use)*(r2_use - r1);
            M_CGM[i] = M_Rmax;
        }
        if(r1 >= Rmax)
        {
            Mcgm_Rmax = M_Rmax;
            break;
        }
    }
    if(!std::isfinite(M_R200))
    {
        logfile << "\n===== CGM Mass =====\n" << "Computed mass is NaN or Inf.\n" << "====================\n";
        throw std::runtime_error("Computed CGM mass is not finite.");
    }
    logfile << "\n===== CGM Mass =====\n" << "MCGM(<R200) = " << M_R200/Msun_to_g << " Msun\n" << "====================\n";
}


void Simulation::integrate()
{
    int N = radius.size();
    radius_max = params.galaxy.Rsonic;
    velocity.resize(N);
    temperature.resize(N);
    density.resize(N);
    Lambda_fctn.resize(N);
    auto reset_solution = [&]()
    {
        std::fill(Lambda_fctn.begin(), Lambda_fctn.end(), 0.0);
        std::fill(velocity.begin(), velocity.end(), 0.0);
        std::fill(temperature.begin(), temperature.end(), 0.0);
        std::fill(density.begin(), density.end(), 0.0);
        velocity[sonic_index] = v_init;
        temperature[sonic_index] = T_init;
        density[sonic_index] = rho_init;
    };
    reset_solution();
    std::vector<double> dr_list(N-1);
    for(int i = 0; i < N-1; i++)
        dr_list[i] = radius[i+1]-radius[i];
    int i = 0;
    auto start = std::chrono::steady_clock::now();
    const double time_limit = params.runtime;
    while(true)
    {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now-start).count();
        if(elapsed > time_limit)
        {
            if(params.validate_galaxies)
            {
                write_catalogue_result(false);
                throw std::runtime_error("Catalogue entry written");
            }
            throw std::runtime_error("No solution found");
        }
        int idx = sonic_index+i;
        if(idx >= N)
        {
            throw std::runtime_error("Integration exceeded radius grid");
        }
        double T = temperature[idx];
        double rho = density[idx];
        double v = velocity[idx];
        double r = radius[idx];
        Lambda_fctn[idx] = cooling.Lambda(rho, T, Z[idx], params.z);
        bool mach = check_Mach(v, T);
        bool bernoulli = check_Bernoulli(v, T, potential[idx]);
        if(!std::isfinite(T) || !std::isfinite(v) || !std::isfinite(rho))
        {
            if(params.validate_galaxies)
            {
                write_catalogue_result(false);
                throw std::runtime_error("Catalogue entry written");
            }
            logfile << "Non-finite state reached\n";
            logfile << "radius = " << radius[idx]/kpc_to_cm << " kpc\n";
            throw std::runtime_error("Integration produced NaN");
        }
        if(elapsed > time_limit/2.0 && r > 3*halo->R200)
        {
            bernoulli = true;
        }
        if(!mach)
        {
            T_init *= 1.05;
            update_initial_state();
            reset_solution();
            i = 0;
            continue;
        }
        if(!bernoulli)
        {
            T_init *= 0.95;
            update_initial_state();
            reset_solution();
            double B = Bernoulli(v, T, potential[idx]);
            i = 0;
            continue;
        }
        try
        {
            auto result = RK4(r, T, rho, v, M_total, radius, dr_list[idx], params.galaxy.Mdot, Z, params.z, cooling, params.tolerance, max_iter);
            velocity[idx+1] = result.first;
            temperature[idx+1] = result.second;
            density[idx+1] = params.galaxy.Mdot/(4.0*M_PI*radius[idx+1]*radius[idx+1]*(-velocity[idx+1]));
        }
        catch(const std::exception& e)
        {
            logfile << "\n===== RK4 error =====\n" << "radius = " << radius[idx]/kpc_to_cm << " kpc\n" << "Error: " << e.what() << "\n" << "=====================\n";
            T_init *= 1.05;
            update_initial_state();
            reset_solution();
            i = 0;
            continue;
        }
        radius_max = std::max(radius_max, radius[sonic_index+i+1]);
        i++;
        if(radius[sonic_index+i] >= Rmax)
        {
            compute_CGM_mass();
            logfile << "Reached fbaryon total = " << (params.galaxy.Mstellar + params.galaxy.Mgas + Mcgm_Rmax)/params.M200 << '\n';
            logfile << "Reached fbaryon <R200 = " << (params.galaxy.Mstellar + params.galaxy.Mgas + Mcgm_R200)/params.M200 << '\n';
            if(params.validate_galaxies)
            {
                write_catalogue_result(true);
                throw std::runtime_error("Catalogue entry written");
            }
            logfile << "\n===== Integration successful =====\n" << "Reached Rmax\n" << "=================================\n";
            break;
        }
    }
}

void Simulation::convert_1D_to_3D()
{
    double R_cmax = params.galaxy.Rsonic;
    int N_R = 256;
    double R_min = params.galaxy.Rgas/1e5;
    double R_max = radius.back();
    std::vector<double> R_new = create_radius_grid(R_min, R_max, N_R);
    int N_Z = 256;
    double Z_min = params.galaxy.Hgas/1e5;
    double Z_max = radius.back();
    std::vector<double> Z_new = create_radius_grid(Z_min, Z_max, N_Z);
    std::vector<std::vector<double>> rho_new(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<std::vector<double>> P_new(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<std::vector<double>> v_new(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<std::vector<double>> M_new(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<double> pressure(temperature.size());
    for(size_t i = 0; i < temperature.size(); i++)
    {
        pressure[i] = density[i]*k_B*temperature[i]/(mu*m_p);
    }
    for(int i = 0; i < N_R; i++)
    {
        double R = R_new[i];
        for(int j = 0; j < N_Z; j++)
        {
            double Z = Z_new[j];
            double r = std::sqrt(R*R + Z*Z);
            rho_new[i][j] = interpolate(r, radius, density);
            P_new[i][j] = interpolate(r, radius, pressure);
            v_new[i][j] = interpolate(r, radius, velocity);
            M_new[i][j] = interpolate(r, radius, M_total);
        }
    }
    std::vector<std::vector<double>> rho_grid(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<std::vector<double>> P_grid(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<std::vector<double>> vr_grid(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<std::vector<double>> vt_grid(N_R, std::vector<double>(N_Z, 0.0));
    std::vector<std::vector<double>> vp_grid(N_R, std::vector<double>(N_Z, 0.0));
    std::cerr << "Expected Cartesian cells = " << N_R*N_Z << "\n";
    rho_2D(R_new, rho_new, Z_new, rho_grid, R_cmax);
    P_2D(R_new, P_new, Z_new, P_grid, R_cmax);
    vr_2D(R_new, v_new, Z_new, vr_grid, R_cmax);
    vt_2D(R_new, v_new, Z_new, vt_grid, R_cmax);
    vp_2D(R_new, M_new, Z_new, vp_grid, R_cmax);
    std::cerr << "Created 2D profiles \n";
    logfile << params.sampling;
    for(int i = 0; i < params.NBox; i++)
    {
        logfile << "Box size: ";
        logfile << params.LBox[i];
        logfile << ", dx =  ";
        logfile << params.dx[i] << "\n";
    }
    if(params.sampling == "equalmass" || params.sampling == "both")
    {
        logfile << "Using equal mass sampling\n";
        int NDisk = Mcgm_Rmax/params.mtarget;
        gas = sample_equal_mass(rho_grid, Z_new, R_new, params.mtarget, NDisk);
        std::cerr << "sampled equal mass \n";
    }
    if(params.sampling == "cartesian" || params.sampling == "both")
    {
        logfile << "Using cartesian sampling\n";
        SampledPositions cart = sample_cartesian(params.LBox, params.dx, &gas);
        gas.x.insert(gas.x.end(), cart.x.begin(), cart.x.end());
        gas.y.insert(gas.y.end(), cart.y.begin(), cart.y.end());
        gas.z.insert(gas.z.end(), cart.z.begin(), cart.z.end());
        std::cerr << "sampled cartesian boxes \n";
    }
    double xmin = -params.boxsize;
    double xmax = params.boxsize;
    std::vector<double> x_new;
    std::vector<double> y_new;
    std::vector<double> z_new;
    x_new.reserve(gas.x.size());
    y_new.reserve(gas.y.size());
    z_new.reserve(gas.z.size());
    for(size_t i = 0; i < gas.x.size(); i++)
    {
        if(gas.x[i] >= xmin && gas.x[i] <= xmax && gas.y[i] >= xmin && gas.y[i] <= xmax && gas.z[i] >= xmin && gas.z[i] <= xmax)
        {
            x_new.push_back(gas.x[i]);
            y_new.push_back(gas.y[i]);
            z_new.push_back(gas.z[i]);
        }
    }
    gas.x.swap(x_new);
    gas.y.swap(y_new);
    gas.z.swap(z_new);
    gas.rho.resize(gas.x.size());
    gas.P.resize(gas.x.size());
    gas.v.resize(gas.x.size(), {0.0, 0.0, 0.0});
    gas.B.resize(gas.x.size(), {0.0, 0.0, 0.0});
    double B0_ISM = 0.0;
    double B0_CGM = 0.0;
    double rho_max = 0.0;
    if(params.magnetic)
    {
        B0_ISM = 1e-1*muG*std::pow(params.M200/Msun_to_g/1e12, 1.0/3.0);
        B0_CGM = 1e-3*muG*std::pow(params.M200/Msun_to_g/1e12, 1.0/3.0);
        for(const auto& row : rho_grid)
        {
            if(!row.empty())
            {
                rho_max = std::max(rho_max, *std::max_element(row.begin(), row.end()));
            }
        }
    }
    double vp, vt, vr;
    gas.T.resize(gas.x.size());
    for(size_t i = 0; i < gas.x.size(); i++)
    {
        double R = std::sqrt(gas.x[i]*gas.x[i] + gas.y[i]*gas.y[i]);
        double Z = std::abs(gas.z[i]);
        gas.rho[i] = interpolate_cylindrical_grid(R, Z, R_new, Z_new, rho_grid);
        gas.P[i] = interpolate_cylindrical_grid(R, Z, R_new, Z_new, P_grid);
        gas.T[i] = gas.P[i]*mu*m_p/(gas.rho[i]*k_B);
        vp = interpolate_cylindrical_grid(R, Z, R_new, Z_new, vp_grid);
        vr = interpolate_cylindrical_grid(R, Z, R_new, Z_new, vr_grid);
        vt = interpolate_cylindrical_grid(R, Z, R_new, Z_new, vt_grid);
        double r = std::sqrt(gas.x[i]*gas.x[i] + gas.y[i]*gas.y[i] + gas.z[i]*gas.z[i]);
        double theta = std::acos(gas.z[i]/r);
        double phi = std::atan2(gas.y[i], gas.x[i]);
        gas.v[i][0] = std::sin(theta)*std::cos(phi)*vr + std::cos(theta)*std::cos(phi)*vt - std::sin(phi)*vp;
        gas.v[i][1] = std::sin(theta)*std::sin(phi)*vr + std::cos(theta)*std::sin(phi)*vt + std::cos(phi)*vp;
        gas.v[i][2] = std::cos(theta)*vr - std::sin(theta)*vt;
        if(gas.z[i] > 0.0)
        {
            gas.v[i][2] = -std::abs(gas.v[i][2]);
        }
        else if(gas.z[i] < 0.0)
        {
            gas.v[i][2] = std::abs(gas.v[i][2]);
        }
        else
        {
            gas.v[i][2] = 0.0;
        }
        if(params.magnetic)
        {
            double factor = std::min(1.0, std::sqrt(gas.rho[i]/rho_max));
            double Bmag = B0_ISM/std::sqrt(3.0)*factor;
            if(R > 0.0)
            {
                gas.B[i][0] = Bmag*gas.y[i]/R;
                gas.B[i][1] = -Bmag*gas.x[i]/R;
            }
            else
            {
                gas.B[i][0] = 0.0;
                gas.B[i][1] = 0.0;
            }
            gas.B[i][2] = B0_CGM*factor;
        }
    }
    std::cerr << "Interpolated on cartesian grid \n";
    logfile << "1D to 3D conversion done.\n";
}

SampledPositions Simulation::create_turbulent_box(double k_inj, double LBox, double dxBox)
{
    logfile << ": LBox = " << LBox/kpc_to_cm << " kpc, dx = " << dxBox/kpc_to_cm << " kpc\n";
    const int N = static_cast<int>(LBox/dxBox);
    logfile << "N = " << N << "\n";
    const double L_turb = LBox;
    SampledPositions result_Box;
    const size_t total = static_cast<size_t>(N)*static_cast<size_t>(N)*static_cast<size_t>(N);
    result_Box.x.resize(total);
    result_Box.y.resize(total);
    result_Box.z.resize(total);
    result_Box.v.resize(total, {0.0, 0.0, 0.0});
    result_Box.B.resize(total, {0.0, 0.0, 0.0});
    result_Box.rho.resize(total, 0.0);
    result_Box.T.resize(total, 0.0);
    size_t index = 0;
    for(int ix = 0; ix < N; ix++)
    {
        for(int iy = 0; iy < N; iy++)
        {
            for(int iz = 0; iz < N; iz++)
            {
                result_Box.x[index] = -LBox/2.0 + (ix + 0.5)*dxBox;
                result_Box.y[index] = -LBox/2.0 + (iy + 0.5)*dxBox;
                result_Box.z[index] = -LBox/2.0 + (iz + 0.5)*dxBox;
                index++;
            }
        }
    }
    const int N_expected = static_cast<int>(params.LBox.back()/dxBox);
    const double norm = static_cast<double>(N)*N*N/(static_cast<double>(N_expected)*N_expected*N_expected);
    logfile << "normalization = " << norm << "\n";
    std::vector<std::array<double, 3>> e_rho = ::create_turbulence(L_turb, N, k_inj, -5.0/3.0, false, false);
    std::vector<std::array<double, 3>> e_v = ::create_turbulence(L_turb, N, k_inj, -5.0/3.0, true, false);
    const size_t Nf = e_v.size();
    if(Nf != total)
    {
        throw std::runtime_error("Turbulence field size does not match N^3");
    }
    std::vector<std::array<double, 3>> e_B;
    if(params.magnetic)
    {
        e_B = ::create_turbulence(L_turb, N, k_inj, -5.0/3.0, true, true);
    }
    for(auto& value : e_rho)
    {
        value[0] *= norm;
    }
    for(auto& value : e_v)
    {
        value[0] *= norm;
        value[1] *= norm;
        value[2] *= norm;
    }
    if(params.magnetic)
    {
        for(auto& value : e_B)
        {
            value[0] *= norm;
            value[1] *= norm;
            value[2] *= norm;
        }
    }
    for(size_t i = 0; i < total; i++)
    {
        const double E_mag_v = std::sqrt(e_v[i][0]*e_v[i][0] + e_v[i][1]*e_v[i][1] + e_v[i][2]*e_v[i][2]);
        if(E_mag_v > 0.0)
        {
            const double v_mag = std::sqrt(E_mag_v);
            result_Box.v[i][0] = v_mag*e_v[i][0]/E_mag_v;
            result_Box.v[i][1] = v_mag*e_v[i][1]/E_mag_v;
            result_Box.v[i][2] = v_mag*e_v[i][2]/E_mag_v;
        }
        result_Box.rho[i] = e_rho[i][0];
        if(params.magnetic)
        {
            const double E_mag_B = std::sqrt(e_B[i][0]*e_B[i][0] + e_B[i][1]*e_B[i][1] + e_B[i][2]*e_B[i][2]);
            if(E_mag_B > 0.0)
            {
                const double B_mag = std::sqrt(E_mag_B);
                result_Box.B[i][0] = B_mag*e_B[i][0]/E_mag_B;
                result_Box.B[i][1] = B_mag*e_B[i][1]/E_mag_B;
                result_Box.B[i][2] = B_mag*e_B[i][2]/E_mag_B;
            }
        }
    }
    std::cerr << "Created one turbulent box. \n";
    return result_Box;
}

void Simulation::create_nested_turbulent_boxes()
{
    if(params.LBox.empty())
    {
        logfile << "No turbulence boxes specified\n";
        return;
    }
    const size_t N_gas = gas.x.size();
    turbulence = SampledPositions{};
    turbulence.x = gas.x;
    turbulence.y = gas.y;
    turbulence.z = gas.z;
    turbulence.v.resize(N_gas, {0.0, 0.0, 0.0});
    turbulence.B.resize(N_gas, {0.0, 0.0, 0.0});
    turbulence.rho.resize(N_gas, 0.0);
    turbulence.T.resize(N_gas, 0.0);
    const int NBox = params.NBox;
    for(int i = NBox - 1; i >= 0; i--)
    {
        const double LBox = params.LBox[i];
        const double dxBox = params.dx[i];
        std::cerr << "Start with turbulent box " << i + 1 << "\n";
        logfile << "Turbulence box " << i + 1 << ": LBox = " << LBox/kpc_to_cm << " kpc, dx = " << dxBox/kpc_to_cm << " kpc\n";
        double k_inj;
        if(i == NBox - 1)
        {
            k_inj = 2.0*M_PI/params.L_inj;
        }
        else
        {
            k_inj = 2.0*M_PI/params.dx[i + 1];
        }
        logfile << "k_inj = " << k_inj << "\n";
        SampledPositions box = create_turbulent_box(k_inj, LBox, dxBox);
        SampledPositions box_on_gas;
        box_on_gas.x = gas.x;
        box_on_gas.y = gas.y;
        box_on_gas.z = gas.z;
        box_on_gas.v.resize(N_gas, {0.0, 0.0, 0.0});
        box_on_gas.B.resize(N_gas, {0.0, 0.0, 0.0});
        box_on_gas.rho.resize(N_gas, 0.0);
        box_on_gas.T.resize(N_gas, 0.0);
        interpolate_turbulence(box, box_on_gas);
        save_3D(box_on_gas, "box_" + std::to_string(i) + ".hdf5");
        for(size_t n = 0; n < N_gas; n++)
        {
            turbulence.rho[n] += box_on_gas.rho[n];
            turbulence.v[n][0] += box_on_gas.v[n][0];
            turbulence.v[n][1] += box_on_gas.v[n][1];
            turbulence.v[n][2] += box_on_gas.v[n][2];
            if(params.magnetic)
            {
                turbulence.B[n][0] += box_on_gas.B[n][0];
                turbulence.B[n][1] += box_on_gas.B[n][1];
                turbulence.B[n][2] += box_on_gas.B[n][2];
            }
        }
        save_3D(turbulence, "turbulence_" + std::to_string(i) + ".hdf5");
        logfile << "Added turbulence from box " << i + 1 << "\n";
    }
}

void Simulation::normalize_turbulence()
{
    const size_t N = gas.x.size();

    if (turbulence.x.size() != N || turbulence.y.size() != N || turbulence.z.size() != N || turbulence.rho.size() != N || turbulence.T.size() != N || turbulence.v.size() != N || (params.magnetic && turbulence.B.size() != N))
    {
        throw std::runtime_error("normalize_turbulence: turbulence field sizes do not match gas");
    }

    std::vector<double> thermal_energy_density(N);

    for (size_t i = 0; i < N; i++)
    {
        thermal_energy_density[i] = gas.P[i] / (GAMMA - 1.0);
    }

    const double sigma_rho = calculate_sigma(turbulence.rho);

    std::vector<double> vx(N);
    std::vector<double> vy(N);
    std::vector<double> vz(N);
    std::vector<double> Bx;
    std::vector<double> By;
    std::vector<double> Bz;

    if (params.magnetic)
    {
        Bx.resize(N);
        By.resize(N);
        Bz.resize(N);
    }

    for (size_t i = 0; i < N; i++)
    {
        vx[i] = turbulence.v[i][0];
        vy[i] = turbulence.v[i][1];
        vz[i] = turbulence.v[i][2];

        if (params.magnetic)
        {
            Bx[i] = turbulence.B[i][0];
            By[i] = turbulence.B[i][1];
            Bz[i] = turbulence.B[i][2];
        }
    }

    const double sigma_vx = calculate_sigma(vx);
    const double sigma_vy = calculate_sigma(vy);
    const double sigma_vz = calculate_sigma(vz);

    const double sigma_v = std::sqrt((sigma_vx * sigma_vx + sigma_vy * sigma_vy + sigma_vz * sigma_vz) / 3.0);

    double sigma_B = 0.0;

    if (params.magnetic)
    {
        const double sigma_Bx = calculate_sigma(Bx);
        const double sigma_By = calculate_sigma(By);
        const double sigma_Bz = calculate_sigma(Bz);

        sigma_B = std::sqrt((sigma_Bx * sigma_Bx + sigma_By * sigma_By + sigma_Bz * sigma_Bz) / 3.0);
    }

    logfile << "Density turbulence sigma = " << sigma_rho << "\n";
    logfile << "Velocity turbulence sigma = " << sigma_v << "\n";

    if (params.magnetic)
    {
        logfile << "Magnetic turbulence sigma = " << sigma_B << "\n";
    }

    if (params.turbulence_scaling)
    {
        logfile << "Allowed number of negative cells = " << params.Nneg << "\n";

        const size_t Np = gas.x.size();

        params.turbulence_rho_percentage = estimate_max_turbulence_percentage(Np, params.Nneg);

        logfile << "Turbulence scaling enabled\n" << "Number of cells = " << Np << "\n" << "Maximum turbulence percentage = " << params.turbulence_rho_percentage << "\n";

        params.turbulence_v_percentage = params.turbulence_rho_percentage / 4.0;

        if (params.magnetic)
        {
            params.turbulence_B_percentage = params.turbulence_v_percentage;
        }

        const double percentage = params.turbulence_rho_percentage;

        const double expected_negative = Np * 0.5 * std::erfc(1.0 / (std::sqrt(2.0) * std::sqrt(percentage)));

        logfile << "Expected negative cells = " << expected_negative << "\n";
        logfile << "Velocity turbulence percentage = " << params.turbulence_v_percentage << "\n";

        if (params.magnetic)
        {
            logfile << "Magnetic turbulence percentage = " << params.turbulence_B_percentage << "\n";
        }
    }

    double rho_factor = 0.0;

    if (sigma_rho > 0.0)
    {
        rho_factor = std::sqrt(params.turbulence_rho_percentage / (sigma_rho * sigma_rho));
    }

    double B_sigma_factor = 0.0;

    if (params.magnetic && sigma_B > 0.0)
    {
        B_sigma_factor = std::sqrt(params.turbulence_B_percentage / (sigma_B * sigma_B));
    }

    for (size_t i = 0; i < N; i++)
    {
        if (sigma_rho > 0.0)
        {
            turbulence.rho[i] *= gas.rho[i] * rho_factor;
        }

        if (gas.rho[i] > 0.0)
        {
            turbulence.T[i] = -gas.T[i] * turbulence.rho[i] / gas.rho[i];
        }
        else
        {
            turbulence.T[i] = 0.0;
        }

        if (sigma_v > 0.0)
        {
            const double factor = std::sqrt(2.0 * params.turbulence_v_percentage * thermal_energy_density[i] / (3.0 * gas.rho[i] * sigma_v * sigma_v));

            turbulence.v[i][0] *= factor;
            turbulence.v[i][1] *= factor;
            turbulence.v[i][2] *= factor;
        }

        if (params.magnetic && sigma_B > 0.0)
        {
            const double factor = std::sqrt(8.0 * M_PI * thermal_energy_density[i] / 3.0) * B_sigma_factor;

            turbulence.B[i][0] *= factor;
            turbulence.B[i][1] *= factor;
            turbulence.B[i][2] *= factor;
        }

        const double smooth_factor = ellipsoid_smoothing_factor(turbulence.x[i], turbulence.y[i], turbulence.z[i], params.galaxy.Rgas, params.galaxy.Hgas, params.m_cutoff);

        const double R_outer = params.LBox.back() * 4 / 5 / 2.0;
        const double delta_outer = params.LBox.back() * 1 / 10 / 2.0;

        const double outer_factor = outer_smoothing_factor(turbulence.x[i], turbulence.y[i], turbulence.z[i], R_outer, delta_outer);

        const double total_factor = smooth_factor * outer_factor;

        turbulence.T[i] *= total_factor;
        turbulence.rho[i] *= total_factor;
        turbulence.v[i][0] *= total_factor;
        turbulence.v[i][1] *= total_factor;
        turbulence.v[i][2] *= total_factor;

        if (params.magnetic)
        {
            turbulence.B[i][0] *= total_factor;
            turbulence.B[i][1] *= total_factor;
            turbulence.B[i][2] *= total_factor;
        }
    }
}

void Simulation::add_turbulence_to_gas()
{
    save_3D(gas, "gas_before_adding.hdf5");
    const size_t N = gas.x.size();

    if (N == 0)
    {
        logfile << "No gas cells\n";
        return;
    }

    const double T_floor = params.Tfloor;

    for (size_t i = 0; i < N; i++)
    {
        double T_tot = gas.T[i] + turbulence.T[i];

        if (T_tot < 0.0)
        {
            T_tot = T_floor;
        }

        const double u_tot = k_B * T_tot / ((GAMMA - 1.0) * mu * m_p);

        double rho_tot = gas.rho[i] + turbulence.rho[i];

        if (rho_tot < 0.0)
        {
            rho_tot = gas.P[i] / ((GAMMA - 1.0) * u_tot);
        }

        gas.T[i] = T_tot;
        gas.rho[i] = rho_tot;

        gas.v[i][0] += turbulence.v[i][0];
        gas.v[i][1] += turbulence.v[i][1];
        gas.v[i][2] += turbulence.v[i][2];

        if (params.magnetic)
        {
            gas.B[i][0] += turbulence.B[i][0];
            gas.B[i][1] += turbulence.B[i][1];
            gas.B[i][2] += turbulence.B[i][2];
        }
    }

    save_3D(gas, "gas_after_adding.hdf5");
    logfile << "Added turbulence to gas fields\n";
}

void Simulation::run()
{
    logfile.open("./output/out_" + params.name + ".txt");

    if(!logfile)
    {
        throw std::runtime_error("Could not open logfile");
    }

    std::cerr << "opened logfile \n";

    create_profiles();

    std::cerr << "created 1D profiles \n";

    initialize_conditions();

    std::cerr << "Initialized ICs\n";

    print_initial_conditions();

    integrate();

    std::cerr << "Integrated 1D \n";

    if(M_total.size() != M_CGM.size())
    {
        throw std::runtime_error("M_total and M_CGM have different sizes");
    }

    for(size_t i = 0; i < M_total.size(); ++i)
    {
        M_total[i] += M_CGM[i];
    }

    logfile << "\n===== Simulation 1D finished =====\n";

    convert_1D_to_3D();

    std::cerr << "Converted to 3D \n";

    logfile << "\n===== Simulation 3D finished =====\n";

    if(params.turbulence)
    {
        create_nested_turbulent_boxes();

        std::cerr << "created nested turbulent boxes \n";

        save_3D(turbulence, "turbulence_after_nestedBoxes.hdf5");

        normalize_turbulence();

        std::cerr << "normalized turbulence \n";

        save_3D(turbulence, "turbulence_after_normalization.hdf5");

        add_turbulence_to_gas();

        std::cerr << "added turbulence to gas \n";
    }

    logfile.close();
}

void Simulation::save_1D(std::string filename)
{
    H5::H5File file(filename, H5F_ACC_TRUNC);

    auto write_vector = [&](const std::string& name, const std::vector<double>& data)
    {
        hsize_t dim = data.size();
        H5::DataSpace space(1, &dim);
        H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    };

    size_t end_index = sonic_index;

    while(end_index < radius.size() && radius[end_index] <= Rmax)
    {
        ++end_index;
    }

    auto slice = [&](const std::vector<double>& v)
    {
        return std::vector<double>(v.begin() + sonic_index, v.begin() + end_index);
    };

    write_vector("Radius", slice(radius));
    write_vector("Velocity", slice(velocity));
    write_vector("Temperature", slice(temperature));
    write_vector("Density", slice(density));
    write_vector("Mass", slice(M_total));
    write_vector("MCGM", slice(M_CGM));
    write_vector("Potential", slice(potential));
    write_vector("Metallicity", slice(Z));
    write_vector("Lambda", slice(Lambda_fctn));

    auto write_attribute = [&](const std::string& name, double value)
    {
        H5::DataSpace space(H5S_SCALAR);
        H5::Attribute attr = file.createAttribute(name, H5::PredType::NATIVE_DOUBLE, space);
        attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    };

    write_attribute("M200", params.M200);
    write_attribute("Mcgm_Rmax", Mcgm_Rmax);
    write_attribute("Mcgm_R200", Mcgm_R200);
    write_attribute("R200", halo->R200);
    write_attribute("boxsize", params.boxsize);
    write_attribute("Rmax", Rmax);
    write_attribute("concentration", params.galaxy.concentration);
    write_attribute("fbaryon_tot", (params.galaxy.Mstellar + params.galaxy.Mgas + Mcgm_Rmax) / params.M200);
    write_attribute("fbaryon_R200", (params.galaxy.Mstellar + params.galaxy.Mgas + Mcgm_R200) / params.M200);
    write_attribute("redshift", params.z);
    write_attribute("Mstellar", params.galaxy.Mstellar);
    write_attribute("Mgas", params.galaxy.Mgas);
    write_attribute("Rgas", params.galaxy.Rgas);
    write_attribute("Hgas", params.galaxy.Hgas);
    write_attribute("Rstar", params.galaxy.Rstar);
    write_attribute("Hstar", params.galaxy.Hstar);
    write_attribute("Rsonic", params.galaxy.Rsonic);
    write_attribute("Mdot", params.galaxy.Mdot);
    write_attribute("Mach_init", MACH_INIT);
    write_attribute("epsilon", epsilon);
}


void Simulation::save_3D(
    const SampledPositions& sample,
    const std::string& filename
)
{
    using namespace H5;

    H5File file(filename, H5F_ACC_TRUNC);

    const hsize_t N = sample.x.size();

    if (sample.y.size() != N ||
        sample.z.size() != N ||
        sample.rho.size() != N ||
        sample.v.size() != N)
    {
        throw std::runtime_error(
            "save_3D: sample arrays have inconsistent sizes"
        );
    }

    std::vector<double> P;
    std::vector<double> T;

    const bool has_P = (sample.P.size() == N);
    const bool has_T = (sample.T.size() == N);

    if (!has_P && !has_T)
    {
        throw std::runtime_error(
            "save_3D: neither pressure nor temperature is available"
        );
    }

    if (has_P)
        P = sample.P;

    if (has_T)
        T = sample.T;

    if (!has_P)
    {
        P.resize(N);

        for (size_t i = 0; i < N; i++)
        {
            const double n = sample.rho[i] / (mu * m_p);

            P[i] = n * k_B * T[i];
        }
    }

    if (!has_T)
    {
        T.resize(N);

        for (size_t i = 0; i < N; i++)
        {
            const double n = sample.rho[i] / (mu * m_p);

            if (n > 0.0)
                T[i] = P[i] / (n * k_B);
            else
                T[i] = 0.0;
        }
    }

    // ========================================================
    // Coordinates
    // ========================================================

    std::vector<double> coordinates(3 * N);

    for (size_t i = 0; i < N; i++)
    {
        coordinates[3*i + 0] = sample.x[i];
        coordinates[3*i + 1] = sample.y[i];
        coordinates[3*i + 2] = sample.z[i];
    }

    {
        hsize_t dims[2] = {N, 3};
        DataSpace space(2, dims);

        DataSet dataset = file.createDataSet(
            "Coordinates",
            PredType::NATIVE_DOUBLE,
            space
        );

        dataset.write(
            coordinates.data(),
            PredType::NATIVE_DOUBLE
        );
    }

    // ========================================================
    // Density
    // ========================================================

    {
        hsize_t dims[1] = {N};
        DataSpace space(1, dims);

        DataSet dataset = file.createDataSet(
            "Density",
            PredType::NATIVE_DOUBLE,
            space
        );

        dataset.write(
            sample.rho.data(),
            PredType::NATIVE_DOUBLE
        );
    }

    // ========================================================
    // Pressure
    // ========================================================

    {
        hsize_t dims[1] = {N};
        DataSpace space(1, dims);

        DataSet dataset = file.createDataSet(
            "Pressure",
            PredType::NATIVE_DOUBLE,
            space
        );

        dataset.write(
            P.data(),
            PredType::NATIVE_DOUBLE
        );
    }

    // ========================================================
    // Temperature
    // ========================================================

    {
        hsize_t dims[1] = {N};
        DataSpace space(1, dims);

        DataSet dataset = file.createDataSet(
            "Temperature",
            PredType::NATIVE_DOUBLE,
            space
        );

        dataset.write(
            T.data(),
            PredType::NATIVE_DOUBLE
        );
    }

    // ========================================================
    // Velocity
    // ========================================================

    std::vector<double> velocity(3 * N);

    for (size_t i = 0; i < N; i++)
    {
        velocity[3*i + 0] = sample.v[i][0];
        velocity[3*i + 1] = sample.v[i][1];
        velocity[3*i + 2] = sample.v[i][2];
    }

    {
        hsize_t dims[2] = {N, 3};
        DataSpace space(2, dims);

        DataSet dataset = file.createDataSet(
            "Velocity",
            PredType::NATIVE_DOUBLE,
            space
        );

        dataset.write(
            velocity.data(),
            PredType::NATIVE_DOUBLE
        );
    }

    // ========================================================
    // Magnetic field
    // ========================================================

    if (params.magnetic)
    {
        if (sample.B.size() != N)
        {
            throw std::runtime_error(
                "save_3D: magnetic field size does not match sample"
            );
        }

        std::vector<double> Bfield(3 * N);

        for (size_t i = 0; i < N; i++)
        {
            Bfield[3*i + 0] = sample.B[i][0];
            Bfield[3*i + 1] = sample.B[i][1];
            Bfield[3*i + 2] = sample.B[i][2];
        }

        hsize_t dims[2] = {N, 3};
        DataSpace space(2, dims);

        DataSet dataset = file.createDataSet(
            "MagneticField",
            PredType::NATIVE_DOUBLE,
            space
        );

        dataset.write(
            Bfield.data(),
            PredType::NATIVE_DOUBLE
        );
    }
}


Simulation::Simulation(
    const Parameters& params,
    const std::string& cooling_file
)
:
params(params),
cooling(cooling_file)
{

}