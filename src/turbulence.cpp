#include "turbulence.hpp"

#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <H5Cpp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double interpolate_scalar(
    double x,
    double y,
    double z,
    const SampledPositions& source,
    const std::vector<double>& field
)
{
    const size_t total = source.x.size();
    const size_t N = static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(total))));

    if(N < 2 || N * N * N != total)
    {
        throw std::runtime_error("interpolate_scalar: source grid is not a valid N^3 grid");
    }

    if(source.y.size() != total || source.z.size() != total)
    {
        throw std::runtime_error("interpolate_scalar: coordinate arrays have different sizes");
    }

    if(field.size() != total)
    {
        throw std::runtime_error("interpolate_scalar: field size does not match flattened grid");
    }

    std::vector<double> x_axis(N);
    std::vector<double> y_axis(N);
    std::vector<double> z_axis(N);

    for(size_t ix = 0; ix < N; ix++)
    {
        const size_t index = (ix * N) * N;
        x_axis[ix] = source.x[index];
    }

    for(size_t iy = 0; iy < N; iy++)
    {
        const size_t index = iy * N;
        y_axis[iy] = source.y[index];
    }

    for(size_t iz = 0; iz < N; iz++)
    {
        z_axis[iz] = source.z[iz];
    }

    const double xmin = *std::min_element(x_axis.begin(), x_axis.end());
    const double xmax = *std::max_element(x_axis.begin(), x_axis.end());
    const double ymin = *std::min_element(y_axis.begin(), y_axis.end());
    const double ymax = *std::max_element(y_axis.begin(), y_axis.end());
    const double zmin = *std::min_element(z_axis.begin(), z_axis.end());
    const double zmax = *std::max_element(z_axis.begin(), z_axis.end());

    if(x < xmin || x > xmax ||
       y < ymin || y > ymax ||
       z < zmin || z > zmax)
    {
        return 0.0;
    }

    std::sort(x_axis.begin(), x_axis.end());
    std::sort(y_axis.begin(), y_axis.end());
    std::sort(z_axis.begin(), z_axis.end());

    x_axis.erase(std::unique(x_axis.begin(), x_axis.end()), x_axis.end());
    y_axis.erase(std::unique(y_axis.begin(), y_axis.end()), y_axis.end());
    z_axis.erase(std::unique(z_axis.begin(), z_axis.end()), z_axis.end());

    if(x_axis.size() != N ||
       y_axis.size() != N ||
       z_axis.size() != N)
    {
        throw std::runtime_error("interpolate_scalar: coordinate axes contain duplicate values");
    }

    auto get_index = [](double value, const std::vector<double>& axis)
    {
        auto it = std::upper_bound(axis.begin(), axis.end(), value);

        int index = static_cast<int>(std::distance(axis.begin(), it)) - 1;

        index = std::max(0, std::min(index, static_cast<int>(axis.size()) - 2));

        return index;
    };

    const int i = get_index(x, x_axis);
    const int j = get_index(y, y_axis);
    const int k = get_index(z, z_axis);

    const double x1 = x_axis[i];
    const double x2 = x_axis[i + 1];
    const double y1 = y_axis[j];
    const double y2 = y_axis[j + 1];
    const double z1 = z_axis[k];
    const double z2 = z_axis[k + 1];

    const double fx = (x - x1) / (x2 - x1);
    const double fy = (y - y1) / (y2 - y1);
    const double fz = (z - z1) / (z2 - z1);

    auto index = [N](int i, int j, int k)
    {
        return (static_cast<size_t>(i) * N + j) * N + k;
    };

    const double f000 = field[index(i, j, k)];
    const double f100 = field[index(i + 1, j, k)];
    const double f010 = field[index(i, j + 1, k)];
    const double f110 = field[index(i + 1, j + 1, k)];
    const double f001 = field[index(i, j, k + 1)];
    const double f101 = field[index(i + 1, j, k + 1)];
    const double f011 = field[index(i, j + 1, k + 1)];
    const double f111 = field[index(i + 1, j + 1, k + 1)];

    return
          f000 * (1.0 - fx) * (1.0 - fy) * (1.0 - fz)
        + f100 * fx         * (1.0 - fy) * (1.0 - fz)
        + f010 * (1.0 - fx) * fy         * (1.0 - fz)
        + f110 * fx         * fy         * (1.0 - fz)
        + f001 * (1.0 - fx) * (1.0 - fy) * fz
        + f101 * fx         * (1.0 - fy) * fz
        + f011 * (1.0 - fx) * fy         * fz
        + f111 * fx         * fy         * fz;
}

void interpolate_turbulence(
    const SampledPositions& source,
    SampledPositions& target
)
{
    const size_t Ns = source.x.size();
    const size_t Nt = target.x.size();


    if (source.y.size() != Ns || source.z.size() != Ns)
    {
        throw std::runtime_error("interpolate_turbulence: source grid sizes differ");
    }



    if (target.y.size() != Nt || target.z.size() != Nt)
    {
        throw std::runtime_error("interpolate_turbulence: target grid sizes differ");
    }


    if (source.rho.size() != Ns)
        throw std::runtime_error("interpolate_turbulence: source rho size mismatch");

    if (source.T.size() != Ns)
        throw std::runtime_error("interpolate_turbulence: source T size mismatch");

    if (source.v.size() != Ns)
        throw std::runtime_error("interpolate_turbulence: source v size mismatch");

    if (source.B.size() != Ns)
        throw std::runtime_error("interpolate_turbulence: source B size mismatch");



    target.rho.resize(Nt);
    target.T.resize(Nt);
    target.v.resize(Nt, {0.0, 0.0, 0.0});
    target.B.resize(Nt, {0.0, 0.0, 0.0});


    std::vector<double> vx_source(Ns);
    std::vector<double> vy_source(Ns);
    std::vector<double> vz_source(Ns);

    std::vector<double> Bx_source(Ns);
    std::vector<double> By_source(Ns);
    std::vector<double> Bz_source(Ns);

    for (size_t i = 0; i < Ns; i++)
    {
        vx_source[i] = source.v[i][0];
        vy_source[i] = source.v[i][1];
        vz_source[i] = source.v[i][2];

        Bx_source[i] = source.B[i][0];
        By_source[i] = source.B[i][1];
        Bz_source[i] = source.B[i][2];
    }


    for (size_t i = 0; i < Nt; i++)
    {
        const double x = target.x[i];
        const double y = target.y[i];
        const double z = target.z[i];

        target.rho[i] = interpolate_scalar(x, y, z, source, source.rho);
        target.T[i] = interpolate_scalar(x, y, z, source, source.T);

        target.v[i][0] = interpolate_scalar(x, y, z, source, vx_source);
        target.v[i][1] = interpolate_scalar(x, y, z, source, vy_source);
        target.v[i][2] = interpolate_scalar(x, y, z, source, vz_source);

        target.B[i][0] = interpolate_scalar(x, y, z, source, Bx_source);
        target.B[i][1] = interpolate_scalar(x, y, z, source, By_source);
        target.B[i][2] = interpolate_scalar(x, y, z, source, Bz_source);
    }
}


std::vector<std::array<double, 3>> create_turbulence(double L, int N, double k_inj, double exponent, bool dimensions, bool divergence_free)
{
    double A = 1.0;

    if (N <= 0)
        throw std::runtime_error("Turbulence grid size N must be positive");

    if (L <= 0.0)
        throw std::runtime_error("Turbulence box size L must be positive");

    const size_t total = static_cast<size_t>(N) * static_cast<size_t>(N) * static_cast<size_t>(N);

    SampledPositions result;

    result.x.resize(N);
    result.y.resize(N);
    result.z.resize(N);

    for (int i = 0; i < N; i++)
    {
        double coordinate = -L / 2.0 + static_cast<double>(i) * L / (N - 1);

        result.x[i] = coordinate;
        result.y[i] = coordinate;
        result.z[i] = coordinate;
    }

    const double dx = L / static_cast<double>(N);
    const double k_visc = 2.0 * M_PI / dx;
    const double dk = 2.0 * M_PI / L;

    std::vector<double> k(N);

    for (int i = 0; i < N; i++)
    {
        if (i < N / 2)
            k[i] = 2.0 * M_PI * static_cast<double>(i) / L;
        else
            k[i] = 2.0 * M_PI * static_cast<double>(i - N) / L;
    }

    fftw_complex* ehatX = nullptr;
    fftw_complex* ehatY = nullptr;
    fftw_complex* ehatZ = nullptr;

    fftw_complex* realX = nullptr;
    fftw_complex* realY = nullptr;
    fftw_complex* realZ = nullptr;

    if (dimensions)
    {
        ehatX = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));
        ehatY = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));
        ehatZ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));

        realX = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));
        realY = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));
        realZ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));

        if (!ehatX || !ehatY || !ehatZ || !realX || !realY || !realZ)
            throw std::runtime_error("Could not allocate FFT arrays");
    }
    else
    {
        ehatX = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));
        realX = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * total));

        if (!ehatX || !realX)
            throw std::runtime_error("Could not allocate FFT arrays");
    }

    std::mt19937 rng(1234);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int l = 0; l < N; l++)
            {
                const size_t index = (static_cast<size_t>(i) * N + static_cast<size_t>(j)) * N + static_cast<size_t>(l);

                const double kx = k[i];
                const double ky = k[j];
                const double kz = k[l];
                const double kmag = std::sqrt(kx*kx + ky*ky + kz*kz);

                double sigma = 0.0;

                if (kmag >= k_inj && kmag <= k_visc && kmag > 0.0)
                {
                    sigma = A / std::pow(k_inj, exponent - 2.0) * std::pow(kmag, exponent - 2.0);
                }
                else if (kmag <= k_inj && kmag <= k_visc && kmag > 0.0)
                {
                    sigma = A;
                }

                double X1 = std::max(uniform(rng), 1e-15);
                double X2 = std::max(uniform(rng), 1e-15);

                const double theta = 2.0 * M_PI * X1;
                const double r = std::sqrt(-2.0 * sigma * std::log(X2));

                fftw_complex random_field;
                random_field[0] = r * std::cos(theta);
                random_field[1] = r * std::sin(theta);

                if (!dimensions)
                {
                    ehatX[index][0] = random_field[0];
                    ehatX[index][1] = random_field[1];
                }
                else
                {
                    const double thetaX = 2.0 * M_PI * uniform(rng);
                    const double thetaY = 2.0 * M_PI * uniform(rng);
                    const double thetaZ = 2.0 * M_PI * uniform(rng);

                    ehatX[index][0] = random_field[0] * std::cos(thetaX) - random_field[1] * std::sin(thetaX);
                    ehatX[index][1] = random_field[0] * std::sin(thetaX) + random_field[1] * std::cos(thetaX);

                    ehatY[index][0] = random_field[0] * std::cos(thetaY) - random_field[1] * std::sin(thetaY);
                    ehatY[index][1] = random_field[0] * std::sin(thetaY) + random_field[1] * std::cos(thetaY);

                    ehatZ[index][0] = random_field[0] * std::cos(thetaZ) - random_field[1] * std::sin(thetaZ);
                    ehatZ[index][1] = random_field[0] * std::sin(thetaZ) + random_field[1] * std::cos(thetaZ);

                    if (divergence_free)
                    {
                        std::complex<double> ex(ehatX[index][0], ehatX[index][1]);
                        std::complex<double> ey(ehatY[index][0], ehatY[index][1]);
                        std::complex<double> ez(ehatZ[index][0], ehatZ[index][1]);

                        if (kmag > 0.0)
                        {
                            const std::complex<double> dot = ex*kx + ey*ky + ez*kz;

                            ex -= kx * dot / (kmag * kmag);
                            ey -= ky * dot / (kmag * kmag);
                            ez -= kz * dot / (kmag * kmag);
                        }

                        ehatX[index][0] = ex.real();
                        ehatX[index][1] = ex.imag();

                        ehatY[index][0] = ey.real();
                        ehatY[index][1] = ey.imag();

                        ehatZ[index][0] = ez.real();
                        ehatZ[index][1] = ez.imag();
                    }
                }
            }
        }
    }

    fftw_plan planX = nullptr;
    fftw_plan planY = nullptr;
    fftw_plan planZ = nullptr;

    if (dimensions)
    {
        planX = fftw_plan_dft_3d(N, N, N, ehatX, realX, FFTW_BACKWARD, FFTW_ESTIMATE);
        planY = fftw_plan_dft_3d(N, N, N, ehatY, realY, FFTW_BACKWARD, FFTW_ESTIMATE);
        planZ = fftw_plan_dft_3d(N, N, N, ehatZ, realZ, FFTW_BACKWARD, FFTW_ESTIMATE);

        if (!planX || !planY || !planZ)
            throw std::runtime_error("Could not create FFTW plans");

        fftw_execute(planX);
        fftw_execute(planY);
        fftw_execute(planZ);
    }
    else
    {
        planX = fftw_plan_dft_3d(N, N, N, ehatX, realX, FFTW_BACKWARD, FFTW_ESTIMATE);

        if (!planX)
            throw std::runtime_error("Could not create FFTW plan");

        fftw_execute(planX);
    }

    std::vector<std::array<double, 3>> field(total, {0.0, 0.0, 0.0});

    if (!dimensions)
    {
        for (size_t i = 0; i < total; i++)
        {
            field[i][0] = realX[i][0];
            field[i][1] = 0.0;
            field[i][2] = 0.0;
        }

        fftw_destroy_plan(planX);
        fftw_free(ehatX);
        fftw_free(realX);
    }
    else
    {
        for (size_t i = 0; i < total; i++)
        {
            field[i][0] = realX[i][0];
            field[i][1] = realY[i][0];
            field[i][2] = realZ[i][0];
        }

        fftw_destroy_plan(planX);
        fftw_destroy_plan(planY);
        fftw_destroy_plan(planZ);

        fftw_free(ehatX);
        fftw_free(ehatY);
        fftw_free(ehatZ);

        fftw_free(realX);
        fftw_free(realY);
        fftw_free(realZ);
    }

    return field;
}

double calculate_sigma(const std::vector<double>& values)
{
    if (values.empty())
        return 0.0;

    double mean = 0.0;

    for (double x : values)
        mean += x;

    mean /= static_cast<double>(values.size());

    double variance = 0.0;

    for (double x : values)
    {
        const double dx = x - mean;
        variance += dx * dx;
    }

    variance /= static_cast<double>(values.size());

    return std::sqrt(variance);
}

double outer_smoothing_factor(double x, double y, double z, double R_outer, double delta_outer)
{
    const double r = std::sqrt(x*x + y*y + z*z);

    return 0.5 * (1.0 - std::tanh((r - R_outer) / delta_outer));
}

using namespace H5;
