#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <utility>
#include "parameters.hpp"
#include "functions.hpp"
#include "sampled_positions.hpp"
#include "turbulence.hpp"
extern std::ofstream logfile;
#include <vector>
#include <array>

class Simulation
{
public:
    Simulation(const Parameters& params, const std::string& cooling_file);
    void run();
    void save_1D(std::string filename);
    void save_3D(const SampledPositions& sample, const std::string& filename);
    double Mcgm_R200 = 0.0;
    double Mcgm_Rmax = 0.0;
    SampledPositions gas;

private:
    std::vector<double> create_radius_grid(double rmin, double rmax, int N);
    std::vector<double> Lambda_fctn;
    CoolingTable cooling;
    double v_init;
    double T_init;
    double rho_init;
    double Rmax;
    double boxsize;
    int sonic_index;
    double radius_max = 0.0;
    Parameters params;
    void create_profiles();
    void integrate();
    void initialize_conditions();
    void update_initial_state();
    void convert_1D_to_3D();
    SampledPositions create_turbulent_box(double L_inj, double LBox, double dxBox);
    void create_nested_turbulent_boxes();
    void normalize_turbulence();
    void add_turbulence_to_gas();
    SampledPositions turbulence;
    std::vector<double> radius;
    std::vector<double> M_total;
    std::vector<double> potential;
    std::vector<double> Z;
    std::vector<double> velocity;
    std::vector<double> temperature;
    std::vector<double> density;
    std::vector<double> M_CGM;
    std::unique_ptr<HernquistHalo> halo;
    std::unique_ptr<DoubleExponentialDisk> gas_disk;
    std::unique_ptr<DoubleExponentialDisk> stellar_disk;
    std::unique_ptr<HernquistBulge> bulge;
    std::unique_ptr<OuterProfile> outer;
    static constexpr double MACH_INIT = 1.0 - 1e-3;
    static constexpr int max_iter = 100;
    double epsilon = MACH_INIT * MACH_INIT - 1;
    void print_initial_conditions() const;
    void write_catalogue_result(bool reached_Rmax);
    void update_CGM_mass(size_t i, double r_curr, double dr, double rho);
    void compute_CGM_mass();
};