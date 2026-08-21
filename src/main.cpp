#include "parameters.hpp"
#include "simulation.hpp"
#include "galaxy_catalogue.hpp"
#include <iostream>
#include <filesystem>

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <param_file> <cooling_table>" << std::endl;
        return 1;
    }

    Parameters params = read_parameters(argv[1]);
    std::string cooling_table = argv[2];
    std::unique_ptr<GalaxyCatalogue> catalogue;

    if(params.create_random_galaxy_catalogue)
    {
        catalogue = std::make_unique<GalaxyCatalogue>();
        catalogue->initialize();
        params.galaxy = catalogue->lookup(params.M200, params.z, params.galaxy_variance);
    }

    params.finalize();

    namespace fs = std::filesystem;

    std::string solution_name = "./output/solution1D_" + params.name + ".hdf5";

    Simulation sim(params, cooling_table);
    sim.run();
    sim.save_1D(solution_name);

    std::string command1 = "python ./output/plot_python.py " + solution_name;
    int result1 = std::system(command1.c_str());

    if(result1 != 0)
        std::cerr << "Python plotting script failed\n";

    sim.save_3D(sim.gas, "./output/solution3D_" + params.name + ".hdf5");

    std::string command2 = std::string("python ./output/plot_python3D.py ./output/solution3D_") + params.name + ".hdf5";
    int result2 = std::system(command2.c_str());

    if(result2 != 0)
        std::cerr << "Python plotting script failed\n";

    return 0;
}