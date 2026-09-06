#include "functions.hpp"
#include "constants.hpp"

#include <H5Cpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <random>

// Cooling table: Ploeckinger&Schaye 2020 https://arxiv.org/pdf/2006.14322
// heating and cooling are divided in 'total primodial' and 'total metal', both need to be added

// If you want to use your own cooling table, you should implement it here. The unit of the lambda function is erg cm^3/s
// You should make sure that the cooling table is suitable for the desired halo mass and redshift (especially CR ionization rates and ISRF ionization rates)
// However, IF a solution can be found, the solution seems to differ not too much between different cooling prescriptions

std::vector<double> read_1D_dataset(H5::H5File& file, const std::string& name)
{
    H5::DataSet dataset = file.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    hsize_t size;
    space.getSimpleExtentDims(&size, nullptr);
    std::vector<double> data(size);
    dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    return data;
}

std::vector<std::string> read_string_dataset(const H5::H5File& file, const std::string& name)
{
    H5::DataSet dataset = file.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    hsize_t size;
    space.getSimpleExtentDims(&size, nullptr);
    H5::StrType type = dataset.getStrType();
    size_t len = type.getSize();
    std::vector<char> buffer(size * len);
    dataset.read(buffer.data(), type);
    std::vector<std::string> result;
    for(size_t i = 0; i < size; i++)
    {
        std::string s(&buffer[i * len], len);
        s.erase(s.find_last_not_of(" \0") + 1);
        result.push_back(s);
    }
    return result;
}



std::vector<double> read_nd_dataset(H5::H5File& file, const std::string& name, std::vector<hsize_t>& shape)
{
    H5::DataSet dataset = file.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    int ndim = space.getSimpleExtentNdims();
    shape.resize(ndim);
    space.getSimpleExtentDims(shape.data(), nullptr);
    size_t total = 1;
    for(auto s : shape) total *= s;
    std::vector<double> data(total);
    dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    return data;
}

int CoolingTable::find_nearest(const std::vector<double>& array, double value)
{
    int index = 0;
    double min_dist = std::abs(array[0] - value);
    for(int i = 1; i < array.size(); i++)
    {
        double dist = std::abs(array[i] - value);
        if(dist < min_dist)
        {
            min_dist = dist;
            index = i;
        }
    }
    return index;
}

// CoolingTable::CoolingTable(
//     const std::string& filename
// )
// {


//     H5::H5File file(
//         filename,
//         H5F_ACC_RDONLY
//     );




//     nH_tab =
//         read_1D_dataset(
//             file,
//             "nH_bins"
//         );


//     // T_tab =
//     //     read_1D_dataset(
//     //         file,
//     //         "temperature"
//     //     );

//     // std::cout << ">>> read T " << std::endl;


//     Z_tab =
//         read_1D_dataset(
//             file,
//             "Z_bins"
//         );



//     U_tab =
//         read_1D_dataset(
//             file,
//             "U_bins"
//         );



//     // z_tab =
//     //     read_1D_dataset(
//     //         file,
//     //         "redshift_bins"
//     //     );

//     z_tab.resize(32);

//     double z_min = 0.0;
//     double z_max = 16.0;
//     int Nz = 32;

//     for (int i = 0; i < Nz; ++i)
//     {
//         double fraction = static_cast<double>(i) / (Nz - 1);

//         z_tab[i] =
//             std::pow(1.0 + z_max, fraction) - 1.0;
//     }



//     cooling =
//         read_nd_dataset(
//             file,
//             "cool",
//             cooling_shape
//         );





//     heating =
//         read_nd_dataset(
//             file,
//             "heat",
//             heating_shape
//         );




//     if(cooling_shape.size() != 4)
//     {
//         throw std::runtime_error(
//             "Cooling table must be 4-dimensional"
//         );
//     }


//     if(heating_shape.size() != 4)
//     {
//         throw std::runtime_error(
//             "Heating table must be 4-dimensional"
//         );
//     }


//     if(
//         cooling_shape[0] != z_tab.size() ||
//         cooling_shape[1] != Z_tab.size() ||
//         cooling_shape[2] != nH_tab.size() ||
//         cooling_shape[3] != U_tab.size()
//     )
//     {
//         throw std::runtime_error(
//             "Cooling table dimensions do not match "
//             "z, Z, nH and U bins"
//         );
//     }


//     if(
//         heating_shape[0] != z_tab.size() ||
//         heating_shape[1] != Z_tab.size() ||
//         heating_shape[2] != nH_tab.size() ||
//         heating_shape[3] != U_tab.size()
//     )
//     {
//         throw std::runtime_error(
//             "Heating table dimensions do not match "
//             "z, Z, nH and U bins"
//         );
//     }


// }

// double CoolingTable::Lambda(
//     double rho,
//     double T,
//     double Z,
//     double z
// )
// {


//     const double u_th =
//         1.5
//         * k_B
//         * T
//         / (mu * m_p);




//     const double nH =
//         rho
//         / (mu * m_p)
//         * 0.76;




//     const int idz =
//         find_nearest(
//             z_tab,
//             z
//         );


//     const int idZ =
//         find_nearest(
//             Z_tab,
//             Z
//         );


//     const int idnH =
//         find_nearest(
//             nH_tab,
//             nH
//         );


//     const int idu =
//         find_nearest(
//             U_tab,
//             u_th
//         );




//     const int NZ =
//         cooling_shape[1];

//     const int NnH =
//         cooling_shape[2];

//     const int NU =
//         cooling_shape[3];




//     auto index =
//         [NZ, NnH, NU](
//             int iz,
//             int iZ,
//             int inH,
//             int iU
//         )
//     {
//         return
//             (((static_cast<size_t>(iz)
//                 * NZ + iZ)
//                 * NnH + inH)
//                 * NU + iU);
//     };




//     const double cool_log =
//         cooling[
//             index(
//                 idz,
//                 idZ,
//                 idnH,
//                 idu
//             )
//         ];


//     const double heat_log =
//         heating[
//             index(
//                 idz,
//                 idZ,
//                 idnH,
//                 idu
//             )
//         ];




//     const double cool =
//         std::pow(
//             10.0,
//             cool_log
//         );


//     const double heat =
//         std::pow(
//             10.0,
//             heat_log
//         );




//     return cool - heat;
// }

// We interpolate the 4D (Temperature, metallicity, hydrogen number density, redshift) to ensure a smooth transition

struct InterpIndex
{
    int i0;
    int i1;
    double f;
};

InterpIndex get_interp_index(const std::vector<double>& table, double x)
{
    if (x <= table.front())
        return {0, 0, 0.0};

    if (x >= table.back())
    {
        int i = table.size() - 1;
        return {i, i, 0.0};
    }

    auto it = std::lower_bound(table.begin(), table.end(), x);
    int i1 = std::distance(table.begin(), it);
    int i0 = i1 - 1;

    double f = (x - table[i0]) / (table[i1] - table[i0]);

    return {i0, i1, f};
}

double interp4D(
    const std::vector<double>& table,
    int i0z, int i1z, double fz,
    int i0T, int i1T, double fT,
    int i0Z, int i1Z, double fZ,
    int i0n, int i1n, double fn,
    int species,
    int NT, int NZ, int Nn, int Ns,
    const std::function<size_t(int,int,int,int,int)>& index)
{
    double result = 0.0;

    for (int az = 0; az <= 1; az++)
    {
        double wz = az ? fz : (1.0 - fz);
        int iz = az ? i1z : i0z;

        for (int aT = 0; aT <= 1; aT++)
        {
            double wT = aT ? fT : (1.0 - fT);
            int iT = aT ? i1T : i0T;

            for (int aZ = 0; aZ <= 1; aZ++)
            {
                double wZ = aZ ? fZ : (1.0 - fZ);
                int iZ = aZ ? i1Z : i0Z;

                for (int an = 0; an <= 1; an++)
                {
                    double wn = an ? fn : (1.0 - fn);
                    int in = an ? i1n : i0n;

                    result += wz * wT * wZ * wn
                            * table[index(iz, iT, iZ, in, species)];
                }
            }
        }
    }

    return result;
}


CoolingTable::CoolingTable(const std::string& filename)
{
    H5::H5File file(filename, H5F_ACC_RDONLY);

    nH_tab = read_1D_dataset(file, "TableBins/DensityBins");
    T_tab  = read_1D_dataset(file, "TableBins/TemperatureBins");
    Z_tab  = read_1D_dataset(file, "TableBins/MetallicityBins");

    z_tab = read_1D_dataset(file, "TableBins/RedshiftBins");
    cooling = read_nd_dataset(file, "Tdep/Cooling", cooling_shape);
    heating = read_nd_dataset(file, "Tdep/Heating", heating_shape);

    auto cool_ids = read_string_dataset(file, "IdentifierCooling");
    auto heat_ids = read_string_dataset(file, "IdentifierHeating");

    for(int i = 0; i < cool_ids.size(); i++)
    {
        if(cool_ids[i].find("TotalPrim") != std::string::npos)
            cool_TotalPrim = i;

        if(cool_ids[i].find("TotalMetal") != std::string::npos)
            cool_TotalMetal = i;
    }

    for(int i = 0; i < heat_ids.size(); i++)
    {
        if(heat_ids[i].find("TotalPrim") != std::string::npos)
            heat_TotalPrim = i;

        if(heat_ids[i].find("TotalMetal") != std::string::npos)
            heat_TotalMetal = i;
    }
}


double CoolingTable::Lambda(double rho, double T, double Z, double z)
{
    if(cool_TotalPrim < 0 || cool_TotalMetal < 0 || heat_TotalPrim < 0 || heat_TotalMetal < 0)
        throw std::runtime_error("Cooling identifiers not found");

    double nH = rho / (mu * m_p) * 0.76;
    double log_nH = std::log10(nH), log_T = std::log10(T), log_Z = std::log10(Z);
    InterpIndex iz = get_interp_index(z_tab, z), in = get_interp_index(nH_tab, log_nH), iT = get_interp_index(T_tab, log_T), iZ = get_interp_index(Z_tab, log_Z);
    int NT = cooling_shape[1], NZ = cooling_shape[2], Nn = cooling_shape[3], Ns = cooling_shape[4];

    auto index = [&](int iz_, int iT_, int iZ_, int in_, int ispecies)
    {
        return ((((iz_ * NT + iT_) * NZ + iZ_) * Nn + in_) * Ns + ispecies);
    };

    auto interpolate = [&](const std::vector<double>& table, int species)
    {
        double result = 0.0;

        for(int az = 0; az <= 1; az++)
        {
            int iz_ = az ? iz.i1 : iz.i0;
            double wz = az ? iz.f : 1.0 - iz.f;

            for(int aT = 0; aT <= 1; aT++)
            {
                int iT_ = aT ? iT.i1 : iT.i0;
                double wT = aT ? iT.f : 1.0 - iT.f;

                for(int aZ = 0; aZ <= 1; aZ++)
                {
                    int iZ_ = aZ ? iZ.i1 : iZ.i0;
                    double wZ = aZ ? iZ.f : 1.0 - iZ.f;

                    for(int an = 0; an <= 1; an++)
                    {
                        int in_ = an ? in.i1 : in.i0;
                        double wn = an ? in.f : 1.0 - in.f;
                        double weight = wz * wT * wZ * wn;
                        result += weight * table[index(iz_, iT_, iZ_, in_, species)];
                    }
                }
            }
        }

        return result;
    };

    double log_cool_prim = interpolate(cooling, cool_TotalPrim), log_cool_metal = interpolate(cooling, cool_TotalMetal);
    double log_heat_prim = interpolate(heating, heat_TotalPrim), log_heat_metal = interpolate(heating, heat_TotalMetal);
    double cool_value = std::pow(10.0, log_cool_prim) + std::pow(10.0, log_cool_metal);
    double heat_value = std::pow(10.0, log_heat_prim) + std::pow(10.0, log_heat_metal);

    return cool_value - heat_value;
}