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





CoolingTable::CoolingTable(const std::string& filename)
{
    H5::H5File file(filename, H5F_ACC_RDONLY);
    auto density_log = read_1D_dataset(file, "TableBins/DensityBins");
    for(auto x : density_log) nH_tab.push_back(pow(10.0, x));
    auto temp_log = read_1D_dataset(file, "TableBins/TemperatureBins");
    for(auto x : temp_log) T_tab.push_back(pow(10.0, x));
    auto Z_log = read_1D_dataset(file, "TableBins/MetallicityBins");
    for(auto x : Z_log) Z_tab.push_back(pow(10.0, x));
    z_tab = read_1D_dataset(file, "TableBins/RedshiftBins");
    cooling = read_nd_dataset(file, "Tdep/Cooling", cooling_shape);
    heating = read_nd_dataset(file, "Tdep/Heating", heating_shape);
    auto cool_ids = read_string_dataset(file, "IdentifierCooling");
    auto heat_ids = read_string_dataset(file, "IdentifierHeating");
    for(int i = 0; i < cool_ids.size(); i++)
    {
        if(cool_ids[i].find("TotalPrim") != std::string::npos) cool_TotalPrim = i;
        if(cool_ids[i].find("TotalMetal") != std::string::npos) cool_TotalMetal = i;
    }
    for(int i = 0; i < heat_ids.size(); i++)
    {
        if(heat_ids[i].find("TotalPrim") != std::string::npos) heat_TotalPrim = i;
        if(heat_ids[i].find("TotalMetal") != std::string::npos) heat_TotalMetal = i;
    }
}

double CoolingTable::Lambda(double rho, double T, double Z, double z)
{
    if(cool_TotalPrim < 0 || cool_TotalMetal < 0 || heat_TotalPrim < 0 || heat_TotalMetal < 0)
    {
        throw std::runtime_error("Cooling identifiers not found");
    }
    int idz = find_nearest(z_tab, z);
    double nH = rho / (mu * m_p) * 0.76;
    int idx = find_nearest(nH_tab, nH);
    int idT = find_nearest(T_tab, T);
    int idZ = find_nearest(Z_tab, Z);
    int Nz = cooling_shape[0];
    int NT = cooling_shape[1];
    int NZ = cooling_shape[2];
    int Nn = cooling_shape[3];
    int Ns = cooling_shape[4];
    auto index = [&](int iz, int iT, int iZ, int in, int ispecies)
    {
        return ((((iz * NT + iT) * NZ + iZ) * Nn + in) * Ns + ispecies);
    };
    double cool_value = pow(10.0, cooling[index(idz, idT, idZ, idx, cool_TotalPrim)]) + pow(10.0, cooling[index(idz, idT, idZ, idx, cool_TotalMetal)]);
    double heat_value = pow(10.0, heating[index(idz, idT, idZ, idx, heat_TotalPrim)]) + pow(10.0, heating[index(idz, idT, idZ, idx, heat_TotalMetal)]);
    return cool_value - heat_value;
}