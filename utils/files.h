#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include "constants.h"

/* Save all mesh to ["../data/output.vtk" = default] */
void saveToVTK(std::vector<double> &kx,
               std::vector<double> &ky,
               std::vector<double> &kz,
               std::vector<double> &phiArray,
               int shapeX = Nx, int shapeY = Ny, int shapeZ = Nz,
               std::string filename = "../data/output.vtk") {
    std::ofstream output_file;
    auto begin = std::chrono::steady_clock::now();
    std::cout << "Saving as:\t\t\t" << filename << std::endl;
    output_file.open(filename);
    output_file << "# vtk DataFile Version 3.0" << std::endl;
    output_file << "written by Donskoi Andrei" << std::endl;
    output_file << "ASCII" << std::endl;
    output_file << "DATASET STRUCTURED_GRID" << std::endl;
    output_file << "DIMENSIONS " << (shapeX + 1) << " " << (shapeY + 1) << " " << (shapeZ + 1) << std::endl;
    output_file << "POINTS " << (shapeX + 1) * (shapeY + 1) * (shapeZ + 1) << " double" << std::endl;

    for (int z = 0; z < shapeZ + 1; ++z)
        for (int y = 0; y < shapeY + 1; ++y)
            for (int x = 0; x < shapeX + 1; ++x) {
                output_file << x << " " << y << " " << z << std::endl;
            }

    output_file << "POINT_DATA " << (shapeX + 1) * (shapeY + 1) * (shapeZ + 1) << std::endl;

    output_file << "CELL_DATA " << shapeX * shapeY * shapeZ << std::endl;
    output_file << "SCALARS k double 3" << std::endl;
    output_file << "LOOKUP_TABLE default" << std::endl;

    for (int i = 0; i < shapeX * shapeY * shapeZ; ++i)
        output_file << kx[i] << " " << ky[i] << " " << kz[i] << std::endl;

    output_file << "POINT_DATA " << (shapeX + 1) * (shapeY + 1) * (shapeZ + 1) << std::endl;
    output_file << "CELL_DATA " << shapeX * shapeY * shapeZ << std::endl;
    output_file << "SCALARS S double 1" << std::endl;
    output_file << "LOOKUP_TABLE default" << std::endl;

    for (int i = 0; i < shapeX * shapeY * shapeZ; ++i)
        output_file << S_by_oil << std::endl;

    output_file << "POINT_DATA " << (shapeX + 1) * (shapeY + 1) * (shapeZ + 1) << std::endl;
    output_file << "CELL_DATA " << shapeX * shapeY * shapeZ << std::endl;
    output_file << "SCALARS Phi double 1" << std::endl;
    output_file << "LOOKUP_TABLE default" << std::endl;

    for (int i = 0; i < shapeX * shapeY * shapeZ; ++i)
        output_file << phiArray[i] << std::endl;

    auto end = std::chrono::steady_clock::now();
    auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Time for creating VTK-file:\t" << (double) reading_time.count() / 1000 << " s" << std::endl;
    output_file.close();
};


/* Separate layer of mesh with z = 50 and save to ["../data/sep_output.vtk" = default] */
void separateData(
        std::vector<double> &kx,
        std::vector<double> &ky,
        std::vector<double> &kz,
        std::vector<double> &phiArray,
        std::vector<double> &kx_s,
        std::vector<double> &ky_s,
        std::vector<double> &kz_s,
        std::vector<double> &phiArray_s,
        std::string filename = "../data/sep_output.vtk") {
    const int z = 50;
    for (int y = 0; y < Ny + 1; ++y)
        for (int x = 0; x < Nx + 1; ++x) {
            kx_s.push_back(kx[x + y * Nx + z * Ny * Nx]);
            ky_s.push_back(ky[x + y * Nx + z * Ny * Nx]);
            kz_s.push_back(kz[x + y * Nx + z * Ny * Nx]);
            phiArray_s.push_back(phiArray[x + y * Nx + z * Ny * Nx]);
        }

    saveToVTK(kx_s, ky_s, kz_s, phiArray_s, 60, 220, 1, filename);
};



/* Fill kx, ky, kz vectors with values from .dat file */
void readData(std::ifstream &file_perm,
              std::ifstream &file_phi,
              std::vector<double> &kx,
              std::vector<double> &ky,
              std::vector<double> &kz,
              std::vector<double> &phiArray
) {
    auto begin = std::chrono::steady_clock::now();
    try {
        std::string x, y, z;
        while (kx.size() != Nx * Ny * Nz) {
            file_perm >> x;
            kx.push_back(std::stod(x));
        }

        while (ky.size() != Nx * Ny * Nz) {
            file_perm >> y;
            ky.push_back(std::stod(y));
        }

        while (kz.size() != Nx * Ny * Nz) {
            file_perm >> z;
            kz.push_back(std::stod(z));
        }
    } catch (...) {
        std::cerr << "Data Error" << std::endl;
        throw;
    }

    try {
        std::string value;
        while (phiArray.size() != Nx * Ny * Nz) {
            file_phi >> value;
            phiArray.push_back(std::stod(value));
        }
    } catch (...) {
        std::cerr << "Data Error" << std::endl;
        throw;
    }

    auto end = std::chrono::steady_clock::now();
    auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Time for reading file:\t\t" << (double) reading_time.count() / 1000 << " s" << std::endl;
    std::cout << "Size of mesh:\t\t\t" << Nx << "x" << Ny << "x" << Nz << std::endl;
    file_perm.close();
    file_phi.close();
};