#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/files.h"
#include "../utils/constants.h"
#include "../utils/solver.h"


/* Run all functions */
int main(int argc, char **argv) {
    std::ifstream file;
    std::ofstream output;
    std::string filename;

    filename = (argc < 2) ? "../data/spe_perm.dat" : argv[1];
    file.open(filename);

    if (!file) {
        std::cerr << "Error:\tFile couldn't be opened" << std::endl;
        return -1;
    }
    std::cout << "Read file:\t\t\t" << filename << std::endl;

    std::vector<double> kx, ky, kz;
    readData(file, kx, ky, kz);

#if GET_SEPARATED_MESH
    std::cout << "\nSave separated layer with z = 50" << std::endl;
    std::vector<double> kx_s, ky_s, kz_s;
    separateData(kx, ky, kz, kx_s, ky_s, kz_s);
#endif

#if CREATE_SEPARATED_MATRIX
    COO A;
    double b[Nx*Ny];
    get_SLAE(A, b, kx_s, ky_s, kz_s);
#endif

#if SAVE_ALL_MESH_AS_VTK
    std::cout << "\nSave all mesh" << std::endl;
    saveToVTK(kx, ky, kz);
#endif

    return 0;
}