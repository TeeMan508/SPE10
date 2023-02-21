#include <iostream>
#include <vector>
#include <cmath>

#include "../utils/files.h"
#include "../utils/constants.h"
#include "../utils/solver.h"


/* Run all functions */
int main(int argc, char **argv) {
    std::ifstream file_perm, file_phi;
    std::ofstream output;
    std::string filename_perm, filename_phi;

    filename_perm = (argc < 2) ? "../data/spe_perm.dat" : argv[1];
    filename_phi = (argc < 3) ? "../data/spe_phi.dat" : argv[2];
    file_perm.open(filename_perm);
    file_phi.open(filename_phi);

    if (!file_perm or !file_phi) {
        std::cerr << "Error:\tFile couldn't be opened" << std::endl;
        return -1;
    }
    std::cout << "Read file:\t\t\t" << filename_perm << std::endl;
    std::cout << "Read file:\t\t\t" << filename_phi << std::endl;

    std::vector<double> kx, ky, kz, phiArray;
    readData(file_perm, file_phi, kx, ky, kz, phiArray);

#if GET_SEPARATED_MESH
    std::cout << "\nSave separated layer with z = 50" << std::endl;
    std::vector<double> kx_s, ky_s, kz_s, phiArray_s;
    separateData(kx, ky, kz, phiArray, kx_s, ky_s, kz_s, phiArray_s);
#endif

#if CREATE_SEPARATED_MATRIX
    COO A;
    double b[Nx*Ny];
    get_SLAE(A, b, kx_s, ky_s,kz_s);
#endif

#if SAVE_ALL_MESH_AS_VTK
    std::cout << "\nSave all mesh" << std::endl;
    saveToVTK(kx, ky, kz, phiArray);
#endif

    system("python3 ../notebooks/Matrix.py frame_name");


    return 0;
}