#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/files.h"
#include "../utils/constants.h"


#define ind(x,y,z) x+y*Nx+z*Nx*Ny

#define SAVE_SEPARATED_MESH true
#define SAVE_ALL_MESH true

/* Calculate transitivity between 2 neighbouring cells */
double computeTransability(
        int x, int y, int z, std::string axis,
        std::vector<double> &kx,
        std::vector<double> &ky,
        std::vector<double> &kz) {
    double prevK, currK;
    if (axis == "x") {
        prevK = kx[ind(x-1,y,z)]; //(x - 1) + y * Nx + z * Nx * Ny
        currK = kx[ind(x,y,z)];  //x + y * Nx + z * Nx * Ny
    } else if (axis == "y") {
        prevK = ky[ind(x,y-1,z)];   //x + (y - 1) * Nx + z * Nx * Ny
        currK = ky[ind(x,y,z)];
    } else if (axis == "z") {
        prevK = kz[ind(x,y,z-1)];
        currK = kz[ind(x,y,z)];
        std::cerr << "Only 2D examples now!" << std::endl;
        throw;
    } else {
        std::cerr << "Axis Error";
        throw;
    }
    return (2*prevK*currK/(prevK+currK));
}



/* Save the A-matrix as ["../data/A.mtx" = default]*/
void save2DMatrix(
        std::vector<double> &A,
        int shapeX = Nx, int shapeY = Ny, int shapeZ = Nz,
        std::string filename = "../data/A.mtx"
) {
    std::ofstream file;
    file.open(filename);

    for (int y = 0; y < shapeX; ++y) // Is it shapeX? mb no?
    {
        for (int x = 0; x < shapeX; ++x) {
            file << A[x] << " ";
        }
        file << std::endl;
    }
    file.close();
}

/* Main algorithm that constructs the A-matrix*/
void createMatrix(
        std::vector<double> &A,
        std::vector<double> &b,
        std::vector<double> &kx,
        std::vector<double> &ky,
        std::vector<double> &kz) {
//    for (int z = 0; z < Nz + 1; ++z)    //<---uncomment this when move on 3d model
        int z=50;  //our layer            //<---  comment this when move on 3d model
        for (int y = 0; y < Ny + 1; ++y)
            for (int x = 0; x < Nx + 1; ++x) {

                if ((x == 0) || (y == 0) || (x == Nx) || (y == Ny))
                    continue;      // boundary condition (have already been filled with 0)

                double T_xp1_x = computeTransability(x + 1, y, z, "x", kx, ky, kz);  // T (x + 1, x)
                double T_x_xm1 = computeTransability(x, y, z, "x", kx, ky, kz);      // T (x, x - 1)

                double T_yp1_y = computeTransability(x, y + 1, z, "y", kx, ky, kz);  // T (y + 1, y)
                double T_y_ym1 = computeTransability(x, y, z, "y", kx, ky, kz);      // T (y, y - 1)


                /*TODO*/

            }
}

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

#if SAVE_SEPARATED_MESH
    std::cout << "\nSave separated layer with z = 50" << std::endl;
    std::vector<double> kx_s, ky_s, kz_s;
    separateData(kx, ky, kz, kx_s, ky_s, kz_s);
#endif

#if SAVE_ALL_MESH
    std::cout << "\nSave all mesh" << std::endl;
    saveToVTK(kx, ky, kz);
#endif

//    std::vector <double> A((Nx+1) * (Ny+1) * (Nz+1) *
//                           (Nx+1) * (Ny+1) * (Nz+1), 0), b((Nx+1) * (Ny+1) * (Nz+1), 0);
//
//    createMatrix(A, b, kx, ky, kz);
    return 0;
}