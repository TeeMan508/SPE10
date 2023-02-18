#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include "constants.h"
#include "COO.h"

#define ind(x, y, z) x+y*Nx+z*Nx*Ny


/* Calculate transitivity between 2 neighbouring cells */
double computeTransability(
        int x, int y, int z, std::string axis,
        std::vector<double> &kx,
        std::vector<double> &ky,
        std::vector<double> &kz) {
    double prevK, currK;
    if (axis == "x") {
        prevK = kx[ind(x - 1, y, z)]; //(x - 1) + y * Nx + z * Nx * Ny
        currK = kx[ind(x, y, z)];  //x + y * Nx + z * Nx * Ny
    } else if (axis == "y") {
        prevK = ky[ind(x, y - 1, z)];   //x + (y - 1) * Nx + z * Nx * Ny
        currK = ky[ind(x, y, z)];
    } else if (axis == "z") {
        prevK = kz[ind(x, y, z - 1)];
        currK = kz[ind(x, y, z)];
        std::cerr << "Only 2D examples now!" << std::endl;
        throw;
    } else {
        std::cerr << "Axis Error";
        throw;
    }
    return (2 * prevK * currK / (prevK + currK));
};


std::array<double,3> d2p(double const &T, double const &T_prev, double const &T_next, std::string &axis) {
    std::array<double ,3> sheme;
    double h;
    if (axis == "x") {
        h = hx;
        std::cout << "axis: x" << std::endl;
    } else if (axis == "y") {
        std::cout << "axis: y" << std::endl;
        h = hy;
    } else if (axis == "z") {
        h = hz;
        std::cerr << "Only 2D examples now! we dont know how to do 3d   if u are too smart you can write it by yourself" << std::endl;
        throw;
    } else {
        std::cerr << "Axis Error";
        throw;
    }
    sheme[0] = T_prev / (h * h);
    sheme[1] = -2 * T / (h * h);
    sheme[2] = T_next / (h * h);
    return sheme;
};

std::array<double,Nx*Ny> get_A(){
    std::array<double,Nx*Ny> A;
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nz; ++j) {
            /*TODO*/
        }
    }
}
