#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/constants.h"

#define ind(x,y,z) x+y*Nx+z*Nx*Ny


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
};


double * d2p(double const& T, double const& T_prev, double const& T_next, std::string &axis){
    double sheme[3];
    double h;
    if (axis=="x") {h=hx;}
    else if (axis=="y"){h=hy;}
    else if (axis=="z"){
        h=hz;
        std::cerr << "Only 2D examples now!" << std::endl;
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

