#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include "constants.h"
#include "COO.h"

#define ind(x, y, z) x+y*Nx+z*Nx*Ny


/*
         *
         |
         |
        [2]
         |
         |
*-[3]---[0]----[1]-*
         |
         |
        [4]
         |
         |
         *


for i=0:n*m
 {

 i=6
    [1]: 2*Kx[i] * Kx[i+1]/((hx^2)*(Kx[i] + Kx[i+1]))=Tau1
    [3]: 2*Kx[i] * Kx[i-1]/((hx^2)*(Kx[i] + Kx[i-1]))=Tau3
    [2]: 2*Ky[i] * Ky[i-Nx]/((hy^2)*(Ky[i] + Ky[i-Nx]))=Tau2
    [4]: 2*Ky[i] * Ky[i+Nx]/((hy^2)*(Ky[i] + Ky[i+Nx]))=Tau4
    [0]: -(Tau1+Tau2+Tau3+Tau4) = Tau0

    A=(Nx*Ny)x(Nx*Ny)

    A[i,i+1] = Tau1
    A[i,i] = Tau0
    A[i,i-1] = Tau3
    A[i,i-Nx] = Tau2
    A[i,i+Nx] = Tau4
    b=0

    if i<Nx:
      [2]: b=-2Ky[i]/(hy**2) * (dirichlet_up)
            A[i,i] += 2Ky[i]/(hy**2)
       [3]: b+=-2Kx[i]/(hx**2) * (dirichlet_left)
            A[i,i] += 2Kx[i]/(hx**2)

 */




COO get_SLAE(std::vector<double> kx,std::vector<double> ky,std::vector<double> kz){
    COO A;
    double Tau1,Tau2,Tau3,Tau4,Tau0,b[Nx*Ny];
    for (int i = 0; i < Nx*Ny; ++i) {

        b[i]=0;
        //if upper boundary
        if (i<Nx){
              b[i]+=-2*ky[i]/(hy*hy)*dirichlet_up;
              Tau2=2*ky[i]/(hy*hy);
        }
        else{
            Tau2 = 2*ky[i] * ky[i-Nx]/((hy^2)*(ky[i] + ky[i-Nx]));
            A.insert_val(i,i-Nx,Tau2);
//            A(i,i-Nx) = Tau2
        }

        //if left boundary
        if (i%Nx==0){
            b[i] += -2 * kx[i] / (hx * hx) * (dirichlet_left);
            Tau3 = 2 * kx[i] / (hx * hx);
        }
        else{
            Tau3 = 2*kx[i] * kx[i-1]/((hx^2)*(kx[i] + kx[i-1]));
            A.insert_val(i,i-1,Tau3);
        }

        //if right boundary
        if ((i+1)%Nx==0){
            b[i] += -2 * kx[i] / (hx * hx) * (dirichlet_right);
            Tau1 = 2 * kx[i] / (hx * hx);
        }
        else{
            Tau1 =2*kx[i] * kx[i+1]/((hx^2)*(kx[i] + kx[i+1]));
            A.insert_val(i,i+1,Tau1);
        }
        //lowet
        if (i>(Ny-1)*Nx){
            b[i]+=-2*ky[i]/(hy*hy)*dirichlet_down;
            Tau4=2*ky[i]/(hy*hy);
        }
        else{
            Tau4 = 2*ky[i] * ky[i+Nx]/((hy^2)*(ky[i] + ky[i+Nx]));
            A.insert_val(i,i+Nx,Tau4);
        }
        Tau0 = Tau1+Tau2+Tau3+Tau4;
        A.insert_val(i,i,-Tau0);
    }
    std::string filename = "../data/b.txt";
    std::ofstream file;
    file.open(filename);
    file<<Nx*Ny<<std::endl;
    for (int i = 0; i < Nx*Ny; ++i) {
        file<<b[i]<<std::endl;
    }
    file.close();
    return A;
}


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

std::array<double,Nx*Ny> get_A(){
    std::array<double,Nx*Ny> A;
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nz; ++j) {
            /*TODO*/
        }
    }
}
