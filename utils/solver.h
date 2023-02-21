#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include "constants.h"
#include "COO.h"

#define ind(x, y, z) ((x)+(y)*Nx+(z)*Nx*Ny)


/*-----------------------------------------------------------------------------
                    S   C   H   E   M   E
         *
         |
         |
        [2]
         |
         |
*-[3]---[0]---[1]-*
         |
         |
        [4]
         |
         |
         *


for i = 0 : n * m
    i = 6
    [1]: 2 * Kx[i] * Kx[i+1] / ((hx^2) * (Kx[i] + Kx[i+1])) = Tau1
    [3]: 2 * Kx[i] * Kx[i-1] / ((hx^2) * (Kx[i] + Kx[i-1])) = Tau3
    [2]: 2 * Ky[i] * Ky[i-Nx] / ((hy^2) * (Ky[i] + Ky[i-Nx])) = Tau2
    [4]: 2 * Ky[i] * Ky[i+Nx] / ((hy^2) * (Ky[i] + Ky[i+Nx])) = Tau4
    [0]: -(Tau1 + Tau2 + Tau3 + Tau4) = Tau0

    A = (Nx * Ny) x (Nx * Ny)

    A[i, i + 1] = Tau1
    A[i, i] = Tau0
    A[i, i - 1] = Tau3
    A[i, i - Nx] = Tau2
    A[i, i + Nx] = Tau4
    b = 0

 ------------------------------------------------------------------------
    if i < Nx:
        [2]: b = -2Ky[i] / (hy**2) * (dirichlet_up)
            A[i,i] += 2Ky[i] / (hy**2)
        [3]: b += -2Kx[i] / (hx**2) * (dirichlet_left)
            A[i,i] += 2Kx[i] / (hx**2)

                    S   C   H   E   M   E
--------------------------------------------------------------------------






*/


COO get_SLAE(
    COO A,
    double b[Nx*Ny],
    std::vector<double> kx,
    std::vector<double> ky,
    std::vector<double> kz,
    std::vector<double> s,
    std::vector<double> phi,
    std::vector<double> p,
    std::vector<double> s_initial,
    std::vector<double> r_o,
    std::vector<double> r_w,
    )
{
    std::cout << "\nStart creating matrix A and vector b" << std::endl;
    auto begin = std::chrono::steady_clock::now();
    double Tau0, Tau1, Tau2, Tau3, Tau4;
    double tpso1, tpso2, tpso3, tpso4;
    double tpsw1, tpsw2, tpsw3, tpsw4;



    for (int i = 0; i < Nx*Ny; ++i) {

        b[i] = 0;

        // Upper boundary
        if (i < Nx){
//             b[i] += -2 * ky[i] / (hy * hy) * dirichlet_up;
            b[i] += Neumann_up * ky[i];
            Tau2 = 2 * ky[i] / (hy * hy);
        }
        else {
            Tau2 = 2 * ky[i] * ky[i-Nx] / ((hy^2) * (ky[i] + ky[i-Nx]));
            A.insert_val(i, i - Nx, Tau2);

        }

        // Left boundary
        if (i % Nx == 0){
//             b[i] += -2 * kx[i] / (hx * hx) * (Neumann_left);
            b[i] += Neumann_left * kx[i];
            Tau3 = 2 * kx[i] / (hx * hx);
        }
        else {
            Tau3 = 2 * kx[i] * kx[i-1] / ((hx^2) * (kx[i] + kx[i-1]));
            A.insert_val(i, i - 1, Tau3);

        }

        // Right boundary
        if ((i + 1) % Nx == 0){
//             b[i] += -2 * kx[i] / (hx * hx) * (dirichlet_right);
            b[i] += Neumann_right * kx[i];
            Tau1 = 2 * kx[i] / (hx * hx);
        }
        else {
            Tau1 = 2 * kx[i] * kx[i+1] / ((hx^2) * (kx[i] + kx[i+1]));
            A.insert_val(i, i + 1, Tau1);

        }

        // Bottom boundary
        if (i >= (Ny - 1) * Nx){
//             b[i] += -2 * ky[i] / (hy * hy) * dirichlet_down;
            b[i] += Neumann_down * ky[i];
            Tau4 = 2 * ky[i] / (hy * hy);
        }
        else {
            Tau4 = 2 * ky[i] * ky[i+Nx] / ((hy^2) * (ky[i] + ky[i+Nx]));
            A.insert_val(i, i+Nx, Tau4);

        }
        Tau0 = Tau1 + Tau2 + Tau3 + Tau4;

        double WI = 0;


        A.insert_val(i, i, -Tau0+WI);

    }
    auto end = std::chrono::steady_clock::now();
    auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    std::cout << "Time for creating matrix A and vector b:\t" << (double) reading_time.count() / 1000 << " s" << std::endl;

    begin = std::chrono::steady_clock::now();
    std::string output_path_A = "../data/A.mtx";
    std::cout << "\nSave matrix A as:\t" << output_path_A << std::endl;
    A.write_to_file(output_path_A);

    std::string output_path_b = "../data/b.txt";
    std::ofstream file;
    std::cout << "Save vector b as:\t" << output_path_b << std::endl;
    file.open(output_path_b);
    file << Nx*Ny << std::endl;
    for (int i = 0; i < Nx*Ny; ++i) {
        file << b[i] << std::endl;
    }
    file.close();
    end = std::chrono::steady_clock::now();
    reading_time =std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    std::cout << "Time for saving matrix A and vector b:\t\t" << (double) reading_time.count() / 1000 << " s" << std::endl;
    return A;
}

