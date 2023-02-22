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



           |                   |                  |
           |                   |                  |
           |       A1          |       A3         |
           |                   |                  |
           |                   |                  |
           |                   |                  |
 MATRIX =  |-------------------|------------------|
           |                   |                  |
           |                   |                  |
           |         A2        |     A4           |
           |                   |                  |
           |                   |                  |



 now on zero time layer all p = 0


 ---------------------
 func S(s[i],s[j], p[i], p[j])
 if p[i]>p[j] return s[i]
 else return s[j]
 ----------------------
 k(s[i]) = 1

 ---------------------
check = true
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
    -*------------------*---------------------*-----------------*-------------------
    CALCULATING RESIDUAL:
    r_0[i] = phi[i] * (s[i]-s_prev[i])/dt + { Tau1*(p[i]-p[i+1])*S(i,i+1) + Tau2*(p[i]-p[i-nx])*S(i,i-nx)+
        +(Tau3*(p[i]-p[i-1])*S(i,i-1)+(Tau4*(p[i]-p[i+nx])*S(i,i+nx) }
    if i==well_index {r_0[i] -= k(s[i])*WI*(p_well - p[i]) }


    r_w[i] =  phi[i] * (-s[i]+s_initial[i])/dt - (+) { Tau1*(p[i]-p[i+1])*(1-S(i,i+1)) + Tau2*(p[i]-p[i-nx])*(1-S(i,i-nx))+
        +(Tau3*(p[i]-p[i-1])*(1-S(i,i-1))+(Tau4*(p[i]-p[i+nx])*(1-S(i,i+nx)) }
    if i==well_index {r_0[i] -= k(s[i])*WI*(p_well - p[i]) }

    if (abs(r_o[i])+abs(r_w[i]))>eps{
        check = false;
    }
    ---------------------------------------------------------------------
    MAKING J-MATRIX







 ------------------------------------------------------------------------
    if i < Nx:
        [2]: b = -2Ky[i] / (hy**2) * (dirichlet_up)
            A[i,i] += 2Ky[i] / (hy**2)
        [3]: b += -2Kx[i] / (hx**2) * (dirichlet_left)
            A[i,i] += 2Kx[i] / (hx**2)

                    S   C   H   E   M   E
--------------------------------------------------------------------------






*/



double k_o(double p_well, double p_i, double s_i) {
    if (p_well > p_i) { return 0; }
    else { return s_i; }
}

double k_w(double p_well, double p_i, double s_i) {
    if (p_well > p_i) { return 1; }
    else { return 1 - s_i; }
}

double S(double s_i, double s_j, double p_i, double p_j) {
    if (p_i > p_j) { return s_i; }
    else { return s_j; }
}

/* Calculate well index for cell */
double computeWellIndex(
        double kx, double ky
) {
    double value;
    value =
            (2 * M_PI * hz * std::sqrt(kx * ky)) /
            (std::log2(
                    (0.28 / rw) *
                    (std::sqrt(hx * hx * std::sqrt(ky / kx) + hy * hy * std::sqrt(kx / ky))) /
                    (std::pow((ky / kx), 1 / 4) + std::pow((kx / ky), 1 / 4) + skinFactor)
            )
            );
    return value;
}


COO get_SLAE(
        COO A,
        double b[Nx * Ny],
        std::vector<double> kx,
        std::vector<double> ky,
        std::vector<double> kz,
        std::vector<double> s,
        std::vector<double> phi,
        std::vector<double> p,
        std::vector<double> r_o,
        std::vector<double> r_w,
) {
    std::cout << "\nStart creating matrix A and vector b" << std::endl;
    auto begin = std::chrono::steady_clock::now();
    double Tau0_o, Tau1_o, Tau2_o, Tau3_o, Tau4_o;
    double Tau0, Tau1, Tau2, Tau3, Tau4;
    double Tau0_w, Tau1_w, Tau2_w, Tau3_w, Tau4_w;
    double tpso1, tpso2, tpso3, tpso4;
    double tpsw1, tpsw2, tpsw3, tpsw4;
    double dro_dsi0, dro_dsi1,dro_dsi2,dro_dsi3,dro_dsi4;

    double A2i, A2j, A3i, A3j, A4i, A4j;


    bool check_eps = true;
    double t = 0;
    double T = 10;
    double norm_o;
    double norm_w;
    std::vector<double> s_prev = s;
    std::vector<double> p_prev = p;
    while (t < T) {
        norm_o = 0;
        norm_w = 0;
        for (int i = 0; i < Nx * Ny; ++i) {

            Tau0_o = 0;
            Tau0_w = 0;

            b[i] = 0;

            // Upper boundary
            if (i < Nx) {
//             b[i] += -2 * ky[i] / (hy * hy) * Neumann_up;
                b[i] += Neumann_up * ky[i];
                Tau2_o = 0;
                Tau2 = 0;
                Tau2_w = 0;
                tpso2 = 0;
                tpsw2 = 0;
                dro_dsi2=0;
            } else {
                Tau2 = 2 * ky[i] * ky[i - Nx] / ((hy * hx) * (ky[i] + ky[i - Nx]));
                Tau2_o = Tau2 * S(s[i], s[i - Nx], p[i], p[i - Nx]); //fill A1
                Tau2_w = Tau2 * (1 - S(s[i], s[i - Nx], p[i], p[i - Nx])); //fill A2


                dro_dsi2 = Tau2 * (p[i]-p[i-Nx])*(p[i]>p[i-Nx]);


                A.insert_val(i, i - Nx, -Tau2_o);  //fill A1

                A.insert_val(i + Nx * Ny, i - Nx, -Tau2_w);  //fill A2


                tpso2 = Tau2_o * (p[i] - p[i - Nx]) * S(s[i], s[i - Nx], p[i], p[i - Nx]);
                tpsw2 = Tau2_w * (p[i] - p[i - Nx]) * (1 - S(s[i], s[i - Nx], p[i], p[i - Nx]));
            }

            // Left boundary
            if (i % Nx == 0) {
//             b[i] += -2 * kx[i] / (hx * hx) * (Neumann_left);
                b[i] += Neumann_left * kx[i];
                Tau3_o = 0;
                Tau3 = 0;
                Tau3_w = 0;
                tpso3 = 0;
                tpsw3 = 0;
                dro_dsi3=0;

            } else {
                Tau3 = 2 * kx[i] * kx[i - 1] / ((hx * hx) * (kx[i] + kx[i - 1]));
                Tau3_o = Tau3 * S(s[i], s[i - 1], p[i], p[i - 1]);
                Tau3_w = Tau3 * (1 - S(s[i], s[i - 1], p[i], p[i - 1]));

                dro_dsi3 = Tau3 * (p[i]-p[i-1])*(p[i]>p[i-1]);

                A.insert_val(i, i - 1, -Tau3_o); //fill A1
                A.insert_val(i + Nx*Ny, i - 1, -Tau3_w);  //fill A2

                tpso3 = Tau3_o * (p[i] - p[i - 1]) * S(s[i], s[i - 1], p[i], p[i - 1]);
                tpsw3 = Tau3_o * (p[i] - p[i - 1]) * (1 - S(s[i], s[i - 1], p[i], p[i - 1]));
            }

            // Right boundary
            if ((i + 1) % Nx == 0) {
//             b[i] += -2 * kx[i] / (hx * hx) * (Neumann_right);
                b[i] += Neumann_right * kx[i];
                Tau1_o = 0;
                Tau1_w = 0;
                Tau1 = 0;
                tpso1 = 0;
                tpsw1 = 0;
                dro_dsi1=0;
            } else {
                Tau1 = 2 * kx[i] * kx[i + 1] / ((hx * hx) * (kx[i] + kx[i + 1]));
                Tau1_o = Tau1 * S(s[i], s[i + 1], p[i], p[i + 1]);
                Tau1_w = Tau1 * (1 - S(s[i], s[i + 1], p[i], p[i + 1]));

                dro_dsi1 = Tau1 * (p[i]-p[i+1])*(p[i]>p[i+1]);

                A.insert_val(i, i + 1, -Tau1_o);   //fill A1
                A.insert_val(i+Nx*Ny, i + 1, -Tau1_w);  //fill A2

                tpso1 = Tau1_o * (p[i] - p[i + 1]) * S(s[i], s[i + 1], p[i], p[i + 1]);
                tpsw1 = Tau1_o * (p[i] - p[i + 1]) * (1 - S(s[i], s[i + 1], p[i], p[i + 1]));
            }

            // Bottom boundary
            if (i >= (Ny - 1) * Nx) {
//             b[i] += -2 * ky[i] / (hy * hy) * Neumann_down;
                b[i] += Neumann_down * ky[i];
                Tau4_o = 0;
                Tau4_w = 0;
                Tau4 = 0;
                tpso4 = 0;
                tpsw4 = 0;
                dro_dsi4=0;
            } else {
                Tau4 = 2 * ky[i] * ky[i + Nx] / ((hy * hy) * (ky[i] + ky[i + Nx]));
                Tau4_o = Tau4 * S(s[i], s[i + Nx], p[i], p[i + Nx]);
                Tau4_w = Tau4 * (1-S(s[i], s[i + Nx], p[i], p[i + Nx]));


                dro_dsi4 = Tau4 * (p[i]-p[i+Nx])*(p[i]>p[i+Nx]);

                A.insert_val(i, i + Nx, -Tau4_o);  //fill A1
                A.insert_val(i+Nx*Ny, i + Nx, -Tau4_w);  //fill A2

                tpso4 = Tau4_o * (p[i] - p[i + Nx]) * S(s, i, i + Nx, p[i], p[i + Nx]);
                tpsw4 = Tau4_o * (p[i] - p[i + Nx]) * (1 - S(s, i, i + Nx, p[i], p[i + Nx]));
            }

            dro_dsi0 = dro_dsi1+dro_dsi2+dro_dsi3+dro_dsi4 + phi[i]/dt;

            double WI = 0;
            if (i == well1_index) {
                WI = computeWellIndex(kx[i], ky[i]);
//            b[i] += WI * WellPressure1;
                r_o[i] -= k_o(WellPressure1, p[i], s[i]) * WI * (WellPressure1 - p[i]);
                r_w[i] -= k_w(WellPressure1, p[i], s[i]) * WI * (WellPressure1 - p[i]);
                Tau0_o = k_o(WellPressure1, p[i], s[i]) * WI;
                Tau0_w = k_w(WellPressure1, p[i], s[i]) * WI;
            }

            if (i == well2_index) {
                WI = computeWellIndex(kx[i], ky[i]);

                Tau0_o = k_o(WellPressure2, p[i], s[i]) * WI;
                Tau0_w = k_w(WellPressure2, p[i], s[i]) * WI;
//            b[i] += WI * WellPressure2;
                r_o[i] -= Tau0_o * (WellPressure2 - p[i]);
                r_w[i] -= Tau0_w * (WellPressure2 - p[i]);

            }
            Tau0_o += Tau1_o + Tau2_o + Tau3_o + Tau4_o;
            Tau0_w += Tau1_w + Tau2_w + Tau3_w + Tau4_w;

            A.insert_val(i, i, Tau0_o);  //fill A1
            A.insert_val(i+Nx*Ny, i, Tau0_o);  //fill A2

            r_o[i] = phi[i] * (s[i] - s_prev[i]) / dt - (tpso1 + tpso2 + tpso3 + tpso4);
            r_w[i] = phi[i] * (-s[i] + s_prev[i]) / dt - (tpsw1 + tpsw2 + tpsw3 + tpsw4);

            norm_o += r_o * r_o;
            norm_w += r_w * r_w;
        }
        norm_o = sqrt(norm_o) + sqrt(norm_w);
        if (norm_o < eps) {
            s_prev = s;
            p_prev = p;
        } else {
            //making JACOBIAN (JOPOGIKAN)



        }
    }
    auto end = std::chrono::steady_clock::now();
    auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    std::cout << "Time for creating matrix A and vector b:\t" << (double) reading_time.count() / 1000 << " s"
              << std::endl;

    begin = std::chrono::steady_clock::now();
    std::string output_path_A = "../data/A.mtx";
    std::cout << "\nSave matrix A as:\t" << output_path_A << std::endl;
    A.write_to_file(output_path_A);

    std::string output_path_b = "../data/b.txt";
    std::ofstream file;
    std::cout << "Save vector b as:\t" << output_path_b << std::endl;
    file.open(output_path_b);
    file << Nx * Ny << std::endl;
    for (int i = 0; i < Nx * Ny; ++i) {
        file << b[i] << std::endl;
    }
    file.close();
    end = std::chrono::steady_clock::now();
    reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    std::cout << "Time for saving matrix A and vector b:\t\t" << (double) reading_time.count() / 1000 << " s"
              << std::endl;
    return A;
}

