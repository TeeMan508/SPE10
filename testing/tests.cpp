                                            /*  D   R   A   F   T   */

#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/files.h"
#include "../utils/constants.h"
#include "../utils/solver.h"


/* Run all functions */
int main(int argc, char **argv) {
//    double T = 1;
//    double T_prev = 1;
//    double T_next = 1;
//    std::array<double,3> res;
//    std::string axis = "z";
//    res = d2p(T, T_prev, T_next, axis);
//    std::cout<<"res: "<<res[0]<<" "<<res[1]<<" "<<res[2]<<std::endl;

//    COO A;
//    A.insert_val(0,1,5.);
//    A.insert_val(0,0,3.);
//    A.insert_val(1,1,6.);
//    A.insert_val(3,3,4.);
//    A.insert_val(2,3,1.);
//    A.insert_val(0,2,1.5);
//
//    A.print_coo();
//    A.print_mat();
//    std::cout<<std::endl<<A.len_mat();
//    A(0,1);
//    A(1,1);
//    A(2,3);
//    A(3,3);

//    std::cout<<A(0,1);
//    std::cout<<A(2);
//A(5);

    COO A;
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
    std::vector<double> kx_s, ky_s, kz_s;
    separateData(kx, ky, kz, kx_s, ky_s, kz_s);

     A = get_SLAE(kx_s, ky_s, kz_s);
     A.write_to_file();

    return 0;
}

                                            /*  D   R   A   F   T   */