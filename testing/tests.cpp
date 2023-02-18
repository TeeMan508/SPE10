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

    COO A;
    A.insert_val(0,1,5.);
    A.insert_val(0,0,3.);
    A.insert_val(1,1,6.);
    A.insert_val(3,3,4.);
    A.insert_val(2,3,1.);
    A.print_coo();
    A.print_mat();
//    std::cout<<std::endl<<A.len_mat();
//    A(0,1);
//    A(1,1);
//    A(2,3);
//    A(3,3);

//    std::cout<<A(0,1);
    std::cout<<A(2);

    return 0;
}

