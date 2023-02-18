#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/files.h"
#include "../utils/constants.h"
#include "../utils/solver.h"


/* Run all functions */
int main(int argc, char **argv) {
    double T = 1;
    double T_prev = 1;
    double T_next = 1;
    std::array<double,3> res;
    std::string axis = "x";
    res = d2p(T, T_prev, T_next, axis);
    std::cout<<"res: "<<res[0]<<" "<<res[1]<<" "<<res[2]<<std::endl;
    return 0;
}