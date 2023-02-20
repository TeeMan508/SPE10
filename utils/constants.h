#pragma once

/* Options for run */
#define SAVE_SEPARATED_MESH true
#define SAVE_ALL_MESH true
#define CREATE_SEPARATED_MATRIX true

/* Params for problem */
#define Nx (60)
#define Ny (220)
#define Nz (85)
#define hx (1) //step by x axis
#define hy (1) //step by y axis
#define hz (1) //step by z axis
#define dirichlet_up (0)
#define dirichlet_left (0)
#define dirichlet_right (10)
#define dirichlet_down (5)