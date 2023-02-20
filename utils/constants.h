#pragma once

/* Options for run */
#define SAVE_ALL_MESH_AS_VTK true
#define GET_SEPARATED_MESH true
#define CREATE_SEPARATED_MATRIX true

/* Params for problem */
#define Nx (60)
#define Ny (220)
#define Nz (85)
#define hx (1) // Step by x axis
#define hy (1) // Step by y axis
#define hz (1) // Step by z axis
#define dirichlet_up (0)
#define dirichlet_left (0)
#define dirichlet_right (10)
#define dirichlet_down (5)