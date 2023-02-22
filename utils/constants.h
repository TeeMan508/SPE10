#pragma once

/* Options for run */
#define SAVE_ALL_MESH_AS_VTK false
#define GET_SEPARATED_MESH true
#define CREATE_SEPARATED_MATRIX true

/* Params for problem */
#define Nx (60)
#define Ny (220)
#define Nz (85)
#define hx (1) // Step by x axis
#define hy (1) // Step by y axis
#define hz (1) // Step by z axis
#define dt (1) // Step by time

#define dirichlet_up (0)
#define dirichlet_left (0)
#define dirichlet_right (10)
#define dirichlet_down (5)
#define rw (1)
#define S_by_oil (0.75)
#define skinFactor (0)
#define eps (0.01)

#define well1_index (1000)
#define well2_index (Nx*Ny-1000)
#define WellPressure1 (10)
#define WellPressure2 (10)
#define Neumann (1)

#define Neumann_up (0)
#define Neumann_left (0)
#define Neumann_right (0)
#define Neumann_down (0)
