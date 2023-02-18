#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <chrono>

#define Nx (60)
#define Ny (220)
#define Nz (85)

#define SAVE_SEPARATED_MESH true
#define SAVE_ALL_MESH true

/* Save all mesh to ["../data/output.vtk" = default] */
void saveToVTK(std::vector <double> & kx, 
               std::vector <double> & ky, 
               std::vector <double> & kz,
               int shapeX = Nx, int shapeY = Ny, int shapeZ = Nz,
               std::string filename = "../data/output.vtk")
{
    std::ofstream output_file;
    auto begin = std::chrono::steady_clock::now();
    std::cout << "Saving as:\t\t\t" << filename << std::endl;
    output_file.open(filename);
    output_file << "# vtk DataFile Version 3.0" << std::endl; 
    output_file << "written by Donskoi Andrei" << std::endl;
    output_file << "ASCII" << std::endl; 
    output_file << "DATASET STRUCTURED_GRID" << std::endl;
    output_file << "DIMENSIONS " << (shapeX+1) << " " << (shapeY+1) << " " << (shapeZ+1) << std::endl;
    output_file << "POINTS " << (shapeX+1) * (shapeY+1) * (shapeZ+1) << " double" << std::endl;

    for (int z = 0; z < shapeZ+1; ++z) 
        for (int y = 0; y < shapeY+1; ++y)
            for (int x = 0; x < shapeX+1; ++x) {
                output_file << x << " " << y << " " << z << std::endl;
            }

    output_file << "POINT_DATA " << (shapeX+1) * (shapeY+1) * (shapeZ+1) << std::endl;

    output_file << "CELL_DATA " << shapeX * shapeY * shapeZ << std::endl;
    output_file << "SCALARS k double 3" << std::endl;
    output_file << "LOOKUP_TABLE default" << std::endl;

    for (int i = 0; i < shapeX * shapeY * shapeZ; ++i) 
        output_file << kx[i] << " " << ky[i] << " " << kz[i] << std::endl;

    auto end = std::chrono::steady_clock::now();
	auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Time for creating VTK-file:\t" << (double) reading_time.count() / 1000 << " s" << std::endl;
    output_file.close();
}


/* Fill kx, ky, kz vectors with values from .dat file */
void readData(std::ifstream & file,
    std::vector <double> & kx, 
    std::vector <double> & ky, 
    std::vector <double> & kz)
{
    auto begin = std::chrono::steady_clock::now();
    try {
        std::string x, y, z;
        while ( kx.size() != Nx*Ny*Nz ) 
        {   
            file >> x;
            kx.push_back(std::stod(x));
        }

        while ( ky.size() != Nx*Ny*Nz ) 
        {   
            file >> y;
            ky.push_back(std::stod(y));
        }

        while ( kz.size() != Nx*Ny*Nz ) 
        {   
            file >> z;
            kz.push_back(std::stod(z));
        }
    } catch (...) {
        std::cerr << "Data Error" << std::endl;
        throw;
    }

    auto end = std::chrono::steady_clock::now();
	auto reading_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Time for reading file:\t\t" << (double) reading_time.count() / 1000 << " s" << std::endl;
    std::cout << "Size of mesh:\t\t\t" << Nx << "x" << Ny << "x" << Nz << std::endl;
    file.close();
}


/* Separate piece of mesh with z = 50 and save to ["../data/sep_output.vtk" = default] */
void separateData(
    std::vector <double> & kx, 
    std::vector <double> & ky, 
    std::vector <double> & kz,
    std::vector <double> & kx_s, 
    std::vector <double> & ky_s, 
    std::vector <double> & kz_s,
    std::string filename = "../data/sep_output.vtk")
{
    const int z = 50;
    for (int y = 0; y < Ny+1; ++y)
        for (int x = 0; x < Nx+1; ++x) {
            kx_s.push_back(kx[x + y * Nx + z * Ny * Nx]);
            ky_s.push_back(ky[x + y * Nx + z * Ny * Nx]);
            kz_s.push_back(kz[x + y * Nx + z * Ny * Nx]);
        }

    saveToVTK(kx_s, ky_s, kz_s, 60, 220, 1, filename);
}


/* Calculate transitivity between 2 neighbouring cells */
double computeTransability(
    int x, int y, int z, std::string axis,
    std::vector <double> & kx,
    std::vector <double> & ky,
    std::vector <double> & kz) 
{
    double prevK, currK;
    if (axis == "x") 
    {
        prevK = kx[(x-1) + y * Nx + z * Nx * Ny];
        currK = kx[x + y * Nx + z * Nx * Ny];
    } 
    else if (axis == "y") 
    {   
        prevK = ky[x + (y-1) * Nx + z * Nx * Ny];
        currK = ky[x + y * Nx + z * Nx * Ny];
    }
    else if (axis == "z") 
    {
        prevK = kz[x + y * Nx + (z-1) * Nx * Ny];
        currK = kz[x + y * Nx + z * Nx * Ny];
        std::cerr << "Only 2D examples now!" << std::endl;
        throw;
    } else {
        std::cerr << "Axis Error";
        throw;
    }
    return 0;
}

/* Save the A-matrix as ["../data/A.mtx" = default]*/
void save2DMatrix(
    std::vector <double> & A,
    int shapeX = Nx, int shapeY = Ny, int shapeZ = Nz,
    std::string filename = "../data/A.mtx"
) 
{
    std::ofstream file;
    file.open(filename);
    
    for (int y = 0; y < shapeX; ++y) // Is it shapeX? mb no?
    {
        for (int x = 0; x < shapeX; ++x) {
            file << A[x] << " ";
        }
        file << std::endl;
    }
    file.close();
}

/* Main algorithm that constructs the A-matrix*/
void createMatrix(
    std::vector <double> & A,
    std::vector <double> & b,
    std::vector <double> & kx,
    std::vector <double> & ky,
    std::vector <double> & kz) 
{
    for (int z = 0; z < Nz+1; ++z) 
        for (int y = 0; y < Ny+1; ++y)
            for (int x = 0; x < Nx+1; ++x) {
                
                if ((x == 0) || (y == 0) || (x == Nx) || (y == Ny))   continue;      // boundary condition (have already been filled with 0)

                double T_xp1_x = computeTransability(x + 1, y, z, "x", kx, ky, kz);  // T (x + 1, x)
                double T_x_xm1 = computeTransability(x, y, z, "x", kx, ky, kz);      // T (x, x - 1)

                double T_yp1_y = computeTransability(x, y + 1, z, "y", kx, ky, kz);  // T (y + 1, y)
                double T_y_ym1 = computeTransability(x, y, z, "y", kx, ky, kz);      // T (y, y - 1)


                /*TODO*/

            }
}

/* Run all functions */
int main(int argc, char** argv) 
{
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
    
    std::vector <double> kx, ky, kz;
    readData(file, kx, ky, kz);

    #if SAVE_SEPARATED_MESH
        std::cout << "\nSave separated layer with z = 50" << std::endl;
        std::vector <double> kx_s, ky_s, kz_s;
        separateData(kx, ky, kz, kx_s, ky_s, kz_s);
    #endif

    #if SAVE_ALL_MESH
        std::cout << "\nSave all mesh" << std::endl;
        saveToVTK(kx, ky, kz);
    #endif

//    std::vector <double> A((Nx+1) * (Ny+1) * (Nz+1) *
//                           (Nx+1) * (Ny+1) * (Nz+1), 0), b((Nx+1) * (Ny+1) * (Nz+1), 0);
//
//    createMatrix(A, b, kx, ky, kz);
    return 0;
}