#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpstatic2.cpp"

// using namespace std;

int main(){
    boost::random::mt19937 rng;
    struct timeval tv;
    gettimeofday(&tv, 0);
    rng.seed(tv.tv_usec);

    size_t const ydim = 6;
    size_t const xdim = 6;
    size_t const N = (ydim - 2) * (xdim - 2);
    auto grid = new double[ydim][xdim];
    // #pragma omp parallel for collapse(2)
    for (int i=0; i<ydim; i++) {      // Initialize the grid
        for (int j=0;j<xdim;j++) {
            grid[i][j] = 0;
        }
    }
    for (int i=0; i<ydim; i++) {      // And borders
        grid[i][0]      = -273.15;
        grid[i][xdim-1] = -273.15;
    }
    for (int i=0; i<xdim; i++) {
        grid[0][i]      = -273.15;
        grid[ydim-1][i] = 40.0;
    }

    BFPStatic<int32_t, N> north;
    BFPStatic<int32_t, N> south;
    BFPStatic<int32_t, N> center;
    BFPStatic<int32_t, N> west;
    BFPStatic<int32_t, N> east;

    for (int i=1; i<ydim-1; i++) {
        for(int j=1;j<xdim-1;j++) {
            north[(i-1)*(ydim-2)+(j-1)]  = grid[i-1][j];
            south[(i-1)*(ydim-2)+(j-1)]  = grid[i+1][j]; 
            center[(i-1)*(ydim-2)+(j-1)] = grid[i][j]; 
            west[(i-1)*(ydim-2)+(j-1)]   = grid[i][j-1]; 
            east[(i-1)*(ydim-2)+(j-1)]   = grid[i][j+1];
        }
    }

    uint iterations = 0;
    uint max_iterations = 10;
    while (iterations < max_iterations){
        iterations++;
        BFPStatic<int32_t, N> temp = bfp_mul_scalar((north + south + center + west + east), 0.2);
        center = temp;
        // cout << center << endl;
    }
    vector<double> res_grid = center.to_float();
    cout << res_grid << endl;

    cout << "output: " << endl;;
    for (int i=0; i<(ydim-2); ++i) {
        for (int j=0; j<(xdim-2); ++j) {
            cout << res_grid[i*4+j] << "\t, ";
        }
        cout << endl;
    }
    cout << "]";

    return 0;
}
