#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "boost/multi_array.hpp"
#include <cassert>

#include "bfpstatic2.cpp"

// using namespace std;

int main(){
    boost::random::mt19937 rng;
    struct timeval tv;
    gettimeofday(&tv, 0);
    rng.seed(tv.tv_usec);

    size_t const ydim = 6;
    size_t const xdim = 6;
    size_t const N = ydim * xdim;
    auto grid = new double[ydim][xdim];

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

    vector<double> grid_1d;
    for (int i = 0; i < ydim; i++){
        for (int j = 0; j < xdim; j++){
            grid_1d.push_back(grid[i][j]);
        }
    }

    auto bfp_grid = BFPStatic<int32_t, N>(grid_1d);

    double delta = 1;
    double epsilon = 0.01;

    uint iterations = 0;
    uint max_iterations = 20;

    while (iterations < max_iterations && delta > epsilon){
        iterations++;
        auto temp = bfp_heat_iteration(bfp_grid, ydim, xdim);
        delta = bfp_sum_abs(temp - bfp_grid);
        cout << delta << endl;
        bfp_grid = temp;
    }

    return 0;
}
