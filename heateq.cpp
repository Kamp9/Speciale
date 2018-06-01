#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>

#include "boost/multi_array.hpp"
#include <cassert>

#include "bfpstatic2.cpp"

// using namespace std;


template <typename T, size_t N>
void print_grid(BFPStatic<T, N> &A, size_t ydim, size_t xdim){
    for(size_t i=0; i<ydim; i++){
        for(size_t j=0; j<xdim; j++){
            printf("%.7f \t", A[i*ydim+j] * pow(2.0, A.exponent));
        }
        cout << endl;
    }
    cout << endl;
}



int main(){
    boost::random::mt19937 rng;
    struct timeval tv;
    gettimeofday(&tv, 0);
    rng.seed(tv.tv_usec);

    clock_t begin;
    clock_t end;
    double elapsed_secs;

    size_t const ydim = 1000;
    size_t const xdim = 1000;
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

    auto bfp_grid = BFPStatic<int16_t, N>(grid_1d);

    double delta = 1;
    double epsilon = 0.005; //1e-10;

    uint iterations = 0;
    uint max_iterations = 100;

    begin = clock();
    while (iterations < max_iterations && delta > epsilon){
    // while(iterations < max_iterations){
        // print_grid(bfp_grid, ydim, xdim);
        // print_grid(bfp_grid, ydim, xdim);

        iterations++;
        heat_result<int16_t, N> res = bfp_heat_iteration(bfp_grid, ydim, xdim);

        // delta = bfp_sum_abs(temp - bfp_grid, ydim, xdim);
        bfp_grid = res.bfp_grid;
        delta = res.delta;
        cout << delta << endl;
    }
    cout << iterations << endl;
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "elapsed-time: " << elapsed_secs << endl;
    
    return 0;
}
