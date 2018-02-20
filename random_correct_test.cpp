#include <iostream>
#include <cstdlib>
#include <time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "bfpstatic2.cpp"

using namespace std;



int main(){
    // Test af forskellige cases:
    // 1) A og B er disjunkte:
    BFPStatic<int8_t,4> Afp{{{-20.25, 4.5, 29.75, 6.75}}};
    BFPStatic<int8_t,4> Bfp{{{-20.25, 4.5, 29.75, 6.75}}};

    BFPStatic<int8_t,4> A{{3, 18, 119, 27}, 0};
    BFPStatic<int8_t,4> B{{2, -79, 98, -104}, 0};
    BFPStatic<int64_t, 100> ABC = gen_bfp<int64_t, 100>();

    cout << ABC.to_float() << endl;
    // check_add(A,A);
    // check_add(Afp,Afp);
    // check_add(A,B);

    return 0;
}