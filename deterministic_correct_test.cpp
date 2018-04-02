#include "bfpstatic2.cpp"

using namespace std;

template <typename T> void check(const T& A, const T& B){
    assert(A.size() == B.size());
    size_t N = A.size();

    auto Afloat = A.to_float(); 
    auto Bfloat = T(B).to_float();
    for(int i=0;i<N;i++){
        if (Afloat[0] != Bfloat[0]){
            cout << "Test Failed: " << Afloat[0] << " != " << Bfloat[0] << endl;
            return;
        }
    }
    cout << "Test Passed" << endl;
}


int main(){
    // ADDITION
    cout << "ADDITION:" << endl;
    BFPStatic<int8_t, 1> A100{{1}, 0};
    BFPStatic<int8_t, 1> C100{{{2}}};
    check(A100 + A100, C100);

    BFPStatic<int8_t, 1> A101{{1}, 0};
    BFPStatic<int8_t, 1> C101{{{2}}};
    check(A101 + A101, C101);

    BFPStatic<int8_t, 1> A102{{126}, 0};
    BFPStatic<int8_t, 1> C102{{{253}}};
    check(A102 + A102, C102);

    BFPStatic<int8_t, 1> A103{{1}, 0};
    BFPStatic<int8_t, 1> B103{{2}, 0};    
    BFPStatic<int8_t, 1> C103{{{3}}};
    check(A103 + A103, C103);


    // SUBTRACTION
    cout << endl << "SUBTRACTION:" << endl;

    BFPStatic<int8_t, 1> A200{{0}, 0};
    BFPStatic<int8_t, 1> C200{{{0}}};
    check(A200 + A200, C200);

    BFPStatic<int8_t, 1> A201{{1}, 0};
    BFPStatic<int8_t, 1> C201{{{0}}};
    check(A201 + A201, C201);

    BFPStatic<int8_t, 1> A202{{126}, 0};
    BFPStatic<int8_t, 1> C202{{{0}}};
    check(A202 + A202, C202);

    BFPStatic<int8_t, 1> A203{{1}, 0};
    BFPStatic<int8_t, 1> B203{{2}, 0};    
    BFPStatic<int8_t, 1> C203{{{-1}}};
    check(A203 + A203, C203);


    // MULTIPLICATION
    cout << endl << "MULTIPLICATION:" << endl;

    BFPStatic<int8_t, 1> A300{{0}, 0};
    BFPStatic<int8_t, 1> B300{{1}, 0};    
    BFPStatic<int8_t, 1> C300{{{0}}};
    check(A300 * B300, C300);

    BFPStatic<int8_t, 1> A301{{1}, 0};
    BFPStatic<int8_t, 1> B301{{2}, 0};    
    BFPStatic<int8_t, 1> C301{{{2}}};
    check(A301 * B301, C301);

    BFPStatic<int8_t, 1> A302{{126}, 0};
    BFPStatic<int8_t, 1> B302{{126}, 0};
    BFPStatic<int8_t, 1> C302{{{15872}}};
    check(A302 * B302, C302);

    BFPStatic<int8_t, 1> A303{{-1}, 0};
    BFPStatic<int8_t, 1> B303{{-1}, 0};
    BFPStatic<int8_t, 1> C303{{{1}}};
    check(A303 * B303, C303);


    // DIVISION
    cout << endl << "DIVISION:" << endl;

    BFPStatic<int8_t, 1> A400{{1}, 0};
    BFPStatic<int8_t, 1> B400{{1}, 0};
    BFPStatic<int8_t, 1> C400{{{1}}};
    check(A400 / B400, C400);

    BFPStatic<int8_t, 1> A401{{0}, 0};
    BFPStatic<int8_t, 1> B401{{1}, 0};
    BFPStatic<int8_t, 1> C401{{{0}}};
    check(A401 / B401, C401);

    BFPStatic<int8_t, 1> A402{{100}, 0};
    BFPStatic<int8_t, 1> B402{{10}, 0};
    BFPStatic<int8_t, 1> C402{{{10}}};
    check(A402 / B402, C402);

    BFPStatic<int8_t, 1> A403{{10}, 0};
    BFPStatic<int8_t, 1> B403{{20}, 0};
    BFPStatic<int8_t, 1> C403{{{0.5}}};
    check(A403 / B403, C403);

    BFPStatic<int8_t, 1> A404{{10}, 0};
    BFPStatic<int8_t, 1> B404{{40}, 0};
    BFPStatic<int8_t, 1> C404{{{0.25}}};
    check(A404 / B404, C404);

    BFPStatic<int8_t, 2> A405{{10, 20}, 0};
    BFPStatic<int8_t, 2> B405{{10, 20}, 0};
    BFPStatic<int8_t, 2> C405{{{1, 1}}};
    check(A405 / B405, C405);


    // SQUARE ROOT
    cout << endl << "SQUARE ROOT:" << endl;

    BFPStatic<int8_t, 1> A501{{1}, 0};
    BFPStatic<int8_t, 1> C501{{{1}}};
    check(bfp_sqrt(A501), C501);

    BFPStatic<int8_t, 1> A502{{9}, 0};
    BFPStatic<int8_t, 1> C502{{{3}}};
    check(bfp_sqrt(A502), C502);

    BFPStatic<int8_t, 1> A503{{100}, 0};
    BFPStatic<int8_t, 1> C503{{{10}}};
    check(bfp_sqrt(A503), C503);

    BFPStatic<int8_t, 1> A504{{2}, 0};
    BFPStatic<int8_t, 1> C504{{{1.375}}};
    check(bfp_sqrt(A504), C504);

    BFPStatic<int8_t, 1> A505{{2}, 2};
    BFPStatic<int8_t, 1> C505{{{2.75}}};
    check(bfp_sqrt(A505), C505);

    BFPStatic<int8_t, 1> A506{{2}, 2};
    BFPStatic<int8_t, 1> C506{{11}, -2};
    check(bfp_sqrt(A506), C506);

    BFPStatic<int8_t, 1> A507{{1}, 10};
    BFPStatic<int8_t, 1> C507{{32}, 0};
    check(bfp_sqrt(A507), C507);

    BFPStatic<int8_t, 1> A508{{2}, 10};
    BFPStatic<int8_t, 1> C508{{2}, 10};
    check(bfp_sqrt(A508), BFPStatic<int8_t, 1>(Vsqrt(C508)));

    return 0;
}
