#include "bfpdynamic_lazy.cpp"

using namespace std;

template <typename T> void check(const T& A, const T& B){
    assert(A.size() == B.size());
    size_t N = A.size();

    auto Afloat = A.to_float();
    auto Bfloat = T(B).to_float();
    cout << Afloat << endl;
    cout << Bfloat << endl;
    for(int i=0;i<N;i++){
        if (Afloat[0] != Bfloat[0]){
            cout << "Test Failed: " << Afloat[0] << " != " << Bfloat[0] << endl;
            return;
        }
    }
    cout << "Test Passed" << endl;
}

void check_int(int i, int a){
    if (i != a){
        cout << "Test Failed: " << i << " != " << a << endl;
    }else{
        cout << "Test Passed" << endl;
    }
}

int main(){
    // floor_log2
    // cout << "CALC SHIFTS:" << endl;
    // int8_t a = 1;

    // check_int(floor_log2(a), 6);

    // a = 2;
    // check_int(floor_log2(a), 5);

    // a = 3;
    // check_int(floor_log2(a), 5);

    // a = 63;
    // check_int(floor_log2(a), 1);

    // a = 64;
    // check_int(floor_log2(a), 0);

    // a = 65;
    // check_int(floor_log2(a), 0);

    // a = 127;
    // check_int(floor_log2(a), 0);
    
    // a = -1;
    // check_int(floor_log2(a), 6);

    // a = -2;
    // check_int(floor_log2(a), 5);


    // ADDITION
    cout << endl << "ADDITION:" << endl;

	auto a1 = vector<int8_t>{1};
	auto b1 = vector<double>{1};
    BFPDynamic<int8_t> A100{a1, 0};
    BFPDynamic<int8_t> C100{{b1}};
    check(A100, C100);


 //    cout << endl << "SQRT2:" << endl;

	auto a2 = vector<int8_t>{100};
	auto b2 = vector<double>{1};
    BFPDynamic<int8_t> A500{a1, 0};
    BFPDynamic<int8_t> C500{{b1}};
    check(bfp_sqrt2(A500), C500 );



    // BFPDynamic<int8_t> A101{{1}, 0};
    // BFPDynamic<int8_t> C101{{{2}}};
    // check(A101 + A101, C101);

    // BFPDynamic<int8_t> A102{{126}, 0};
    // BFPDynamic<int8_t> C102{{{253}}};
    // check(A102 + A102, C102);

    // BFPDynamic<int8_t> A103{{1}, 0};
    // BFPDynamic<int8_t> B103{{2}, 0};    
    // BFPDynamic<int8_t> C103{{{3}}};
    // check(A103 + A103, C103);


    // SUBTRACTION
    // cout << endl << "SUBTRACTION:" << endl;

    // BFPDynamic<int8_t> A200{{0}, 0};
    // BFPDynamic<int8_t> C200{{{0}}};
    // check(A200 + A200, C200);

    // BFPDynamic<int8_t> A201{{1}, 0};
    // BFPDynamic<int8_t> C201{{{0}}};
    // check(A201 + A201, C201);

    // BFPDynamic<int8_t> A202{{126}, 0};
    // BFPDynamic<int8_t> C202{{{0}}};
    // check(A202 + A202, C202);

    // BFPDynamic<int8_t> A203{{1}, 0};
    // BFPDynamic<int8_t> B203{{2}, 0};    
    // BFPDynamic<int8_t> C203{{{-1}}};
    // check(A203 + A203, C203);


    // // MULTIPLICATION
    // cout << endl << "MULTIPLICATION:" << endl;

    // BFPDynamic<int8_t> A300{{0}, 0};
    // BFPDynamic<int8_t> B300{{1}, 0};    
    // BFPDynamic<int8_t> C300{{{0}}};
    // check(A300 * B300, C300);

    // BFPDynamic<int8_t> A301{{1}, 0};
    // BFPDynamic<int8_t> B301{{2}, 0};    
    // BFPDynamic<int8_t> C301{{{2}}};
    // check(A301 * B301, C301);

    // BFPDynamic<int8_t> A302{{126}, 0};
    // BFPDynamic<int8_t> B302{{126}, 0};
    // BFPDynamic<int8_t> C302{{{15872}}};
    // check(A302 * B302, C302);

    // BFPDynamic<int8_t> A303{{-1}, 0};
    // BFPDynamic<int8_t> B303{{-1}, 0};
    // BFPDynamic<int8_t> C303{{{1}}};
    // check(A303 * B303, C303);


    // // DIVISION
    // cout << endl << "DIVISION:" << endl;

    // BFPDynamic<int8_t> A400{{1}, 0};
    // BFPDynamic<int8_t> B400{{1}, 0};
    // BFPDynamic<int8_t> C400{{{1}}};
    // check(A400 / B400, C400);

    // BFPDynamic<int8_t> A401{{0}, 0};
    // BFPDynamic<int8_t> B401{{1}, 0};
    // BFPDynamic<int8_t> C401{{{0}}};
    // check(A401 / B401, C401);

    // BFPDynamic<int8_t> A402{{100}, 0};
    // BFPDynamic<int8_t> B402{{10}, 0};
    // BFPDynamic<int8_t> C402{{{10}}};
    // check(A402 / B402, C402);

    // BFPDynamic<int8_t> A403{{10}, 0};
    // BFPDynamic<int8_t> B403{{20}, 0};
    // BFPDynamic<int8_t> C403{{{0.5}}};
    // check(A403 / B403, C403);

    // BFPDynamic<int8_t> A404{{10}, 0};
    // BFPDynamic<int8_t> B404{{40}, 0};
    // BFPDynamic<int8_t> C404{{{0.25}}};
    // check(A404 / B404, C404);

    // BFPStatic<int8_t, 2> A405{{10, 20}, 0};
    // BFPStatic<int8_t, 2> B405{{10, 20}, 0};
    // BFPStatic<int8_t, 2> C405{{{1, 1}}};
    // check(A405 / B405, C405);


    // // SQUARE ROOT
    // cout << endl << "SQUARE ROOT:" << endl;

    // BFPDynamic<int8_t> A501{{1}, 0};
    // BFPDynamic<int8_t> C501{{{1}}};
    // check(bfp_sqrt(A501), C501);

    // BFPDynamic<int8_t> A502{{9}, 0};
    // BFPDynamic<int8_t> C502{{{3}}};
    // check(bfp_sqrt(A502), C502);

    // BFPDynamic<int8_t> A503{{100}, 0};
    // BFPDynamic<int8_t> C503{{{10}}};
    // check(bfp_sqrt(A503), C503);

    // BFPDynamic<int8_t> A504{{2}, 0};
    // BFPDynamic<int8_t> C504{{{1.375}}};
    // check(bfp_sqrt(A504), C504);

    // BFPDynamic<int8_t> A505{{2}, 2};
    // BFPDynamic<int8_t> C505{{{2.75}}};
    // check(bfp_sqrt(A505), C505);

    // BFPDynamic<int8_t> A506{{2}, 2};
    // BFPDynamic<int8_t> C506{{11}, -2};
    // check(bfp_sqrt(A506), C506);

    // BFPDynamic<int8_t> A507{{1}, 10};
    // BFPDynamic<int8_t> C507{{32}, 0};
    // check(bfp_sqrt(A507), C507);

    // SQUARE ROOT2


    return 0;
}
