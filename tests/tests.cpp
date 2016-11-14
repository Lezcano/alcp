#include "gtest/gtest.h"

#include <vector>
#include <map>

#include "fpelem.hpp"
#include "zxelem.hpp"
#include "fpxelem.hpp"
#include "modularGCD.hpp"
#include "integerCRA.hpp"
#include "generalPurpose.hpp"
#include "hensel.hpp"

using namespace alcp;

TEST(symForm, randomPolynomial){
    Fp_b f(17);
    std::vector<Fpelem_b> v;
    for(int i = 0; i < 9; ++i)
        v.push_back(f.get(7*i+3));
    Fpxelem_b aux (v);
    Zxelem_b zx (aux);
    Zxelem_b test(std::vector<big_int>({3,-7,0, 7, -3, 4, -6, 1, 8}));

    EXPECT_EQ(test, zx);
}


TEST(moudlarGCD, randomPoly){
    constexpr int n = 3;
    Zxelem_b a[n] = {Zxelem_b(std::vector<big_int>({-360, -171, 145, 25, 1})),
                     Zxelem_b(std::vector<big_int>({-5,2,8,-3,-3,0,1,0,1})),
                     Zxelem_b(std::vector<big_int>({2,2}))};
    Zxelem_b b[n] = {Zxelem_b(std::vector<big_int>({-15,-14,-1,15,14,1})),
                     Zxelem_b(std::vector<big_int>({21,-9,-4,0,5,0,3})),
                     Zxelem_b(std::vector<big_int>({1,2,1}))};
    Zxelem_b res[n] = {Zxelem_b(std::vector<big_int>({15,14,1})),
                       Zxelem_b(std::vector<big_int>({1})),
                       Zxelem_b(std::vector<big_int>({1,1}))};
    for(int i=0;i<n;++i)
        EXPECT_EQ(res[i], modularGCD(a[i], b[i]));
}

TEST(CRA, randomPoly){
    EXPECT_EQ(integerCRA({99,97,95}, {49,-21,-30}) ,-272300);
}

TEST(GCD, zero_zero){
    int a, b;
    EXPECT_EQ(gcd(0,0), 0);
    EXPECT_EQ(eea(0, 0, a, b), 0);
}

TEST(GCD, zero_positive){
    int a, b;
    for(int i=1;i<10;++i){
        EXPECT_EQ(gcd(0,i), i);
        EXPECT_EQ(gcd(i,0), i);
        EXPECT_EQ(eea(0,i,a,b), i);
        EXPECT_EQ(eea(i,0,a,b), i);
    }
}

TEST(GCD, zero_negative){
    int a, b;
    for(int i=-1;i>-10;--i){
        EXPECT_EQ(gcd(0,i), -i);
        EXPECT_EQ(gcd(i,0), -i);
        EXPECT_EQ(eea(0,i,a,b), -i);
        EXPECT_EQ(eea(i,0,a,b), -i);
    }
}

TEST(GCD, random_values){
    int a, b;
    EXPECT_EQ(gcd(42,56), 14);
    EXPECT_EQ(eea(42,56, a, b), 14);
    EXPECT_EQ(gcd(987,1491), 21);
    EXPECT_EQ(eea(987,1491, a, b), 21);
}

TEST(pollardo_factoring, prime_numbers){
    constexpr int n = 3;
    long long num[n] = {3,5,29};
    std::map<long long, short> sol[n] = {
        {{3,1}},
        {{5,1}},
        {{29,1}},
    };
    for(int i=0;i<n;++i)
        EXPECT_EQ(factorInteger(num[i]), sol[i]);
}

TEST(pollardo_factoring, composite_numbers){
    constexpr int n = 5;
    long long num[n] = {4,20,25,100,300};
    std::map<long long, short> sol[n] = {
        {{2,2}},
        {{2,2}, {5,1}},
        {{5,2}},
        {{2,2}, {5,2}},
        {{2,2}, {3,1}, {5,2}}
    };

    for(int i=0;i<n;++i)
        EXPECT_EQ(factorInteger(num[i]), sol[i]);
}

TEST(pollardo_logarithm, composite_numbers){
    long long log;

    // Compute log_2(5) in F_1019
    EXPECT_TRUE(pollardRhoLogarithm(2, 5, 1019, log));
    EXPECT_EQ(fastPowMod(2, log, 1019), 5);
    EXPECT_TRUE(pollardRhoLogarithm(2, 1, 5, log));
    EXPECT_EQ(fastPowMod(2, log, 5), 1);
}

std::vector<std::vector<Zxelem_b>> test_cases_hensel(){
    std::vector<std::vector<Zxelem_b>> tests(4);
    //(x-3) (x+4) (x^2+2) (x+1) (x^2+1) (x^4+x^3+x^2+x+1)
    tests[0] = {Zxelem_b(std::vector<big_int>({-3,1})),
                Zxelem_b(std::vector<big_int>({4,1})),
                Zxelem_b({2,0,1}),
                Zxelem_b({1,1,1,1,1})
    };

    // Cyclotomic polynomials of order 11 and 13
    tests[1] = {Zxelem_b({1,1,1,1,1,1,1,1,1,1,1,1,1}),
                Zxelem_b({1,1,1,1,1,1,1,1,1,1,1})
    };

    //(2x+5)*(6*x^2-10*x+7)
    tests[2] = {Zxelem_b(std::vector<big_int>({5, 2})),
                Zxelem_b({7, -10, 6})
    };

    tests[3] = {Zxelem_b({1,1,1,1,1,1,1,1,1,1,1,1,1}),
                Zxelem_b({1,1,1,1,1,1,1,1,1,1,1}),
                Zxelem_b({1,1,1,1,1,1,1}),
                Zxelem_b({1,1,1,1,1}),
                Zxelem_b({1,1,1}),
                Zxelem_b(std::vector<big_int>({1,1})),
                Zxelem_b(std::vector<big_int>({0,1}))
      };
    return tests;
}

std::vector<std::pair<Zxelem_b, std::size_t>> vectorWithMultiplicity(std::vector<Zxelem_b> v){
    std::vector<std::pair<Zxelem_b, std::size_t>> ret(v.size());
    std::transform(v.begin(), v.end(), ret.begin(),
            [](Zxelem_b p){ return std::make_pair(p, 1);});
    return ret;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Zxelem_b>& v){
    os << "[";
    for(auto e : v)
        os << e << std::endl;
    os << "]" << std::endl;
    return os;
}

void test_hensel_sq_free(int leadingCoeff, bool double_roots = false){
    std::vector<std::vector<Zxelem_b>> pols = test_cases_hensel();
    Zxelem_b z;
    auto greater_fun = [](const std::pair<Zxelem_b, std::size_t>& a,
                          const std::pair<Zxelem_b, std::size_t>& b){
                if(a.first.deg() > b.first.deg()) return true;
                if(a.first.deg() < b.first.deg()) return false;
                for(int i = a.first.deg(); i >= 0; --i){
                    if(a.first[i] > b.first[i]) return true;
                    if(a.first[i] < b.first[i]) return false;
                }
                return false;
            };

    for(auto pp : pols){
        z = 1;
        for(auto p : pp)
            z *= p;

        auto ppmult = vectorWithMultiplicity(pp);

        if(double_roots){
            z*=z;
            for(auto& e : ppmult)
                e.second = 2;
        }


        if(leadingCoeff != 1){
            ppmult.emplace_back(std::make_pair(leadingCoeff, 1));
            z*=leadingCoeff;
        }


        auto sol = factorizationHensel(z);

        std::sort(ppmult.begin(), ppmult.end(), greater_fun);
        std::sort(sol.begin(), sol.end(), greater_fun);

        EXPECT_EQ(ppmult,sol);
    }
}

TEST(hensel_square_free, monic_polynomials){
    test_hensel_sq_free(1);
}

TEST(hensel_square_free, non_monic_polynomials){
    test_hensel_sq_free(12);
}

TEST(hensel_general, monic_polynomials){
    test_hensel_sq_free(1, true);
}

TEST(hensel_general, non_monic_polynomials){
    test_hensel_sq_free(13, true);
}

int main (int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    int returnValue;

    returnValue =  RUN_ALL_TESTS();

    return returnValue;
}
