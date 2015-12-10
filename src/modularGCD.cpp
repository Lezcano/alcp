#include "modularGCD.hpp"
#include "generalPurpose.hpp"
#include "integerCRA.hpp"
#include "zxelem.hpp"
#include "types.hpp"
#include <algorithm> // std::min
#include <chrono>
#include <random>
#include <limits>


big_int intContent (const Zxelem_b &a){
	//Calcula el gcd de los coeficientes con el signo del coeficiente director
}

big_int randomPrime (){
	std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> distr(1, std::numeric_limits<int>::max()/2-2);
    int a;
    do{
        a = distr(generator);
        a = 2*a+1;
    }while(!millerRabin(a));
    return a;
}

Fpxelem_b redModP(const Zxelem_b &a, const Fp_b &f){
    std::vector<Fpelem_b> ret(a.deg()+1);
    for(int i = 0; i <= a.deg(); ++i)
        ret[i] = f.get(a[i]);
    return ret;
}


/* Modular GCD
 *  Given A, B \in Z[X] nonzero, it obtains gcd (A, B) via modular
 *  reduction.
 * */
Zxelem_b modularGCD(Zxelem_b a, Zxelem_b b){
	big_int ia = content(a); a /= ia;
	big_int ib = content(b); b /= ib;
	//Compute coefficient bound of gcd(a, b)
	big_int ic = gcd (ia, ib);
	big_int g = gcd (a.lc(), b.lc());
	big_int q = 0;
    Zxelem_b h = 0;
    Zxelem_b c;
    big_int p;
    std::size_t n = std::min(a.deg(), b.deg());
	big_int limit = (1<<n)*g*std::min(normInf(a), normInf(b));

	while (true){
        do{
            p = randomPrime();
        }while(g % p == 0);

        // These variables ought be defined inside the loop since they are p-dependent
        Fp_b f(p);

        Fpxelem_b ap = redModP(a,f);
        Fpxelem_b bp = redModP(b,f);
        Fpxelem_b cp = gcd(ap, bp);
        if(cp.deg() == 0)
            cp = f.get(1);
        Fpelem_b gp = f.get(g);

        // Normalize so gp = lcoeff(cp)
        cp = gp*cp.lc().inv()*cp;

        // It detects if the previous reductions were unlucky
        if(cp.deg() < n) {
            q = p;
            h = Zxelem_b(cp);
            n = cp.deg();
        }
        else{ // cp.deg() >= n
            for(int i = 0; i <= h.deg(); ++i){
                h[i] = integerCRA({q,p}, {h[i],(big_int)cp[i]});
            }
            q*=p;
        }
        if(q > limit){
            c = h/content(h);
            if(a % c == 0 && b % c == 0)
                return ic*c;
        }
        else if(cp.deg() == 0)
            return ic;
	}
}
