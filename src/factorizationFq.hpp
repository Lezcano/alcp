#ifndef __FACTORIZATION_FQ
#define __FACTORIZATION_FQ

#include<vector>
#include<utility>

//It returns the factors and their multiplicity
template<typename Fxelem>
std::vector< std::pair< Fxelem, unsigned int> > factorizationBerlekamp (const Fxelem & a);

//It returns the factors and their multiplicity
template<typename Fxelem>
std::vector< std::pair< Fpxelem, unsigned int> > factorizationCantorZassenhaus (const Fpxelem & a);

template<typename Fxelem>
std::vector< Fxelem > berlekamp_simple (const Fxelem &pol);

//Part I
template<typename Fxelem>
std::vector< std::pair< Fxelem, unsigned int> > squareFreeFF (Fxelem a);

//Part II
template<typename Fxelem>
std::vector< std::pair< Fxelem, unsigned int> > partialFactorDD (Fxelem pol);

//Part III
template<typename Fxelem>
std::vector< Fxelem > splitFactorsDD (const Fxelem &pol, int n);

#endif // __FACTORIZATION_FQ
