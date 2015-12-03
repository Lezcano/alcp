#ifndef __HENSEL_SUBSETS_HPP
#define __HENSEL_SUBSETS_HPP

#include "zxelem.hpp"
#include "fpxelem.hpp"
#include "types.hpp"


class ord{
public:
  ord() {}
  bool operator() (const DegTag & a, const DegTag &b) const{
	  return a.deg < b.deg;
}};

typedef struct{
	unsigned int deg, tag;
} DegTag;

typedef struct {
	Fpxelem> pol;
	std::vector< std::pair< Fpxelem, unsigned int> > factors;
	std::vector<DegTag> degTag; //Quitar esto, ya no es necesario
	std::vector<unsigned int> sums;
	std::vector<std::set<DegTag, ord> > predecessor;
	unsigned int numOfCases;
} Pri;

class HenselSubsets{
	public:
		HenselSubsets(const Zxelem & poli);

		bool oneMorePrime()const;
		void insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors, const Fpxelem & poli);

		bool bestOption();
		
		bool firstIsIrreducible();
		bool secondIsIrreducible();

		void removeFirstLastOption();
		void removeSecondLastOption();

	private:
		unsigned int howManyPrimes = 2;
		bool firstIrr, secondIrr;
		std::vector<Pri> global;
		std::vector<unsigned int> intersection;
		unsigned int intersectionSize;
		unsigned int sumOfDeg;
		
};

#endif // __HENSEL_SUBSETS_HPP_
