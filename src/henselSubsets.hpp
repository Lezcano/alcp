#ifndef __HENSEL_SUBSETS_HPP
#define __HENSEL_SUBSETS_HPP

#include "zxelem.hpp"
#include "fpxelem.hpp"
#include "types.hpp"


typedef struct{
	unsigned int deg, tag;
} DegTag;

typedef struct {
	std::vector< std::pair< Fpxelem, unsigned int> > factors;
	std::vector<DegTag> degTag;
	std::vector<unsigned int> sums;
	unsigned int numOfCases;
} Pri;

class HenselSubsets{
	public:
		HenselSubsets(const Zxelem & poli);

		bool oneMorePrime()const;
		void insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors);

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
