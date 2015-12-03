#ifndef __HENSEL_SUBSETS_HPP
#define __HENSEL_SUBSETS_HPP

#include "zxelem.hpp"
#include "fpxelem.hpp"
#include "types.hpp"
#include <set>
#include <vector>
#include <utility>
#include <stack>

typedef struct{
	unsigned int deg, tag;
} DegTag;

class ord{
public:
  ord() {}
  bool operator() (const DegTag & a, const DegTag &b) const{
	  return a.deg < b.deg;
}};

typedef struct {
	Fpxelem pol;
	std::vector< std::pair< Fpxelem, unsigned int> > factors;
	std::vector<DegTag> degTag; //TODO: Quitar esto, ya no es necesario
	std::vector<unsigned int> sums;
	std::vector<std::multiset<DegTag, ord> > predecessor;
	unsigned int numOfCases;
} Pri;

class HenselSubsets{
	public:
		HenselSubsets(const Zxelem & poli);

		bool oneMorePrime();
		void insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors, const Fpxelem & poli);

		bool bestOption(Fpxelem & u, Fpxelem & w );
		
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
		int index;
		unsigned int index_intersection;
		std::stack<unsigned int> stackInd;
		std::stack<Fpxelem> stackPol;
		std::stack<std::set<DegTag, ord>::iterator> stackIt;
		unsigned int numOfFactors;
		
};

#endif // __HENSEL_SUBSETS_HPP_
