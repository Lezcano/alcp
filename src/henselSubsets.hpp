#ifndef __HENSEL_SUBSETS_HPP
#define __HENSEL_SUBSETS_HPP

#include "zxelem.hpp"
#include "fpxelem.hpp"
#include "types.hpp"
#include <set>
#include <vector>
#include <utility>
#include <stack>
#include <map>

struct DegTag{
	unsigned int deg, tag;
};

class ord{
public:
  ord() {}
  bool operator() (const DegTag & a, const DegTag &b) const{
	  return a.tag < b.tag;
}};

struct Option {
	bool b;
	Fpxelem u, w;
};

struct Pri{
	Fpxelem pol;
	std::vector< std::pair< Fpxelem, unsigned int> > factors;
	std::vector<unsigned int> sums;
	std::vector<std::set<DegTag, ord> > predecessor;
	std::map<unsigned int, unsigned int>map;
	unsigned int numOfCases;
};

class HenselSubsets{
	public:
		HenselSubsets(unsigned int poliDeg);

		bool oneMorePrime();
		void insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors, const Fpxelem & poli);

		Option bestOption();
		
		void removeFirstLastOption();

	private:
		unsigned int howManyPrimes;
		std::vector<Pri> global;
		Pri globind;
		std::vector<unsigned int> intersection;
		unsigned int intersectionSize;
		unsigned int semiSumOfDeg;
		int index;
		unsigned int index_intersection;
		std::stack<unsigned int> stackInd;
		std::stack<Fpxelem> stackPol;
		std::stack<std::set<DegTag, ord>::iterator> stackIt;
		unsigned int numOfFactors;
		
};

#endif // __HENSEL_SUBSETS_HPP_
