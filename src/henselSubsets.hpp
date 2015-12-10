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
	Fpxelem_b u, w;
};

struct Pri{
	Fpxelem_b pol;
	std::vector< std::pair< Fpxelem_b, unsigned int> > factors;
	std::vector<unsigned int> sums;
	std::vector<std::set<DegTag, ord> > predecessor;
	std::map<unsigned int, unsigned int>map;
	unsigned int numOfCases;
};

class HenselSubsets{
	public:
		HenselSubsets(const Zxelem_b &poli);

		bool oneMorePrime();
		void insert(const std::vector<std::pair<Fpxelem_b, unsigned int> > & factors, const Fpxelem_b & poli);

		Option bestOption();

		void removeFirstLastOption(Zxelem_b w);

		Zxelem_b getLast();


    private:
        unsigned int intersectionSize;
        unsigned int semiSumOfDeg;
        unsigned int sumOfDeg;
        std::vector<unsigned int> intersection;
        const unsigned int howManyPrimes;
        std::vector<Pri> global;
        int index;
        unsigned int index_intersection;
        unsigned int numOfFactors;
        Pri globind;
        bool hadRemoved;
        Zxelem_b last;
        std::stack<unsigned int> stackInd;
        std::stack<Fpxelem_b> stackPol;
        std::stack<std::set<DegTag, ord>::iterator> stackIt;
};

#endif // __HENSEL_SUBSETS_HPP_
