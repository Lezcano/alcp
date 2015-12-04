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
	std::vector<unsigned int> sums;
	std::vector<std::set<DegTag, ord> > predecessor;
	std::map<unsigned int, unsigned int>map,
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
		


		Pri globind;
		std::vector<unsigned int> intersection;
		unsigned int intersectionSize;
		unsigned int semiSumOfDeg;//hay que actualizarlo al hacer split
		int index;
		unsigned int index_intersection;
		//Las pilas tienen que ser nuevas. De hecho, no hay que olvidarse de que al hacer split solo queremos usar los elementos de la interseccion que son menores que la mitad del grado actual
		std::stack<unsigned int> stackInd;
		std::stack<Fpxelem> stackPol;
		std::stack<std::set<DegTag, ord>::iterator> stackIt;
		unsigned int numOfFactors;

};

#endif // __HENSEL_SUBSETS_HPP_
