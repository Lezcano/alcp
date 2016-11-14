#ifndef __HENSEL_SUBSETS_HPP
#define __HENSEL_SUBSETS_HPP

#include <set>
#include <vector>
#include <utility>
#include <stack>
#include <map>

#include "zxelem.hpp"
#include "fpxelem.hpp"
#include "types.hpp"

namespace alcp {
	struct DegTag {
		std::size_t deg, tag;
	};

	class ord {
	public:
		ord() { }

		bool operator()(const DegTag &a, const DegTag &b) const {
			return a.tag < b.tag;
		}
	};

	struct Option {
		bool b;
		Fpxelem_b u, w;
	};

	struct Pri {
		Fpxelem_b pol;
		std::vector<std::pair<Fpxelem_b, std::size_t> > factors;
		std::vector<std::size_t> sums;
		std::vector<std::set<DegTag, ord> > predecessor;
		std::map<std::size_t, std::size_t> map;
		std::size_t numOfCases;
	};

	class HenselSubsets {
	public:
		HenselSubsets(const Zxelem_b &poli);

		bool oneMorePrime();

		void insert(const std::vector<std::pair<Fpxelem_b, std::size_t> > &factors, const Fpxelem_b &poli);

		Option bestOption();

		void removeFirstLastOption(Zxelem_b w);

		Zxelem_b getLast();

	private:
		std::size_t semiSumOfDeg;
		std::size_t sumOfDeg;
		std::vector<std::size_t> intersection;
		const std::size_t howManyPrimes;
		std::vector<Pri> global;
		int index;
		std::size_t index_intersection;
		Pri globind;
		bool hadRemoved;
		Zxelem_b last;
		std::stack<std::size_t> stackInd;
		std::stack<Fpxelem_b> stackPol;
		std::stack<std::set<DegTag, ord>::iterator> stackIt;
	};
}

#endif // __HENSEL_SUBSETS_HPP_
