#ifndef __FP_HPP
#define __HP_HPP

#include "fpelem.hpp"

class Fp{
	Fp(ll p): _p(p){}

	Fpelem getZero(){return Fpelem(0,_p);}
	Fpelem getOne()	{return Fpelem(1,_p);}
	
	vector<Fpelem> getElems(){
		vector<Fpelem> ret;
		for(ll i=0;i<_p;++i) 
			ret.push_back(Fpelem(i,_p));
		return ret;
	}

	ll _p;
};

#endif // __FP_HPP
