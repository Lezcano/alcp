#include "henselSubsets.hpp"
#include "types.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <utility>
/*
unsigned int howManyPrimes = 2;
bool firstIrr, secondIrr;
std::vector<Pri> global;
std::vector<unsigned int> intersection;
unsigned int intersectionSize;
unsigned int sumOfDeg;*/

HenselSubsets::HenselSubsets(const Zxelem & poli){//TODO: De momento solo uso el grado de poli, mirar
	intersectionSize = poli.deg();
	sumOfDeg = poli.deg();
	intersection.assign(intersectionSize, 1);//TODO: Mirar que esté bien
	firstIrr = 0;
	secondIrr = 0;
	howManyPrimes = 2;
}

bool ord (DegTag a, DegTag b){
	return a.deg < b.deg;
}

bool oneMorePrime(){
	return howManyPrimes-- != 0;
}

void insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors){
	size_t ind = global.size();
	global.push_back({factors, std::vector<DegTag>(), std::vector<unsigned int>(), 0}); //Quiero guardarme una referencia a factors
	unsigned int tag = 0, sumOfDeg = 0 ;
	for (unsigned int j = 0; j < factors.size(); j++){
		unsigned int deg = factors[j].first.deg();
		for (unsigned int i = 0; i < factors[j].second; i++){
			global[ind].degTag.push_back({deg, tag++});
		}
	}
	sort(global[ind].degTag.begin(), global[ind].degTag.end(), ord);//TODO: Mirar si es necesario
	std::vector<unsigned int> aux[2] = {std::vector<unsigned int>(tag, 0), std::vector<unsigned int>(tag, 0)};
	int which = 0;
	for (unsigned int i = 0; i< tag; i++){ //Iterate over the factors
		int deg = global[ind].degTag[i].deg;
		for (unsigned int j = 1 ; j < sumOfDeg; j++){//No cuento el caso en el que la suma sea todos, en ese caso no hay que elevar
			if (aux[which][j] != 0){
				aux[1-which][j+deg] += aux[which][j];
			}
		}
		aux[1-which][deg]++;
		aux[which].assign(tag, 0);
		which = 1 - which;
	} //Ahora tengo en aux[which] un vector con las multiplicidades de las posibles sumas
	
	//Ahora vamos a hacer la intersección con intersection
	for (unsigned int i = 0; i < sumOfDeg; i++){
		if (aux[which][i] != 0){
			if (intersection[i] == 0)
				aux[which][i] = 0;
			else
				global[ind].numOfCases += aux[which][i];
		}
		else{
			if (intersection[i] != 0){
				intersection[i] = 0;
				intersectionSize--;
				for (unsigned int j = 0; j < ind; j++ ){
					global[j].numOfCases -= global[j].sums[i];
					global[j].sums[i] = 0;
				}
			}
		}
	}
	global[ind].sums = aux[which];
}

bool bestOption(Fpxelem ){
	return true;
}

bool firstIsIrreducible(){return firstIrr;}
bool secondIsIrreducible(){return secondIrr;}

void removeFirstLastOption(){}
void removeSecondLastOption(){}


