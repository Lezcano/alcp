#include "henselSubsets.hpp"
#include "types.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <utility>
#include <stack>
#include <iostream> //TODO quitar
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
	howManyPrimes = 2;
	index = -1;
	index_intersection = -1;
	numOfFactors = 0;
}

bool HenselSubsets::oneMorePrime(){
	return howManyPrimes-- != 0;
}
//TODO: considerar sólo hasta la mitad del grado del polinomio
void HenselSubsets::insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors, const Fpxelem & poli){
	size_t ind = global.size();
	global.push_back({
		poli,
		factors,
		std::vector<unsigned int>(),
		std::vector<std::multiset<DegTag, ord> >(sumOfDeg+1, std::multiset<DegTag, ord>()),
		0}
	); //Quiero guardarme una referencia a factors
	std::vector<unsigned int> aux[2] = {std::vector<unsigned int>(sumOfDeg+1, 0), std::vector<unsigned int>(sumOfDeg+1, 0)};
	int which = 0;
	//Voy creando un vector que en la posición i contiene el número de maneras que hay de sumar i con los grados de factors
	for (unsigned int i = 0; i < factors.size(); i++){
		unsigned int deg = factors[i].first.deg();
		for (unsigned int j = 0; j < factors[i].second; j++){
			numOfFactors++;
			for (unsigned int k = 1 ; k <= sumOfDeg; k++){//Creo que es necesario el caso de sumOfDeg por si el polinomio es irreducible
				if (aux[which][k] != 0){
					aux[1-which][k+deg] += aux[which][k];
					global[ind].predecessor[k+deg].insert({deg, i});
				}
			}
			aux[1-which][deg]++;
			aux[which] = aux[1-which];
			which = 1 - which;
		}
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
					//global[j].sums[i] = 0; //TODO Esto se puede poner a cero o no, depende de lo que haga después, ahora mismo da igual.
				}
			}
		}
	}
	global[ind].sums = aux[which];

	for (unsigned int i = 1; i<=sumOfDeg; i++){
		if (intersection[i] != 0) std::cout << i << " ";
	}
	std::cout << std::endl;
}

bool HenselSubsets::bestOption(Fpxelem & u, Fpxelem & w ){
	if (index == -1){
		if (global.size() == 0) return false;
		unsigned int min = global[0].numOfCases;
		index = 0; 
		for (unsigned int i = 1; i < global.size(); i++) {
			if (min > global[i].numOfCases){
				min = global[i].numOfCases;
				index = i;
			}
		}
	}
	if (stackIt.empty()){
		while (++index_intersection < sumOfDeg && intersection[index_intersection] == 0 );
		if (index_intersection == sumOfDeg) return false;
		stackIt.push(global[index].predecessor[index_intersection].begin());
		stackInd.push(index_intersection);
		stackPol.push(global[index].factors[stackIt.top()->tag].first);
		while (stackIt.top()->deg != stackInd.top()){//Esto es facil que falle si no se programa bien, si da un error al depurar busca aquí
			stackIt.push(global[index].predecessor[stackInd.top() - stackIt.top()->deg].begin());
			stackPol.push(stackPol.top() * global[index].factors[stackIt.top()->tag].first);
			stackInd.push(stackInd.top() - stackIt.top()->deg);
		} 
		u = stackPol.top();
		w = global[index].pol / u;
		return true;
	}
	else{
		while(stackIt.top()++ == global[index].predecessor[stackInd.top()].end()){
			stackIt.pop(); stackPol.pop(); stackInd.pop();
		}
		if (stackIt.empty()) return bestOption(u, w);
		auto it = stackIt.top(); 
		stackIt.pop(); stackPol.pop(); stackInd.pop();
		stackIt.push(it++);
		stackPol.push(stackPol.top() * global[index].factors[stackIt.top()->tag].first);
		stackInd.push(stackInd.top() - stackIt.top()->deg);
		while (stackIt.top()->deg != stackInd.top()){//Esto es facil que falle si no se programa bien, si da un error al depurar busca aquí
			stackIt.push(global[index].predecessor[stackInd.top() - stackIt.top()->deg].begin());
			stackPol.push(stackPol.top() * global[index].factors[stackIt.top()->tag].first);
			stackInd.push(stackInd.top() - stackIt.top()->deg);
		} 
		u = stackPol.top();
		w = global[index].pol / u;
		return true;
	}
}

bool HenselSubsets::firstIsIrreducible(){
	return stackIt.size() == 1;
}
bool HenselSubsets::secondIsIrreducible(){
	return numOfFactors - stackIt.size() == 1;
}

void HenselSubsets::removeFirstLastOption(){}
void HenselSubsets::removeSecondLastOption(){}
