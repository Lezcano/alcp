//#include "henselSubsets.hpp"
//#include "types.hpp"
//#include <vector>
//#include <set>
//#include <algorithm>
//#include <utility>
//#include <stack>
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
	index = -1;
}

bool oneMorePrime(){
	return howManyPrimes-- != 0;
}

void insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors, const Fpxelem & poli){
	size_t ind = global.size();
	global.push_back({poli, factors, std::vector<DegTag>(), std::vector<unsigned int>(), 0}); //Quiero guardarme una referencia a factors
	unsigned int tag = 0, sumOfDeg = 0 ;
	for (unsigned int j = 0; j < factors.size(); j++){
		unsigned int deg = factors[j].first.deg();
		for (unsigned int i = 0; i < factors[j].second; i++){
			global[ind].degTag.push_back({deg, tag});
			numOfFactors++;
		}
		tag++;//TODO IMPORTANTE!! Al poner así el tag los factores los pares deg tag no son necesariamente unicos
	}
	std::vector<unsigned int> aux[2] = {std::vector<unsigned int>(tag, 0), std::vector<unsigned int>(tag, 0)};
	int which = 0;
	for (unsigned int i = 0; i< tag; i++){ //Iterate over the factors
		int deg = global[ind].degTag[i].deg;
		int tag = global[ind].degTag[i].tag;
		for (unsigned int j = 1 ; j < sumOfDeg; j++){//No cuento el caso en el que la suma sea todos, en ese caso no hay que elevar
			if (aux[which][j] != 0){
				aux[1-which][j+deg] += aux[which][j];
				global[ind].predecessor[j+deg].insert({deg, tag});
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
					//global[j].sums[i] = 0; //TODO Esto se puede poner a cero o no, depende de lo que haga después, ahora mismo da igual.
				}
			}
		}
	}
	global[ind].sums = aux[which];
}

bool bestOption(Fpxelem & u, Fpxelem & w ){
	if (index == -1){
		if (global.size == 0) return false;
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
		stackPol.push(global[index].factors[stackIt.top()->tag]);
		while (stackIt.top()->deg != stackInd.top()){//Esto es facil que falle si no se programa bien, si da un error al depurar busca aquí
			stackIt.push(global[index].predecessor[stackInd.top - stackIt.top()->deg].begin());
			stackPol.push(stackPol.top() * global[index].factors[stackIt.top()->tag]);
			stackInd(stackInd.top - stackIt.top()->deg);
		} 
		u = stackPol.top();
		w = global[index].pol / u;
		return true;
	}
	else{
		//Avanzar si se puede
	}
	
	return true;
}

bool firstIsIrreducible(){
	return stackIt.size() == 1;
}
bool secondIsIrreducible(){
	return numOfFactors - stackIt.size() == 1;
}

void removeFirstLastOption(){}
void removeSecondLastOption(){}

int index = -1;
unsigned int index_intersection = -1;
stack<Fpxelem> stackPol;
stack<set<DegTag, ord>::iterator> stackIt;
unsigned int numOfFactors = 0;
