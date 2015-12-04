//#include "henselSubsets.hpp"
//#include "types.hpp"
//#include <vector>
//#include <set>
//#include <algorithm>
//#include <utility>
//#include <stack>
//#include <map> 
//#include <iostream> //TODO quitar

HenselSubsets::HenselSubsets(const Zxelem & poli){//TODO: De momento solo uso el grado de poli, mirar
	intersectionSize = poli.deg()/2+1;
	semiSumOfDeg = poli.deg()/2;
	intersection.assign(intersectionSize, 1);
	howManyPrimes = 2;
	index = -1;
	index_intersection = -1;
	numOfFactors = 0;
}

bool HenselSubsets::oneMorePrime(){
	return howManyPrimes-- != 0;
}


/* Insertar lo que hace es crear un nuevo elemento en el array "global" al cual le añadimos el polinomio total módulo p, sus factores, un vector sums que en la posición i tiene el número k sii hay k maneras de sumar i combinando los grados de los polinomios de los factores (para calcular esto se hace por programación dinámica sobre el array aux[2]) y un vector de sets en el que la posición i tiene un set. Para cada posible combinación \sum_{j=0}^{k_0} deg(f_j) = i el set guarda el degTag de f_{k_0}.
 *
 * Guardarse este set sirve para reconstruir de manera eficiente las posibles sumas, funciona de la siguiente manera. Para generar la primera suma, supongamos que el vector intersection nos dice que la suma vale i. En ese caso cogeremos el primer elemento del set de i y nos fijamos en qué valor tiene su grado (sea d éste), si el grado es exactamente i entonces ya tenemos una configuración. Si no, nos fijamos en el set de i-d y cogemos el primer elemento y comprobamos si el grado es exactamente i-d y si no lo es procedemos recursivamente. La implementación de esto se hace guardando en una pila cada uno de los elementos de los sets que vamos recorriendo y guardamos también los índices que determinan de qué set eran los elementos. Llevamos también una pila con polinomios en la que la base de la pila tiene el polinomio indicado por el tag de la base de la pila de elementos del set, y en general en una posición cualquiera n de la pila llevamos el producto de la posición anterior multiplicado por el polinomio que determina el tag de la posición n de la pila de elementos de los sets. Esto se hace para generar eficientemente los polinomios que devolvemos en cada iteración (cada iteración será amortizada cte + una division de polinomios). ¿Qué se hace en cada iteración? Se coge el elemento de la pila de los elementos del set y se quita la cima y se sustituye por el siguiente elemento del último set. Si éste se nos ha acabado, hacemos otro pop y cogemos el siguiente elemento del set que nos ha quedado y seguimos hasta reconstruir otra solución.
 * */

//TODO: Por otro lado los sets los tiro a la basura salvo el del bueno, por lo que quizá sea mejor calcularlo una vez que ya sepamos cuál es el bueno


void HenselSubsets::insert(const std::vector<std::pair<Fpxelem, unsigned int> > & factors, const Fpxelem & poli){
	if (index != -1) return; //Si ya he llamado una vez a bestOption no puedo instertar más
	size_t ind = global.size();
	global.push_back({
		poli,
		factors,
		std::vector<unsigned int>(),
		std::vector<std::set<DegTag, ord> >(semiSumOfDeg+1, std::set<DegTag, ord>()),
		std::map<unsigned int, unsigned int>(),
		0}
	); //Quiero guardarme una referencia a factors
	std::vector<unsigned int> aux[2] = {std::vector<unsigned int>(semiSumOfDeg+1, 0), std::vector<unsigned int>(semiSumOfDeg+1, 0)};
	int which = 0;
	unsigned int tag = 0;
	//Voy creando un vector que en la posición i contiene el número de maneras que hay de sumar i con los grados de factors
	for (unsigned int i = 0; i < factors.size(); i++){
		unsigned int deg = factors[i].first.deg();
		for (unsigned int j = 0; j < factors[i].second; j++){
			numOfFactors++;
			for (unsigned int k = 1 ; k <= semiSumOfDeg; k++){//Sólo nos hace falta calcular las cosas hasta la mitad del grado
				if (aux[which][k] != 0 && k + deg <= semiSumOfDeg){
					aux[1-which][k+deg] += aux[which][k];
					global[ind].predecessor[k+deg].insert({deg, tag});
					global[ind].map[tag++]= i;
				}
			}
			aux[1-which][deg]++;
			aux[which] = aux[1-which];
			which = 1 - which;
		}
	} //Ahora tengo en aux[which] un vector con las multiplicidades de las posibles sumas
	
	//Ahora vamos a hacer la intersección con intersection
	for (unsigned int i = 0; i <= semiSumOfDeg; i++){
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

	for (unsigned int i = 1; i<=semiSumOfDeg; i++){
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
		globind = global[index];//move
		global.clear();//TODO: Gestionar bien este destructor. Destruyo todo, es memoria que ya no necesito. En teoría así debería valer.
	}
	if (stackIt.empty()){
		while (++index_intersection <= semiSumOfDeg && intersection[index_intersection] == 0 );
		if (index_intersection == semiSumOfDeg+1) return false;
		stackIt.push(globind.predecessor[index_intersection].begin());
		stackInd.push(index_intersection);
		stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
		while (stackIt.top()->deg != stackInd.top()){//Esto es facil que falle si no se programa bien, si da un error al depurar busca aquí
			stackIt.push(globind.predecessor[stackInd.top() - stackIt.top()->deg].begin());
			stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);
			stackInd.push(stackInd.top() - stackIt.top()->deg);
		} 
		u = stackPol.top();
		w = globind.pol / u;
		return true;
	}
	else{
		while(stackIt.top()++ == globind.predecessor[stackInd.top()].end()){
			stackIt.pop(); stackPol.pop(); stackInd.pop();
		}
		if (stackIt.empty()) return bestOption(u, w);
		auto it = stackIt.top(); 
		stackIt.pop(); stackPol.pop(); stackInd.pop();
		stackIt.push(it++);

		if (stackPol.empty())
			stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
		else
			stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);

		stackInd.push(stackInd.top() - stackIt.top()->deg);
		while (stackIt.top()->deg != stackInd.top()){//Esto es facil que falle si no se programa bien, si da un error al depurar busca aquí
			stackIt.push(globind.predecessor[stackInd.top() - stackIt.top()->deg].begin());
			stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);
			stackInd.push(stackInd.top() - stackIt.top()->deg);
		} 
		u = stackPol.top();
		w = globind.pol / u;
		return true;
	}
}


void HenselSubsets::removeFirstLastOption(){
	while(stackIt.size() > 1){
	//Borrar los elementos usados	

	}
	if (stackIt.top()++ == globind.predecessor[stackInd.top()].end()){
		stackIt.pop(); stackPol.pop(); stackInd.pop();
	}
	else{
		auto aux = stackIt.top();
		stackIt.pop(); stackPol.pop(); stackInd.pop();
		stackIt.push(aux++);
		stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
		stackInd.push(stackInd.top() - stackIt.top()->deg);
	}
}
//NO HAY QUE HACER SPLIT NI RECURSION!!!! EL PROCESAR EL INTERSECTION DE MENOR A MAYOR NOS ASEGURA QUE EL PRIMERO DE LOS DOS FACTORES VA A SER IRREDUCIBLE
