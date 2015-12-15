#include "henselSubsets.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <utility> // std::move
#include <stack>
#include <map>
#include <iostream> //TODO quitar

#include "types.hpp"

const bool verbose = false;

//TODO: De momento solo uso el grado de poli, mirar
HenselSubsets::HenselSubsets(const Zxelem_b &poli):
    intersectionSize (poli.deg()/2+1),
    semiSumOfDeg (poli.deg()/2),
    sumOfDeg (poli.deg()),
    intersection(intersectionSize, 1),
    howManyPrimes (6),
    index (-1),
    index_intersection (-1),
    numOfFactors (0),
    hadRemoved (false),
    last(poli) { };


bool HenselSubsets::oneMorePrime(){
    return global.size() != howManyPrimes;
}


/* Insertar lo que hace es crear un nuevo elemento en el array "global" al cual le añadimos el polinomio total módulo p, sus factores, un vector sums que en la posición i tiene el número k sii hay k maneras de sumar i combinando los grados de los polinomios de los factores (para calcular esto se hace por programación dinámica sobre el array aux[2]) y un vector de sets en el que la posición i tiene un set. Para cada posible combinación \sum_{j=0}^{k_0} deg(f_j) = i el set guarda el degTag de f_{k_0}.
 *
 * Guardarse este set sirve para reconstruir de manera eficiente las posibles sumas, funciona de la siguiente manera. Para generar la primera suma, supongamos que el vector intersection nos dice que la suma vale i. En ese caso cogeremos el primer elemento del set de i y nos fijamos en qué valor tiene su grado (sea d éste), si el grado es exactamente i entonces ya tenemos una configuración. Si no, nos fijamos en el set de i-d y cogemos el primer elemento y comprobamos si el grado es exactamente i-d y si no lo es procedemos recursivamente. La implementación de esto se hace guardando en una pila cada uno de los elementos de los sets que vamos recorriendo y guardamos también los índices que determinan de qué set eran los elementos. Llevamos también una pila con polinomios en la que la base de la pila tiene el polinomio indicado por el tag de la base de la pila de elementos del set, y en general en una posición cualquiera n de la pila llevamos el producto de la posición anterior multiplicado por el polinomio que determina el tag de la posición n de la pila de elementos de los sets. Esto se hace para generar eficientemente los polinomios que devolvemos en cada iteración (cada iteración será amortizada cte + una division de polinomios). ¿Qué se hace en cada iteración? Se coge el elemento de la pila de los elementos del set y se quita la cima y se sustituye por el siguiente elemento del último set. Si éste se nos ha acabado, hacemos otro pop y cogemos el siguiente elemento del set que nos ha quedado y seguimos hasta reconstruir otra solución.
 * */

//TODO: Por otro lado los sets los tiro a la basura salvo el del bueno, por lo que quizá sea mejor calcularlo una vez que ya sepamos cuál es el bueno


void HenselSubsets::insert(const std::vector<std::pair<Fpxelem_b, unsigned int> > & factors, const Fpxelem_b & poli){
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
            global[ind].map[tag]= i;
            for (unsigned int k = 1 ; k <= semiSumOfDeg; k++){//Sólo nos hace falta calcular las cosas hasta la mitad del grado
                if (aux[which][k] != 0 && k + deg <= semiSumOfDeg){
                    aux[1-which][k+deg] += aux[which][k];
                    global[ind].predecessor[k+deg].insert({deg, tag});
                }
            }
            if (deg <= semiSumOfDeg){
                aux[1-which][deg]++;
                global[ind].predecessor[deg].insert({deg, tag});
            }
            aux[which] = aux[1-which];
            which = 1 - which;
            tag++;
        }
    } //Ahora tengo en aux[which] un vector con las multiplicidades de las posibles sumas
    if(verbose){
        std::cout << "Posibles sumas para el primo " << factors[0].first.getSize() << std::endl;
        for (unsigned int i = 1; i <= semiSumOfDeg; i++){
                if (aux[which][i] != 0) {
                    std::cout << i << " ";
                }
            }
        std::cout << std::endl;

        for (auto aaa: global[ind].factors){
            if(aaa.second == 1) std::cout << aaa.first << std::endl;
            else std::cout << "("<< aaa.first << ")^" << aaa.second <<std::endl;
        }
    }

    //<end Quitar>

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
                }
            }
        }
    }
    global[ind].sums = aux[which];

    if(verbose){
        std::cout << "Numero de casos con este primo teniendo los resultados anteriores: " << global[ind].numOfCases << std::endl;
        std::cout << "=====" << std::endl;
    }

    /*
     for (auto j : global[ind].predecessor){
        std::cout << ":";
        for (auto k: j){
            std::cout << k.deg << " " << k.tag<< "; ";
        }
        std::cout << ":" << std::endl;
    }
    */
}

Option HenselSubsets::bestOption(){
    if (index == -1){
        if (global.size() == 0) return {false, Fpxelem_b(), Fpxelem_b()};//The polynomials are irrelevant
        unsigned int min = global[0].numOfCases;
        index = 0;
        for (unsigned int i = 1; i < global.size(); i++) {
            if (min > global[i].numOfCases){
                min = global[i].numOfCases;
                index = i;
            }
        }
        globind = std::move(global[index]);
        global.clear();//TODO: Gestionar bien este destructor. Destruyo todo, es memoria que ya no necesito. En teoría así debería valer.
        if(verbose){
            std::cout << "El primo que he escogido es " << globind.pol.getSize() << " y hay " << min << " posibilidades" << std::endl;
            std::cout << "Esta es la factorizacion modulo el primo:";
            for (auto aaa: globind.factors)
                std::cout << "("<< aaa.first << ")^" << aaa.second <<std::endl;
            std::cout << "=====" << std::endl;
        }
    }
    if (stackIt.empty()){
        hadRemoved = false;
        while (++index_intersection <= semiSumOfDeg &&
                (intersection[index_intersection] == 0 || globind.predecessor[index_intersection].empty()) );
        if (index_intersection >= semiSumOfDeg+1) return {false, Fpxelem_b(), Fpxelem_b() };//The polynomials are irrelevant
        stackIt.push(globind.predecessor[index_intersection].begin());
        stackInd.push(index_intersection);
        stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
        while (stackIt.top()->deg != stackInd.top()){//Esto es facil que falle si no se programa bien, si da un error al depurar busca aquí
            if (globind.predecessor[stackInd.top() - stackIt.top()->deg].empty() ||
                globind.predecessor[stackInd.top() - stackIt.top()->deg].begin()->tag >= stackIt.top()->tag){

                auto it = stackIt.top();
                stackIt.pop();
                it++;
                while(!stackIt.empty() &&
                        (it == globind.predecessor[stackInd.top()].end() ||
                         it->tag >= stackIt.top()->tag)){
                    it = stackIt.top();
                    it++;
                    stackIt.pop(); stackPol.pop(); stackInd.pop();
                }
                if (it == globind.predecessor[stackInd.top()].end()){
                    stackPol.pop(); stackInd.pop();
                    return bestOption();
                }
                stackPol.pop();
                stackInd.pop();
                if (stackInd.empty())
                    stackInd.push(index_intersection);
                else
                    stackInd.push(stackInd.top() - stackIt.top()->deg);
                stackIt.push(it);
                if (stackPol.empty())
                    stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
                else
                    stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);
            }
            else{
                stackInd.push(stackInd.top() - stackIt.top()->deg);
                stackIt.push(globind.predecessor[stackInd.top()].begin());
                stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);

            }
        }
        return {true, stackPol.top(), globind.pol / stackPol.top() };
    }
    else{
        if (index_intersection > semiSumOfDeg) return {false, Fpxelem_b(), Fpxelem_b() };//The polynomials are irrelevant
        auto it = stackIt.top();
        stackIt.pop();
        if (hadRemoved)
            hadRemoved = false; //At this point stackIt.size() should always be 1.
        else
            it++;
        while(!stackIt.empty() && (it == globind.predecessor[stackInd.top()].end() || it->tag >= stackIt.top()->tag)){
            it = stackIt.top();
            it++;
            stackIt.pop(); stackPol.pop(); stackInd.pop();
        }
        if (it == globind.predecessor[stackInd.top()].end()){
            stackPol.pop(); stackInd.pop();
            return bestOption();
        }
        stackPol.pop();
        stackInd.pop();
        if (stackInd.empty())
            stackInd.push(index_intersection);
        else
            stackInd.push(stackInd.top() - stackIt.top()->deg);
        stackIt.push(it);
        if (stackPol.empty())
            stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
        else
            stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);

        while (stackIt.top()->deg != stackInd.top()){//Esto es facil que falle si no se programa bien, si da un error al depurar busca aquí
            if (globind.predecessor[stackInd.top() - stackIt.top()->deg].empty() ||
                globind.predecessor[stackInd.top() - stackIt.top()->deg].begin()->tag >= stackIt.top()->tag){

                auto it = stackIt.top();
                stackIt.pop();
                it++;
                while(!stackIt.empty() && (it == globind.predecessor[stackInd.top()].end() || it->tag >= stackIt.top()->tag)){
                    it = stackIt.top();
                    it++;
                    stackIt.pop(); stackPol.pop(); stackInd.pop();
                }
                if (it == globind.predecessor[stackInd.top()].end()){
                    stackPol.pop(); stackInd.pop();
                    return bestOption();
                }
                stackPol.pop();
                stackInd.pop();
                if (stackInd.empty())
                    stackInd.push(index_intersection);
                else
                    stackInd.push(stackInd.top() - stackIt.top()->deg);
                stackIt.push(it);
                if (stackPol.empty())
                    stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
                else
                    stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);
            }
            else{
                stackInd.push(stackInd.top() - stackIt.top()->deg);
                stackIt.push(globind.predecessor[stackInd.top()].begin());
                stackPol.push(stackPol.top() * globind.factors[ globind.map[stackIt.top()->tag] ].first);
            }
        }
        return {true, stackPol.top(), globind.pol / stackPol.top() };
    }
}


void HenselSubsets::removeFirstLastOption(Zxelem_b w){
    globind.pol /= stackPol.top();
    last = w;
    sumOfDeg -= stackPol.top().deg();
    semiSumOfDeg = sumOfDeg/2;
    while(stackIt.size() > 1){
        DegTag dt = *stackIt.top();
        for (unsigned int i = 1; i <= semiSumOfDeg; i++)
            globind.predecessor[i].erase(dt);
        stackIt.pop(); stackPol.pop(); stackInd.pop();
    }
    DegTag dt = *stackIt.top();//El último elemento lo borro al final para poder hacer bien lo de avanzar el puntero

    if (++stackIt.top() == globind.predecessor[stackInd.top()].end()){
        stackIt.pop(); stackPol.pop(); stackInd.pop();
    }
    else{
        stackPol.pop(); stackInd.pop();
        stackPol.push(globind.factors[ globind.map[stackIt.top()->tag] ].first);
        stackInd.push(index_intersection);
    }

    for (unsigned int i = 1; i <= semiSumOfDeg; i++)
        globind.predecessor[i].erase(dt);
    hadRemoved = true; 
}

Zxelem_b HenselSubsets::getLast(){
    return last;
}
