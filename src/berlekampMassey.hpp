#ifndef __BERLEKAMP_MASSEY
#define __BERLEKAMP_MASSEY

#include <vector>
#include <utility>

#include "fqelem.hpp"
#include "fqxelem.hpp"

namespace alcp {
	/*TODO: Se puede en vez de usar polinomios, usar vectores de Felems
	 * que vale igual y la parte de multiplicar por x^m sería más
	 * eficiente entre otras cosas
	 *
	 * f_m: polinomio que había la última vez que el grado cambió.
	 * d_m: su discrepancia
	 * f_k: polinomio actual
	 * d_k: su discrepancia
	 * l:	grado de f_k en cada momento
	 * m:	número de etapas que han pasado desde que se actualizó el grado
	 *
	 *
	 */
	template<typename Fxelem>
	Fxelem berlekampMassey(const std::vector<typename Fxelem::Felem> & s){
		/*/
		for(auto elem: s){
			std::cout << elem << std::endl;
		}
		/**/
		unsigned int n = s.size();
		Fxelem f_m(getOne(s[0])), f_k(getOne(s[0]));
		unsigned int l = 0;
		int m = 1;
		typename Fxelem::Felem d_m = getOne(s[0]);

		for (unsigned int i = 0; i < n; i++){
			typename Fxelem::Felem d_k = getZero(s[0]);
			unsigned int min = f_k.deg();
			if (l < min ) min = l;
			for(size_t j = 0; j <= min; j++ ){
				d_k += s[i-j]*f_k[j];
			}
			//std::cout << "d_k" << d_k << std::endl;
			if (d_k == 0){
				m++;
			}
			else if (2*l <= i ){
				Fxelem	aux(f_k);
				std::vector<typename Fxelem::Felem> x_m(m+1, getZero(s[0]));
				x_m[m] = 1;
				f_k -= (d_k/d_m)*f_m*Fxelem(x_m);
				//std::cout << f_k << std::endl;
				f_m = aux;
				m = 1;
				l = i + 1 - l;
				d_m = d_k;
			}
			else{
				std::vector<typename Fxelem::Felem> x_m(m+1, getZero(s[0]));
				x_m[m] = 1;
				f_k -= (d_k/d_m)*f_m*Fxelem(x_m);
				//std::cout << f_k << std::endl;
				m++;
			}
		}
		return f_k;
	}
}


#endif /* __BERLEKAMP_MASSEY */
