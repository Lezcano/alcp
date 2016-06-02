#ifndef __BCH_CODES
#define __BCH_CODES

#include <vector>
#include <utility>
#include <set>

#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "fqelem.hpp"
#include "fqxelem.hpp"

#include "types.hpp"


namespace alcp {
	class BCH {
	public:
		//TODO: para hacer mensajes random, usar el generador de polinomios aleatorios del archivo de factorización.
		BCH(const Fpxelem_b &primitive_poly, size_t n, size_t l, size_t c, size_t d);

		Fpxelem_b encode(Fpxelem_b message);

		std::pair<std::set<int>, Fpxelem_b > randomErrors(Fpxelem_b v);

		Fpxelem_b decode(Fpxelem_b message);

		Fpxelem_b getG() const;

		size_t getDimension() const;
	private:
		Fqelem_b alpha;
		Fpxelem_b g; //To generalize, this could be an Fqxelem
		Fq_b field_ext;
		size_t dimension;
		size_t distance;
		size_t length;
		size_t maxErrors;
		size_t c;

	};


}



#endif /* BCHC_ODES_HPP */

