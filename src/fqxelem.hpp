// Implementation of a GF(p) field
#ifndef __FQXELEM_HPP
#define __FQXELEM_HPP

#include <vector>

#include "fqelem.hpp"
#include "polRing.hpp"

namespace alcp {
    template<class Integer>
    class Fqxelem : public PolynomialRing<Fqxelem, Fqelem<Integer>, Integer> {
    private:
        using FBase = PolynomialRing<Fqxelem, Fqelem<Integer>, Integer>;
    public:
        // Base field
        using F = Fq<Integer>;
        using Felem = Fqelem<Integer>;

        // Inherit ctor
        using FBase::PolynomialRing;

        Fqxelem() = default;

        Fqxelem (const std::vector<std::vector<Integer>> & v, const Fq<Integer>& f) : FBase(
                [](const std::vector<std::vector<Integer>>& v1, const Fq<Integer>& f1) -> Fqxelem {
                    std::vector<Fqelem<Integer>> ret(v1.size());
                    std::transform(v1.begin(),
                                   v1.end(),
                                   ret.begin(),
                                   [&f1](const std::vector<Integer>& v2){
                                       return f1.get(Fpxelem<Integer>(v2,f1.getP()));
                                   });
                    return ret;
                }(v,f)){}

        Fqxelem (std::vector<std::vector<Integer>> && v, Fq<Integer>&& f) : FBase(
                [](std::vector<std::vector<Integer>>&& v1, Fq<Integer>&& f1) -> Fqxelem {
                    std::vector<Fqelem<Integer>> ret(v1.size());
                    std::transform(std::make_move_iterator(v1.begin()),
                                   std::make_move_iterator(v1.end()),
                                   ret.begin(),
                                   [&f1](std::vector<Integer>&& v2){
                                       return f1.get(Fpxelem<Integer>(std::move(v2),f1.getP()));
                                   });
                    return ret;
                }(std::move(v),std::move(f))){}

        Fqxelem (const std::vector<Fpxelem<Integer> > &e, const Fq<Integer> & f) : FBase(
		   [&]() -> Fqxelem {
			   std::vector<Fqelem<Integer>> ret(e.size());
			   std::transform(e.begin(),
							  e.end(),
							  ret.begin(),
							  [&f](const Fpxelem<Integer>& v){
								  return f.get(v);
							  });
			   return ret;
		   }()){}

        const Fq<Integer> getField() const {
            return this->lc().getField();
        }

        Integer getSize() const {
            return this->getField().getSize();
        }

        friend Fqxelem<Integer> getZero(const Fqxelem<Integer> &e) { return Fqxelem<Integer>(getZero(e.lc())); }

        friend Fqxelem<Integer> getOne(const Fqxelem<Integer> &e) { return Fqxelem<Integer>(getOne(e.lc())); }

        friend Fqxelem<Integer> unit(const Fqxelem<Integer> &e) { return e.lc(); }

        friend bool compatible(const Fqxelem<Integer> &lhs, const Fqxelem<Integer> &rhs) {
            return lhs.getField() == rhs.getField();
        }

        friend bool operator==(const Fqxelem<Integer> &lhs, Integer rhs) {
            return lhs.deg() == 0 && lhs.lc() == lhs.getField().get(rhs);
        }

        friend bool operator==(Integer lhs, const Fqxelem<Integer> &rhs) { return rhs == lhs; }

        friend bool operator!=(const Fqxelem<Integer> &lhs, Integer rhs) { return !(lhs == rhs); }

        friend bool operator!=(Integer lhs, const Fqxelem<Integer> &rhs) { return !(rhs == lhs); }

        /*
        // Este diseño es puta mierda, mirar como hacer bien el toFpxeleme, toZxelemSym y esta función (en polRing?)
        friend Fqxelem<Integer> toFqxelem(const Fpxelem<Integer> &e, const Fq<Integer> & f) {
            std::vector<Fqelem<Integer>> v(e.deg() + 1);
            //TODO we should do this with transform
            for (std::size_t i = 0; i <= e.deg(); ++i)
                v[i] = f.get(e[i]);
            return Fqxelem<Integer>(std::move(v));
        }
        */
    };

    using Fqxelem_b = Fqxelem<big_int>;

    //template class Fqxelem<big_int>;
}

#endif
