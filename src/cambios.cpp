
        // Base field
        using F = Fp;
        using Felem = Fpelem;

Faltan:
    // Externas
    to_string(Felem) // Friend

    bool getZero(const Fpelem &e);
    unit()
    normalForm()
    compatible(lhs, rhs) 
            if(_f.getSize() != rhs._f.getSize())
                throw ENotCompatible("Asignation failed. The vectors "+ to_string()+ " and " + rhs.to_string() + " are not in the same field.");

    - In PolinomialRing<Felem>(const std::vector...)
            ll _p=_f.getSize();
            // Check integrity of v
            for (auto &i : _v)
                if(i.getSize()!=_p)
                    throw ENotCompatible("Not all the elements in the array are in the same field!");

        bool operator==(ll rhs)const{
            return *this == PolinomialRing<Felem>(_f.get(rhs));
        }

        bool operator!=(ll rhs)const{
            return *this != PolinomialRing<Felem>(_f.get(rhs));
        }

        // Prime p of the base field Fp[X]
        ll getSize()const{return _f.getSize();}

        const F getField()const{return _f;}

