// Inteface of a finite field
#ifndef __F_H
#define __F_H

#include"types.h"

class F : public ED{
    public:
        virtual FElem eea(const FElem &e1,const FElem &e2);


    private:
        class Fp{
            public:
                Fp(myInt p){
                    _p = p;
                }
            private:
        }
}

#endif // __FELEM_H
