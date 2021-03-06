#ifndef __EXCEPTIONS_HPP
#define __EXCEPTIONS_HPP

#include <string>
#include <iosfwd>

namespace alcp {
    class ExcepALCP {
    public:
        ExcepALCP() { }

        ExcepALCP(const std::string &msg) : _msg(msg) { }

        const std::string msg() const { return _msg; }

        friend std::ostream &operator<<(std::ostream &out, const ExcepALCP &e);

    protected:
        std::string _msg;
    };

    inline std::ostream &operator<<(std::ostream &out, const ExcepALCP &e) {
        out << e._msg;
        return out;
    }


// Macro defined in order to not repeat
// many times the same structure
#define DECLARE_EXCEPTION(Exception) \
class Exception : public ExcepALCP { \
public: \
Exception() {} \
Exception(const std::string &msg) : ExcepALCP(msg) {} \
}

    /**
     *  Thrown when a Field with a non prime p is declared
     */
    DECLARE_EXCEPTION(EpNotPrime);

    /**
     *  Thrown when two numbers which are not in the same field are added
     */
    DECLARE_EXCEPTION(EOperationUnsupported);

    /**
     *  Yet to be implemented functionality
     */
    DECLARE_EXCEPTION(ENotImplementedYet);

    /**
     *  Thrown when tying to create an element of Fp[X] with a vector with elemnts
     *   not in F_p (but in other Fpelem)
     *  Also thrown when trying to assing to an element in some ED another element
     *   which is not in a compatible ED (e.g. an element assign an element from F_3
     *   an element in F_5;
     */
    DECLARE_EXCEPTION(ENotCompatible);

    /**
     * Thrown when trying to create an element of F[X] with an empty vector
     */
    DECLARE_EXCEPTION(EEmptyVector);
    DECLARE_EXCEPTION(EFPXNotIrreducible);
    DECLARE_EXCEPTION(EDifferentSizeVectorsCRA);
    DECLARE_EXCEPTION(ENotInitializedRing);

    /**
	 * Thrown when trying to create a BCH object and the data provided does not fulfill the requirements
	 */
    DECLARE_EXCEPTION(EBadBCHInitialization);
    /**
	 * Thrown when trying to encode a message longer that it is possible
	 */
    DECLARE_EXCEPTION(EBadFormatMessageBCH);

    /**
   	 * Thrown when trying to decode a message with more errors than those that are allowed
   	 */
    DECLARE_EXCEPTION(ETooManyErrorsBCH);

 
 
}

#endif // __EXCEPTIONS_HPP
