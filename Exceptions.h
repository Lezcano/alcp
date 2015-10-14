#ifndef __EXCEPTIONS_H
#define __EXCEPTIONS_H

#include <string>
#include <iosfwd>

/**
 * Every exception inherits from this class
 *
 * This is not needed in C++17
 */

class ExcepALCP {
public:
	ExcepALCP() {}
	ExcepALCP(const std::string &msg) : _msg(msg) {}

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
Exception() {}; \
Exception(const std::string &msg) : ExcepALCP(msg) {} \
};

/**
 *  Thrown when a Field with a non prime p is declared
 */
DECLARE_EXCEPTION(EpNotPrime);

/**
 *  Thrown when two numbers which are not in the same field are added
 */
DECLARE_EXCEPTION(EOperationUnsupported);

#endif // __EXCEPTIONS_H
