#include "userInterface.hpp"

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include <cstdarg>


#include "types.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "fqelem.hpp"
#include "zxelem.hpp"
#include "exceptions.hpp"
#include "factorizationFq.hpp"
#include "integerCRA.hpp"
#include "hensel.hpp"
#include "modularGCD.hpp"
#include "generalPurpose.hpp"

namespace alcp {
    class Command;
    ///////////////////////////////////////////////////////////////////////////
    // User interface
    ///////////////////////////////////////////////////////////////////////////

    UserInterface::UserInterface() {
        cmds["help"] = std::unique_ptr<Command>(new CommandHelp());
        cmds["factorBerlekamp"] = std::unique_ptr<Command>(new CommandBerlekamp());
        cmds["factorCZ"] = std::unique_ptr<Command>(new CommandCantorZassenhaus());
        cmds["factorHensel"] = std::unique_ptr<Command>(new CommandHensel());
        cmds["modularGCD"] = std::unique_ptr<Command>(new CommandModularGCD());
        cmds["chineseRA"] = std::unique_ptr<Command>(new CommandCRA());
        cmds["eea_ed"] = std::unique_ptr<Command>(new CommandEEA_ED());
        cmds["factorInteger"] = std::unique_ptr<Command>(new CommandPollardFactor());
        cmds["discreteLog"] = std::unique_ptr<Command>(new CommandPollardLog());
        cmds["isPrime"] = std::unique_ptr<Command>(new CommandMillerRabin());
        cmds["isIrreducibleGFp"] = std::unique_ptr<Command>(new CommandIrrGFp());
    }

    UserInterface &UserInterface::instance() {
        static UserInterface s;
        return s;
    }

    bool UserInterface::isCommand(const std::string &s) {
        return cmds.find(s) != cmds.end();
    }

	void unrecognizedOp(){
		std::cout << "Unrecognized option" << std::endl;
		std::cout << "Write \"help\" for help" << std::endl;
	}
	bool alcpScan(std::istringstream &iss, const char *fmt, ...);

    void UserInterface::run() {
        std::string cmdline;
        std::cout << "Computer algebra system by Mario Lezcano and David Martinez" << std::endl;
        std::cout << "Write \"help\" for help" << std::endl;
        while (true) {
            std::cout << ">> ";
            getline(std::cin, cmdline);
            if (cmdline == "quit")
                break;

            std::istringstream ss(cmdline);
			std::string cmd;
			if (!alcpScan(ss, "s", &cmd)){
				unrecognizedOp();
				continue;
			}

            auto it = cmds.find(cmd);
            if (it == cmds.end()){
				unrecognizedOp();
				continue;
            }
            else
                (it->second)->parseAndRun(ss);
        }
    }

    void UserInterface::help() {
        for (auto &pair: cmds) {
            std::cout << pair.first << std::endl;
            std::cout << "    ";
            (pair.second)->help(pair.first);
        }
    }

    void UserInterface::callHelp(const std::string &s) {
        cmds[s]->help(s);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Commands
    /////////////////////////////////////////////////////////////////////////////

    bool isVector(std::istringstream &args, std::vector<big_int> &v);

    bool closedParen(std::istringstream &args);

    bool end(std::istringstream &args);

    void UserInterface::CommandHelp::parseAndRun(std::istringstream &args) {
        std::string s;
        args.sync();
  		std::istringstream args2(args.str().substr(args.tellg()));
        if (alcpScan(args, "(s)$", &s)) {
            if (!UserInterface::instance().isCommand(s))
                throw 1;
            UserInterface::instance().callHelp(s);
        }
        else {
            if (!alcpScan(args2, "$"))
                throw 1;
            UserInterface::instance().help();
        }
    }

    void UserInterface::CommandHelp::help(const std::string &name) {
        std::cout <<
        "Outputs a list with help for all the commands of the system or help for a single command, if specified" <<
        std::endl;
        std::cout << "        " << name << std::endl;
        std::cout << "        " << name << "(command)" << std::endl;
    }

    void UserInterface::CommandBerlekamp::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
            std::vector< std::vector<big_int> > f;
            big_int p;
			std::istringstream args2(args.str().substr(args.tellg()));
			if (alcpScan(args, "(f,v,n)$", &f, &v, &p)){
                Fq_b f(p, v.size() - 1);//v.size()-1 is the degree of the polynomial, thus it is the exponent of the size of the field
                std::vector<Fqelem_b> vf;
                for (unsigned int i = 0; i < v.size(); i++) {
                    vf.push_back(f.get(v[i]));
                }
                Fqxelem_b pol(vf);
                auto factors = factorizationBerlekamp(pol);
                std::cout << "The factors of the polynomial:" << std::endl << "    " << pol << std::endl;
                std::cout << "are the following:" << std::endl;
                for (auto &pair: factors) {
                    if (pair.second != 1)
                        std::cout << "( ";
                    std::cout << pair.first;
                    if (pair.second != 1)
                        std::cout << " )^" << pair.second;
                    std::cout << std::endl;
                }
            }
			else
				if (alcpScan(args2, "(v,n)$", &v, &p)){
                Fp_b f(p);

                std::vector<Fpelem_b> vf;
                for (unsigned int i = 0; i < v.size(); i++) {
                    vf.push_back(f.get(v[i]));
                }
                Fpxelem_b pol(vf);
                auto factors = factorizationBerlekamp(pol);
                std::cout << "The factors of the polynomial:" << std::endl << "    " << pol << std::endl;
                std::cout << "are the following:" << std::endl;
                for (auto &pair: factors) {
                    if (pair.second != 1)
                        std::cout << "( ";
                    std::cout << pair.first;
                    if (pair.second != 1)
                        std::cout << " )^" << pair.second;
                    std::cout << std::endl;
                }
            }
			else
				throw 1;

        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandBerlekamp::help(const std::string &name) {
        std::cout <<
        "Factorizes a polynomial in GF(p^m)[x] using Berlekamp algorithm, being GF(p^m) the finite field with p^m elements, with p prime and m a natural number." <<
        std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), p, m)" << std::endl;
    }

    void UserInterface::CommandCantorZassenhaus::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
            std::vector< std::vector<big_int> > f;
            big_int p;
			std::istringstream args2(args.str().substr(args.tellg()));

			if (alcpScan(args, "(f,n,v)$", &f, &p, &v)){
                Fq_b f(p, v.size() - 1);//v.size()-1 is the degree of the polynomial, thus it is the exponent of the size of the field
                std::vector<Fqelem_b> vf;
                for (unsigned int i = 0; i < v.size(); i++) {
                    vf.push_back(f.get(v[i]));
                }
                Fqxelem_b pol(vf);
                auto factors = factorizationCantorZassenhaus(pol);
                std::cout << "The factors of the polynomial:" << std::endl << "    " << pol << std::endl;
                std::cout << "are the following:" << std::endl;
                for (auto &pair: factors) {
                    if (pair.second != 1)
                        std::cout << "( ";
                    std::cout << pair.first;
                    if (pair.second != 1)
                        std::cout << ")^" << pair.second;
                    std::cout << std::endl;
                }
            }
			else if (alcpScan(args2, "(v,n)$", &v, &p)){
                Fp_b f(p);
                std::vector<Fpelem_b> vf;
                for (unsigned int i = 0; i < v.size(); i++) {
                    vf.push_back(f.get(v[i]));
                }
                Fpxelem_b pol(vf);
                auto factors = factorizationCantorZassenhaus(pol);
                std::cout << "The factors of the polynomial:" << std::endl << "    " << pol << std::endl;
                std::cout << "are the following:" << std::endl;
                for (auto &pair: factors) {
                    if (pair.second != 1)
                        std::cout << "( ";
                    std::cout << pair.first;
                    if (pair.second != 1)
                        std::cout << ")^" << pair.second;
                    std::cout << std::endl;
                }
            }
			else
				throw 1;
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandCantorZassenhaus::help(const std::string &name) {
        std::cout <<
        "Factorizes a polynomial in GF(p^m)[x] using Cantor-Zassenhaus algorithm, being GF(p^m) the finite field with p^m elements, with p prime and m a natural number." <<
        std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
        //TODO
        std::cout << "        " << name << "(((a0_0, a_01, ...a0_n0), ..., (am_0, ..., ak_nk)), p, (b_0, ..., b_m))" <<
        std::endl;
        std::cout << "    " << "The first k+1 vectors will be polynomials the coefficients";
    }

    void UserInterface::CommandHensel::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
			if (!alcpScan(args, "(v)$", &v))
				throw 1;
            Zxelem_b pol(v);
            auto factors = factorizationHensel(pol);
            std::cout << "Factors:" << std::endl;
            for (auto &pair: factors) {
                if (pair.second != 1)
                    std::cout << "( ";
                std::cout << pair.first;
                if (pair.second != 1)
                    std::cout << " )^" << pair.second;
                std::cout << std::endl;
            }
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandHensel::help(const std::string &name) {
        std::cout << "Factorizes a polynomial in Z[x] using Hensel algorithm" << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n))" << std::endl;
    }

    void UserInterface::CommandModularGCD::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v1, v2;
			if (!alcpScan(args, "(v,v)$", &v1, &v2))
                throw 1; //throw new ParseError();
            Zxelem_b pol1(v1), pol2(v2);
            std::cout << "gcd( " << pol1 << " ,  " << pol2 << " ) = " << std::endl;
            std::cout << modularGCD(pol1, pol2) << std::endl;
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandModularGCD::help(const std::string &name) {
        std::cout <<
        "Computes the greatest common divisor of two polynomials in Z[x] using the modular gcd algorithm." << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), (b_0, b_1, ..., b_m))" << std::endl;
    }

    void UserInterface::CommandCRA::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> m, u;
			if (!alcpScan(args, "(v,v)$", &m, &u) || m.size() != u.size())
                throw 1;
            std::cout << integerCRA(m, u) << std::endl;
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandCRA::help(const std::string &name) {
        std::cout <<
        "Given positive moduli m_i in Z (0 <= i <= n) which are relatively prime and given corresponding residues u_i in Z_{m_i} t computes the unique integer u in Z_m (where m = \\prod m_i) such that u = u_i (mod m_i) i = 0,...,n. The behavior is not specified if m_i are not relatively prime." <<
        std::endl;
        std::cout << "        " << name << "((m_0, m_1, ..., m_n), (u_0, u_1, ..., u_n)" << std::endl;
    }

    void UserInterface::CommandEEA_ED::parseAndRun(std::istringstream &args) {
        try {
            big_int a, b, p;
            std::vector<big_int> v1, v2;
            std::vector<big_int> i;
            std::vector<std::vector<big_int>> vv1, vv2;
            std::string s(args.str().substr(args.tellg()));

            // MALLLLL no deberia ser std::istringstream && args??
            // TODO a medias
            //if(alcpScan(args2, "(n,n)$", &a, &b)){ }
            //if(alcpScan(std::istringstream(s), "(n,n)$", &a, &b)){ }
            //else if(alcpScan(args2, "(n,n,n)$", &a, &b, &p)){}
            //else if(alcpScan(args2, "(v,v,n)$", &v1, &v2, &p)){}
            //else if(alcpScan(args2, "(v,v,v)$", &v1, &v2, &i)){}
            //else if(alcpScan(args2, "(f,f,n,v)$", &vv1, &vv2, &p, &i)){}
            //else throw 1;
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandEEA_ED::help(const std::string &name) {
        std::cout << "Given two elements of an Euclidean Domain, it returns the coefficientes of the Bezout identity and their greatest common divisor." << std::endl;
        std::cout << "Currently supported Euclidean Domains are: Z, GF(p), GF(p)[X]), GF(p^m), GF(p^m)[X]." << std::endl;
        std::cout << "        " << name << "(a,b)" << std:: endl;
        std::cout << "        " << name << "(a,b,p)" << std:: endl;
        std::cout << "        " << name << "((a_0, ...,a_n), (b_0, ..., b_n), p)" << std:: endl;
        std::cout << "        " << name << "((a_0, ...,a_n), (b_0, ..., b_n), (i_0, ..., i_{n+1}))" << std:: endl;
        std::cout << "        " << name << "(((a0_0, a_01, ...a0_n0), ..., (am_0, ..., ak_nk)), ((b0_0, b_01, ...b0_n0), ..., (bm_0, ..., bk_nk)), p, (i_0, ..., i_m))" << std::endl;
    }

    void UserInterface::CommandPollardFactor::parseAndRun(std::istringstream &args) {
        try {
            big_int aux;
            if (!alcpScan(args, "(n)$", &aux))
                throw 1; //throw new ParseError();
            std::cout << "WARNING: This algorithm is probabilistic. Some integers might not fully factorize." << std::endl;
            std::cout << "The number " << aux << " factorizes as:" << std::endl;
            //std::cout << factorizationPollardRhoBrent(aux);
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandPollardFactor::help(const std::string &name) {
        std::cout << "Given an integer, it returns its factors" << std::endl;
        std::cout << "        " << name << "(a)" << std::endl;
    }

    void UserInterface::CommandPollardLog::parseAndRun(std::istringstream &args) {
        try {
            long long a, b, p, log;
            if (!alcpScan(args, "(n,n,n)$", &a, &b, &p))
                throw 1; //throw new ParseError();
            pollardRhoLogarithm(2, 5, 1019, log);
            std::cout << "log_" << a << "(" << b << ") = " << log << " (mod " << p << ")" << std::endl;
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandPollardLog::help(const std::string &name) {
        std::cout <<
        "Given a prime number p and two positive integers a, b, with a being a generator of the cyclic multiplicative group of order p-1, it computes log_a(b) in F_p" << std::endl;
        std::cout << "        " << name << "(a)" << std::endl;
    }

    void UserInterface::CommandMillerRabin::parseAndRun(std::istringstream &args) {
        try {
            big_int num;
            if (!alcpScan(args, "(n)$", &num))
                throw 1; //throw new ParseError();
            std::cout << num << " is ";
            if (!millerRabin(num))
                std::cout << "not ";
            std::cout << "prime" << std::endl;
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandMillerRabin::help(const std::string &name) {
        std::cout <<
        "Given an integer, the algorithm determines whether a number is prime. The output is correct with a very high probability." <<
        std::endl;
        std::cout << "        " << name << "(a)" << std::endl;
    }

    void UserInterface::CommandIrrGFp::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
            big_int p;
            if (alcpScan(args, "(v,n)$", &v, &p))
                throw 1; //throw new ParseError();
            Fp_b f(p);
            std::vector<Fpelem_b> vf;
            for (unsigned int i = 0; i < vf.size(); i++) {
                vf.push_back(f.get(v[i]));
            }
            std::cout << "The polynomial is ";
            Fpxelem_b pol(vf);
            if (!pol.irreducible())
                std::cout << "not ";
            std::cout << "irreducible" << std::endl;
        } catch (...) {
            std::cout << "Parse error" << std::endl;
        }
    }

    void UserInterface::CommandIrrGFp::help(const std::string &name) {
        std::cout <<
        "Given a polynomial in GF(p)[x], the algorithm determines whether  it is irreducible. GF(p) refers to the finite field with p elements where p is prime" <<
        std::endl;
        std::cout << "        " << name << "((a_0, a_1, ... a_n), p)" << std::endl;
    }

    bool closedParen(std::istringstream & args){
        char c;
        if (!(args >> c) || c!= ')' )
            return false;
        return true;
    }

    bool end(std::istringstream & args){
        char c;
        if (args >> c)
            return false;
        return true;
    }

	bool inline isNumber(std::istringstream &args, big_int &n) {
		if (!(args >> n))
			return false;
		return true;
	}

	bool isVector(std::istringstream &args, std::vector<big_int> &v) {
		std::string s;
		if (args.get() != '(')
			return false;
		big_int num;
		if (!isNumber(args, num))
			return false;
		v.push_back(num);
		while (args.peek() != ')') {
			if (args.get() != ',' || !isNumber(args, num))
				return false;
			v.push_back(num);
		}
		args.ignore();
		return true;
	}

	bool isVectorOfVectors(std::istringstream &args, std::vector<std::vector<big_int>> &vv) {
		std::string s;
		std::vector<big_int> v;
		if (args.get() != '(')
			return false;
		if (!isVector(args, v))
			return false;
		vv.push_back(v);
		while (args.peek() != ')') {
			v.clear();
			if (args.get() != ',' || !isVector(args, v))
				return false;
			vv.push_back(v);
		}
		args.ignore();
		return true;
	}

	void isString(std::istringstream &iss, std::string & s){
		// We  extract the remaining string in the iss
		s = iss.str().substr(iss.tellg());
		std::size_t lst = s.find_first_of("(), \t");
		if (lst == std::string::npos) {
			iss.str("");
			return;
		}
		s = s.substr(0, lst);
		iss.ignore(lst);
	}
	bool alcpScan(std::istringstream &iss, const char *fmt, ...) {
		va_list args;
		va_start(args, fmt);
        std::size_t fst;

		while (*fmt != '\0' && *fmt != '$') {
            fst = iss.str().substr(iss.tellg()).find_first_not_of(" \t");
            iss.ignore(fst);
			if (*fmt == 'n') {
                big_int n;
				if (!isNumber(iss, n))
					return false;
				auto *num = va_arg(args, big_int*);
				*num = n;
			}
			else if (*fmt == 'v') {
                std::vector<big_int> v;
				if (!isVector(iss, v))
					return false;
				auto v2 = va_arg(args, std::vector<big_int>*);
				*v2 = v;
			}
			else if (*fmt == 'f') {
                std::vector<std::vector<big_int>> vv;
				if (!isVectorOfVectors(iss, vv))
					return false;
				auto vv2 = va_arg(args, std::vector<std::vector<big_int>>*);
				*vv2 = vv;
			}
			else if(*fmt == 's') {
                std::string s;
				isString(iss, s);
				auto s2 = va_arg(args, std::string*);
				*s2 = s;
			}
			else {
				char c;
				iss.get(c);
				if (c != ' ' && c != *fmt)
					return false;
			}
            fst = iss.str().substr(iss.tellg()).find_first_not_of(" \t");
            iss.ignore(fst);
			++fmt;
		}
		if(*fmt == '$' && *(fmt+1) == '\0'){
			std::string aux = iss.str().substr(iss.tellg());
			if(aux.find_first_not_of(" ") != std::string::npos)
				return false;
			return true;
		}
		else if(*fmt == '\0')
			return true;
		else {
			std::cerr << "Debug: Bad format in alcpScan" << std::endl;
			return false;
		}
	}

}
