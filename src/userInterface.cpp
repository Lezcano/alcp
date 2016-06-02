#include "userInterface.hpp"

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include <cstdarg>
#include <algorithm>
#include <utility>
#include <cstring>

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
#include "bchCodes.hpp"

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
        cmds["BCH"] = std::unique_ptr<Command>(new CommandBCH());
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
        if (alcpScan(args, "(s)$", &s)) {
            if (UserInterface::instance().isCommand(s))
				UserInterface::instance().callHelp(s);
			else
                std::cout << "Parse error" << std::endl;
        }
        else if (alcpScan(args, "$"))
            UserInterface::instance().help();
		else
       	    std::cout << "Parse error" << std::endl;
    }

    void UserInterface::CommandHelp::help(const std::string &name) {
        std::cout << "Outputs a list with help for all the commands of the system or help for a single command, if specified" << std::endl;
        std::cout << "        " << name << std::endl;
        std::cout << "        " << name << "(command)" << std::endl;
    }

    void UserInterface::CommandBerlekamp::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
            std::vector< std::vector<big_int> > f;
            big_int p;

			if (alcpScan(args, "(f,v,n)$", &f, &v, &p)){
                Fqxelem_b pol(f, Fq_b(Fpxelem_b(v, p)));
                std::cout << "The factors of the polynomial:" << std::endl << "    " << pol << std::endl;
                std::cout << "are the following:" << std::endl;
                auto factors = factorizationBerlekamp(pol);
                for (auto &pair: factors) {
                	if (pair.second != 1)
                        std::cout << "( ";
                    std::cout << pair.first;
                    if (pair.second != 1)
                        std::cout << " )^" << pair.second;
                    std::cout << std::endl;
                }
            }
			else if (alcpScan(args, "(v,n)$", &v, &p)){
                Fpxelem_b pol(v, p);
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
				std::cout << "Parse error" << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
    }

    void UserInterface::CommandBerlekamp::help(const std::string &name) {
        std::cout <<
        "Factorizes a polynomial in GF(p^m)[x] using Berlekamp algorithm." << std::endl <<
        "    GF(p^m) refers to the finite field of p^m elements, with p prime and m a postive integer." << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), p, m)" << std::endl;
    }

    void UserInterface::CommandCantorZassenhaus::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
            std::vector< std::vector<big_int> > f;
            big_int p;

			if (alcpScan(args, "(f,v,n)$", &f, &v, &p)){
                Fqxelem_b pol(f, Fq_b(Fpxelem_b(v, p)));
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
			else if (alcpScan(args, "(v,n)$", &v, &p)){
                Fpxelem_b pol(v, p);
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
				std::cout << "Parse error" << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
    }

    void UserInterface::CommandCantorZassenhaus::help(const std::string &name) {
        std::cout <<
        "Factorizes a polynomial in GF(p^m)[x] using Cantor-Zassenhaus algorithm." << std::endl <<
        "    GF(p^m) refers to the finite field of p^m elements, with p prime and m a postive integer." << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
        std::cout << "        " << name << "(((a0_0, a_01, ...a0_n0), ..., (am_0, ..., ak_nk)), p, (b_0, ..., b_m))" << std::endl;
        std::cout << "    " << "The first k+1 vectors will be polynomials the coefficients" << std::endl;
    }

    void UserInterface::CommandHensel::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
			if (!alcpScan(args, "(v)$", &v))
				std::cout << "Parse error" << std::endl;
            Zxelem_b pol(v);
            auto factors = factorizationHensel(pol);
            std::cout << "The factors of the polynomial:" << std::endl <<
            			 "    " << pol << std::endl;
            for (auto &pair: factors) {
                if (pair.second != 1)
                    std::cout << "( ";
                std::cout << pair.first;
                if (pair.second != 1)
                    std::cout << " )^" << pair.second;
                std::cout << std::endl;
            }
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
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
                std::cout << "Parse error" << std::endl; 
            Zxelem_b pol1(v1), pol2(v2);
            std::cout << "gcd( " << pol1 << " ,  " << pol2 << " ) = " << std::endl;
            std::cout << modularGCD(pol1, pol2) << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
    }

    void UserInterface::CommandModularGCD::help(const std::string &name) {
        std::cout << "Computes the greatest common divisor of two polynomials in Z[x] using the modular gcd algorithm." << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ..., a_n), (b_0, b_1, ..., b_m))" << std::endl;
    }

    void UserInterface::CommandCRA::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> m, u;
			if (!alcpScan(args, "(v,v)$", &m, &u) || m.size() != u.size())
                std::cout << "Parse error" << std::endl;
            std::cout << integerCRA(m, u) << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
    }

    void UserInterface::CommandCRA::help(const std::string &name) {
        std::cout <<
        "Given positive moduli m_i in Z which are relatively prime and " <<
        " given corresponding residues u_i in Z_{m_i} it computes the unique integer u in Z_m " <<
        " (where m = \\prod m_i) such that u = u_i (mod m_i)" <<  std::endl <<
        "    The behavior is not specified if m_i are not relatively prime." <<  std::endl;
        std::cout << "        " << name << "((m_0, m_1, ..., m_n), (u_0, u_1, ..., u_n)" << std::endl;
    }

    void UserInterface::CommandEEA_ED::parseAndRun(std::istringstream &args) {
        big_int a, b, p;
        std::vector<big_int> v1, v2, i;
        std::vector<std::vector<big_int>> vv1, vv2;

        try {
            if( alcpScan(args, "(n,n)$", &a, &b)){
                big_int x,y;
                auto g = eea(a,b,x,y);
                std::cout << "gcd(" << a << "," << b << ")";
                std::cout << " = "  << g << " = ";
                std::cout << a << "*";
                if(x < 0) std::cout << "(" << x << ")";
                else std::cout << x;
                std::cout << "+";

                if(b < 0) std::cout << "(" << b << ")";
                else std::cout << b;
                std::cout << "*";

                if(y < 0) std::cout << "(" << y << ")";
                else std::cout << y;
                std::cout << std::endl;
            }
            else if(alcpScan(args, "(v,v,n)$", &v1, &v2, &p)){
                Fpxelem_b a1(v1, p), b1(v2, p), x, y;
                auto g = eea(a1,b1,x,y);
                std::cout << "gcd(" << a1 << ", " << b1 << ")" << std::endl <<
                        " = " << g << " (mod  " << p << ") = ";

                if(a1.deg() != 0)   std::cout << "(" << a1 << ")";
                else std::cout << a1;
                std::cout << "*";

                if(x.deg() != 0)   std::cout << "(" << x << ")";
                else std::cout << x;
                std::cout << "+";

                if(b1.deg() != 0)   std::cout << "(" << b1 << ")";
                else std::cout << b1;
                std::cout << "*";

                if(y.deg() != 0)   std::cout << "(" << y << ")";
                else std::cout << y;
                std::cout << std::endl;
            }
            else if(alcpScan(args, "(f,f,n,v)$", &vv1, &vv2, &p, &i)){
                Fpxelem_b mod(i,p);
                Fq_b f(mod);
                Fqxelem_b a1(vv1,f), b1(vv2,f), x, y;

                auto g = eea(a1,b1,x,y);
                std::cout << a1 << "*" << x << "+" << b1 << "*" << y <<
                " = " << g << " (mod  " << mod << ") = gcd(" << a1 << ", " << b1 << ")" << std::endl;
            }
            else std::cout << "Parse error" << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
    }

    void UserInterface::CommandEEA_ED::help(const std::string &name) {
        std::cout << "Given two elements of an Euclidean Domain, it returns the coefficientes of the Bezout identity and their greatest common divisor." << std::endl;
        std::cout << "    Current supported Euclidean Domains are: Z, GF(p)[X], GF(p^m)[X]." << std::endl;
        std::cout << "        " << name << "(a,b)" << std:: endl;
        std::cout << "        " << name << "((a_0, ...,a_n), (b_0, ..., b_n), p)" << std:: endl;
        std::cout << "        " << name << "(((a0_0, a_01, ...a0_n0), ..., (am_0, ..., ak_nk)), ((b0_0, b_01, ...b0_n0), ..., (bm_0, ..., bk_nk)), p, (i_0, ..., i_m))" << std::endl;
    }

    void UserInterface::CommandPollardFactor::parseAndRun(std::istringstream &args) {
        try {
            big_int aux;
            if (!alcpScan(args, "(n)$", &aux))
                std::cout << "Parse error" << std::endl; 
            std::cout << "WARNING: This algorithm is probabilistic. Some integers might not fully factorize." << std::endl;
            std::cout << aux << " = ";
            auto map = factorInteger(static_cast<long long>(aux));

            auto it = map.begin();

            std::cout << it->first;
            if(it->second != 1)
                std::cout << "^" << it->second;

            it++;
            for(; it != map.end(); ++it){
                std::cout << " * ";
                std::cout << it->first;
                if(it->second != 1)
                    std::cout << "^" << it->second;
            }
            std::cout << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
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
                std::cout << "Parse error" << std::endl; 
            pollardRhoLogarithm(2, 5, 1019, log);
            std::cout << "log_" << a << "(" << b << ") = " << log << " (mod " << p << ")" << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
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
                std::cout << "Parse error" << std::endl; 
            std::cout << num << " is ";
            if (!millerRabin(num))
                std::cout << "not ";
            std::cout << "prime" << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
    }

    void UserInterface::CommandMillerRabin::help(const std::string &name) {
        std::cout <<
        "Given an integer, the algorithm determines whether it is prime. " << std::endl <<
         "    The output is correct with a very high probability." << std::endl;
        std::cout << "        " << name << "(a)" << std::endl;
    }

    void UserInterface::CommandIrrGFp::parseAndRun(std::istringstream &args) {
        try {
            std::vector<big_int> v;
            big_int p;
            if (!alcpScan(args, "(v,n)$", &v, &p))
                std::cout << "Parse error" << std::endl; 
            Fpxelem_b pol(v,p);
            std::cout << "The polynomial " << pol << " is ";
            if (!pol.irreducible())
                std::cout << "not ";
            std::cout << "irreducible in " << pol.getField() << std::endl;
        } catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
    }

    void UserInterface::CommandIrrGFp::help(const std::string &name) {
        std::cout <<
        "Given a polynomial in GF(p)[x], the algorithm determines whether it is irreducible." << std::endl <<
		"    GF(p) refers to the finite field of p elements, with p prime." << std::endl;
        std::cout << "        " << name << "((a_0, a_1, ... a_n), p)" << std::endl;
    }

    void UserInterface::CommandBCH::parseAndRun(std::istringstream &args){
    	try {
    		size_t l, c, d;
			std::vector<big_int> v;
			big_int p;
			if (alcpScan(args, "(n, v, m, m, m)$", &p, &v, &l, &c, &d)){
				Fpxelem_b prim_pol(v, p);
				BCH bch = BCH(prim_pol, 1, l, c, d);
				std::cout << "Give me the message(Intro for random message): " << std::endl;
				std::string cmdline;
				std::cout << ">> ";
				getline(std::cin, cmdline);
				Fpxelem_b message;
				if (cmdline == ""){
					message = randomPol(Fp_b(p), bch.getDimension()-1);
					std::cout << "The random message chosen is" << std::endl;
					std::cout << "MESSAGE := " << message << std::endl;
				}
				else{
					std::istringstream ss(cmdline);
					if (alcpScan(args, "v", &v)){
						message = Fpxelem_b(v, p);
					}
					else{
						std::cout << "Parse error" << std::endl; //TODO throw exception
					}
				}
				std::cout << std::endl <<"The codification of your message is:" << std::endl;
				Fpxelem_b sent = bch.encode(message);
				std::cout << "CODE SENT:= " << sent << std::endl;
				std::cout << "Sending your message through an error-prone channel" << std::endl;
				std::set<int> error_indices;
				Fpxelem_b received;
				std::tie(error_indices, received) = bch.randomErrors(sent);
				std::cout << "An oracle tells us that " << error_indices.size() <<
						" ocurred during the transmission at indices ";
				for(auto & elem: error_indices){
					std::cout << elem << ", ";
				}
				std::cout << "The code received is" << std::endl;
				std::cout << "CODE RECEIVED:= " << received << std::endl;
				std::cout << std::endl << std::endl << "Decoding using Berlekamp algorithm."<< std::endl;
				Fpxelem_b decoded = bch.decode(received);
				std::cout << "The corrected code is" << std::endl;
				std::cout << "CODE CORRECTED := " << decoded << std::endl << std::endl;
				if (sent == decoded){
					std::cout << "I have checked it is the same code that was sent" << std:: endl;
					std::cout << "The original message is" << decoded/bch.getG()<< std:: endl;
				}
				else{
					std::cout << "OH! YOUR HAVE A MISTAKE. YOUR CODE HAS NOT WORKED.";
				}
			}
			else
				std::cout << "Parse error" << std::endl;
    	} catch (const ExcepALCP& e){
			std::cout << e.msg() << std::endl;
			throw e;
        }
	}

	void UserInterface::CommandBCH::help(const std::string &name) {
		//TODO
		std::cout << "Create a BCH codification system and then ask for a message to encode. The program changes some random coefficients of the encoded message and then it decodes it";
		std::cout << std::endl << std::endl;
		std::cout << "        " << "Usage:" << std::endl;
		std::cout << "        " << name << "(p, (a_0, a_1, ..., a_n), l, c, d)" << std::endl;
		std::cout << "        " << "Atributes description" << std::endl;
		std::cout << "        " << "p: prime number that defines the message base field" << std::endl;
		std::cout << "        " << "(a_0, a_1, ..., a_n): are the coefficients of a primitive polinomial f in F_p[x] such that there is a primitive l^th root of unity in F_p[t]/<f>, where l is also given." << std::endl;
		std::cout << "        " << "l: see previous description."<< std::endl;
		std::cout << "        " << "c, d: The generating polynomial of the BHC code will be the lcm of the minimum polinomials of t^c, t^{c+1}, ... t^{c+d-2} \\in F_p[t]/<f>" << std::endl;
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

	bool inline isBigNumber(std::istringstream &args, big_int &n) {
		if (!(args >> n))
			return false;
		return true;
	}
	bool inline isNumber(std::istringstream &args, int &n) {
			if (!(args >> n))
				return false;
			return true;
		}

	bool isVector(std::istringstream &args, std::vector<big_int> &v) {
		std::string s;
		if (args.get() != '(')
			return false;
		big_int num;
		if (!alcpScan(args, "n", &num))
			return false;
		v.push_back(num);
		while (args.peek() != ')') {
			if (!alcpScan(args, ",n", &num))
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
		if (!alcpScan(args, "v", &v))
			return false;
		vv.push_back(v);
		while (args.peek() != ')') {
			v.clear();
			if (!alcpScan(args, ",v", &v))
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
        // If The iss is empty, we return false
        if(iss.rdbuf()->in_avail() == 0){
            return (!std::strcmp(fmt, "$") || !std::strcmp(fmt, ""));
        }

        // Get the position of the string to restore it later if there has been any parsing error
        int pos = iss.tellg();

		while (*fmt != '\0' && *fmt != '$') {
            while(iss.peek() == ' ' || iss.peek() == '\t')
                iss.ignore();


			if (*fmt == 'n'){
                big_int n;
				if (!isBigNumber(iss, n)) {
                    iss.clear();
                    iss.seekg(pos);
                    return false;
                }
				auto *num = va_arg(args, big_int*);
				*num = n;
			}
			else if (*fmt == 'm'){
                int m;
				if (!isNumber(iss, m)) {
                    iss.clear();
                    iss.seekg(pos);
                    return false;
                }
				auto *num = va_arg(args, int*);
				*num = m;
			}
			else if (*fmt == 'v') {
                std::vector<big_int> v;
				if (!isVector(iss, v)) {
                    iss.clear();
                    iss.seekg(pos);
                    return false;
                }
				auto v2 = va_arg(args, std::vector<big_int>*);
				*v2 = v;
			}
			else if (*fmt == 'f') {
                std::vector<std::vector<big_int>> vv;
				if (!isVectorOfVectors(iss, vv)) {
                    iss.clear();
                    iss.seekg(pos);
                    return false;
                }
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
				if (c != *fmt) {
                    iss.clear();
                    iss.seekg(pos);
                    return false;
                }
			}

            while(iss.peek() == ' ' || iss.peek() == '\t')
                iss.ignore();

			++fmt;
		}

		if(*fmt == '$' && *(fmt+1) == '\0'){
            // Check if the iss is empty
			if(iss.rdbuf()->in_avail() != 0) {
                iss.clear();
                iss.seekg(pos);
                return false;
            }
			return true;
		}
		else if(*fmt == '\0')
			return true;

		iss.clear();
		iss.seekg(pos);
		return false;

	}

}
