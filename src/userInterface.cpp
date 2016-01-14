#include "userInterface.hpp"

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>

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

namespace alcp{
	class Command;
///////////////////////////////////////////////////////////////////////////
// User interface
///////////////////////////////////////////////////////////////////////////

	UserInterface::UserInterface(){
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
	UserInterface & UserInterface::instance(){
		static UserInterface s;
		return s;
	}
	bool UserInterface::isCommand (const std::string & s){
		return cmds.find(s) != cmds.end();
	}

	void UserInterface::run(){
		std::string cmdline;
		std::cout << "Computer algebra system by Mario Lezcano and David Martinez" << std::endl;
		std::cout << "Write \"help\" for help" << std::endl;
		while(true){
			std::cout << ">> ";
			getline(std::cin, cmdline);
			std::stringstream ss(cmdline);
			getline(ss, cmdline, '(');

			if(cmdline == "quit")
				break;

			auto it = cmds.find(cmdline);

			if (it == cmds.end()){
				std::cout << "Unrecognized option" << std::endl;
				std::cout << "Write \"help\" for help" << std::endl;
				//this->help();
			}
			else{
				(it->second)->parseAndRun(ss);
			}
		}
	}

	void UserInterface::help(){
		for (auto &pair: cmds){
			std::cout << pair.first << std::endl;
			std::cout << "    ";
			(pair.second)->help(pair.first);
		}
	}

	void UserInterface::callHelp(const std::string & s){
		cmds[s]->help(s);
	}

/////////////////////////////////////////////////////////////////////////////
// Commands
/////////////////////////////////////////////////////////////////////////////
	bool comma(std::stringstream & args);
	bool isVector(std::stringstream & args, std::vector<big_int> & v);
	bool closedParen(std::stringstream & args);
	bool end(std::stringstream & args);
	
	void UserInterface::CommandHelp::parseAndRun(std::stringstream & args){
		std::string s;
		char c;
		if ((args >> c) && c == '(' ){
			getline(args, s, ')');
			if (!end(args) || !UserInterface::instance().isCommand(s))
				throw 1;
			UserInterface::instance().callHelp(s);
		}
		else{
			if (!end(args))
				throw 1;
			UserInterface::instance().help();
		}
	}

	void UserInterface::CommandHelp::help( const std::string & name){
		std::cout << "Outputs a list with help for all the commands of the system or help for a single command, if specified" << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name;
		std::cout<< "   " << name << "(command)" << std::endl;
	}

	void UserInterface::CommandBerlekamp::parseAndRun(std::stringstream & args){
		try{
			std::vector<big_int> v;
			big_int p;
			std::size_t exp;
			char c;
			if (!isVector(args, v) || !comma(args) || !(args >> p) || !(args >> c))
				throw 1; //throw new ParseError();
			if ( c == ','){
				if (!(args >> exp) ||  !closedParen(args))
					throw 1; //throw new ParseError();
				Fq_b f(p, exp);
				std::vector<Fqelem_b> vf;
				for (unsigned int i = 0; i < v.size(); i++){
					vf.push_back(f.get(v[i]));
				}
				Fqxelem_b pol(vf);
				auto factors = factorizationBerlekamp(pol);
				std::cout << "Factors:" << std::endl;
				for (auto &pair: factors){
					if (pair.second != 1)
						std::cout << "( ";
					std::cout << pair.first;
					if (pair.second != 1)
						std::cout << ")^" << pair.second;
					std::cout << std::endl;
				}
			}
			else{
				if (c != ')' || !end(args)){
					throw 1; //throw new ParseError();
				}
				Fp_b f(p);

				std::vector<Fpelem_b> vf;
				for (unsigned int i = 0; i < v.size(); i++){
					vf.push_back(f.get(v[i]));
				}
				Fpxelem_b pol(vf);
				auto factors = factorizationBerlekamp(pol);
				std::cout << "The factors of the polynomial:" << std::endl << "    " << pol << std::endl;
				std::cout << "are the following:" << std::endl;
				for (auto &pair: factors){
					if (pair.second != 1)
						std::cout << "( ";
					std::cout << pair.first;
					if (pair.second != 1)
						std::cout << " )^" << pair.second;
					std::cout << std::endl;
				}
			}

		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandBerlekamp::help( const std::string & name){
		std::cout << "Factorizes a polynomial in GF(p^m)[x] using Berlekamp algorithm, being GF(p^m) the finite field with p^m elements, with p prime and m a natural number." << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
		std::cout<< "   " << name << "((a_0, a_1, ..., a_n), p, m)" << std::endl;
	}
	void UserInterface::CommandCantorZassenhaus::parseAndRun(std::stringstream & args){
		try{
			std::vector<big_int> v;
			big_int p;
			std::size_t exp;
			char c;
			if (!isVector(args, v) || !comma(args) || !(args >> p) || !(args >> c))
				throw 1; //throw new ParseError();
			if ( c == ','){
				if (!(args >> exp) ||  !closedParen(args))
					throw 1; //throw new ParseError();
				Fq_b f(p, exp);
				std::vector<Fqelem_b> vf;
				for (unsigned int i = 0; i < v.size(); i++){
					vf.push_back(f.get(v[i]));
				}
				Fqxelem_b pol(vf);
				auto factors = factorizationCantorZassenhaus(pol);
				std::cout << "The factors of the polynomial:" << std::endl << "    " << pol << std::endl;
				std::cout << "are the following:" << std::endl;
				for (auto &pair: factors){
					if (pair.second != 1)
						std::cout << "( ";
					std::cout << pair.first;
					if (pair.second != 1)
						std::cout << ")^" << pair.second;
					std::cout << std::endl;
				}
			}
			else{
				if (c != ')' || !end(args))
					throw 1; //throw new ParseError();
				Fp_b f(p);
				std::vector<Fpelem_b> vf;
				for (unsigned int i = 0; i < v.size(); i++){
					vf.push_back(f.get(v[i]));
				}
				Fpxelem_b pol(vf);
				auto factors =  factorizationCantorZassenhaus(pol);
				std::cout << "Factors:" << std::endl;
				for (auto &pair: factors){
					if (pair.second != 1)
						std::cout << "( ";
					std::cout << pair.first;
					if (pair.second != 1)
						std::cout << ")^" << pair.second;
					std::cout << std::endl;
				}
			}
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandCantorZassenhaus::help( const std::string & name){
		std::cout << "Factorizes a polynomial in GF(p^m)[x] using Cantor-Zassenhaus algorithm, being GF(p^m) the finite field with p^m elements, with p prime and m a natural number." << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
		//TODO
		std::cout << "   " << name << "(((a0_0, a_01, ...a0_n0), ..., (am_0, ..., ak_nk)), p, (b_0, ..., b_m))"<< std::endl;
		std::cout << "      " << "The first k+1 vectors will be polynomials the coefficients";
	}
	void UserInterface::CommandHensel::parseAndRun(std::stringstream & args){
		try{
			std::vector<big_int> v;
			if ( !isVector(args, v) || !closedParen(args) )
				throw 1; //throw new ParseError();
			Zxelem_b pol(v);
			auto factors = factorizationHensel(pol);
			std::cout << "Factors:" << std::endl;
			for (auto &pair: factors){
				if (pair.second != 1)
					std::cout << "( ";
				std::cout << pair.first;
				if (pair.second != 1)
					std::cout << " )^" << pair.second;
				std::cout << std::endl;
			}
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandHensel::help( const std::string & name){
		std::cout << "Factorizes a polynomial in Z[x] using Hensel algorithm" << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "((a_0, a_1, ..., a_n))" << std::endl;
	}
	void UserInterface::CommandModularGCD::parseAndRun(std::stringstream & args){
		try{
			std::vector<big_int> v1, v2;
			if (!isVector(args, v1) || !comma(args) ||
					 !isVector(args, v2)	|| closedParen(args))
				throw 1; //throw new ParseError();
			Zxelem_b pol1(v1), pol2(v2);
			std::cout << modularGCD(pol1, pol2) << std::endl;
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandModularGCD::help( const std::string & name){
		std::cout << "Computes the greatest common divisor of two polynomials in Z[x] using the modular gcd algorithm." << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "((a_0, a_1, ..., a_n), (b_0, b_1, ..., b_m))" << std::endl;
	}
	void UserInterface::CommandCRA::parseAndRun(std::stringstream & args){
		try{
			big_int aux;
			std::vector<big_int> m, u;
			if ( !(args >> aux) || !comma(args))
				throw 1;
			m.push_back(aux);
			if ( !(args >> aux) )
				throw 1;
			u.push_back(aux);
			while (comma(args) && args >> aux ){
				m.push_back(aux);
				if (!comma(args) || !(args >> aux))
					throw 1; //throw new ParseError();
				u.push_back(aux);
			}
			if (!closedParen(args))
				throw 1;
			std::cout << integerCRA(m, u) << std::endl;
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandCRA::help( const std::string & name){
		std::cout << "Given positive moduli m_i in Z (0 <= i <= n) which are relatively prime and given corresponding residues u_i in Z_{m_i} t computes the unique integer u in Z_m (where m = \\prod m_i) such that u = u_i (mod m_i) i = 0,...,n. The behavior is not specified if m_i are not relatively prime." << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "((m_0, u_0, m_1, u_1, ..., m_n, u_n)" << std::endl;
	}

	void UserInterface::CommandEEA_ED::parseAndRun(std::stringstream & args){
		try{
//TODO
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandEEA_ED::help( const std::string & name){
	}
	void UserInterface::CommandPollardFactor::parseAndRun(std::stringstream & args){
		try{
			big_int aux;
			if (!(args >> aux) || !closedParen(args))
				throw 1; //throw new ParseError();
			//TODO	//std::cout << factorizationPollardRhoBrent(aux);
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandPollardFactor::help( const std::string & name){
		std::cout << "Given an integer, it returns its factors" << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "(a)" << std::endl;
	}
	void UserInterface::CommandPollardLog::parseAndRun(std::stringstream & args){
		try{
			//TODO
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandPollardLog::help( const std::string & name){
	}
	void UserInterface::CommandMillerRabin::parseAndRun(std::stringstream & args){
		try{
			big_int num;
			if (!(args >> num) || !closedParen(args))
				throw 1; //throw new ParseError();
			std::cout << num << " is ";
			if(!millerRabin(num))
				std::cout << "not ";
			std::cout << "prime" << std::endl;
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandMillerRabin::help( const std::string & name){
		std::cout << "Given an integer, the algorithm determines whether a number is prime. The output is correct with a very high probability." << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "(a)" << std::endl;
	}
	void UserInterface::CommandIrrGFp::parseAndRun(std::stringstream & args){
		try{
			std::vector<big_int> v;
			big_int p;
			if (!isVector(args, v) || !comma(args) || !(args >> p) || !closedParen(args))
				throw 1; //throw new ParseError();
			Fp_b f(p);
			std::vector<Fpelem_b> vf;
			for (unsigned int i = 0; i < vf.size(); i++){
				vf.push_back(f.get(v[i]));
			}
			std::cout << "The polynomial is ";
			Fpxelem_b pol(vf);
			if (! pol.irreducible())
				std::cout << "not ";
			std::cout << "irreducible" << std::endl;
		}catch(...){
			std::cout << "Parse error" << std::endl;
		}
	}
	void UserInterface::CommandIrrGFp::help( const std::string & name){
		std::cout << "Given a polynomial in GF(p)[x], the algorithm determines whether  it is irreducible. GF(p) refers to the finite field with p elements where p is prime" << std::endl;
		std::cout << "FORMAT" << std::endl;
		std::cout<< "   " << name << "((a_0, a_1, ... a_n), p)" << std::endl;
	}

	bool comma(std::stringstream & args){
		char c;
		if (! (args >> c) || c != ',')
			return false;
		return true;
	}

	bool isVector(std::stringstream & args, std::vector<big_int> & v){
		std::string s;
		char c;
		if (!(args >> c) || c != '(')
			return false;
		getline(args, s, ')');
		std::stringstream aux(s);
		big_int num;
		if (!(aux >> num))
			return false;
		v.push_back(num);
		while (aux >> c >> num){
			if (c != ',')
				return false;
			v.push_back(num);
		}
		return true;
	}
	bool closedParen(std::stringstream & args){
		char c;
		if (!(args >> c) || c!= ')' )
			return false;
		return true;
	}

	bool end(std::stringstream & args){
		char c;
		if (args >> c)
			return false;
		return true;
	}
}

