#include "userInterface.hpp"
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

#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace alcp {


	bool coma(std::stringstream & args);
	bool isVector(std::stringstream & args, std::vector<big_int> & v);
	bool parenCl(std::stringstream & args);

	class CommandHelp: public Command {
		void parseAndRun(std::stringstream & args){
			string s;
			char c;
			if ((args >> c) && c == '(' ){
				getline(args, s, ')');
				if (!args.eof() || !UserInterface.isCommand(s))
					throw 1;
				UserInterface.getCommand(s).help();	
			}
			else{
				if (!args.eof())
					throw 1;
				UserInterface.help();

			}
		}
		void help(string & name){
			std::cout << "Outputs a list with help for all the commands of the system or help for a single command, if specified" << std::endl;
			std::cout << "FORMAT" << std::endl;
			std:cout<< "   " << name;
			std:cout<< "   " << name << "(command)" << std::endl;
		}
	}

	class CommandBerlekamp: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				big_int p;
				big_int exp;
				char c;
				if (!isVector(args, v) || !coma(args) || !(args >> p) || !parenCl(args))
					throw 1; //throw new ParseError();
				if ( (args.str())[0] != ')'){
					if (!(args >> exp) ||  !parenCl(args))
						throw 1; //throw new ParseError();
					Fq_b f(p, exp);
					std::vector<Fqelem_b> vf;
					for (unsigned int i = 0; i < vf.size(); i++){
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
					Fp_b f(p);
					std::vector<Fpelem_b> vf;
					for (unsigned int i = 0; i < vf.size(); i++){
						vf.push_back(f.get(v[i]));
					}
					Fpxelem_b pol(vf);
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

			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
			std::cout << "Factorizes a polynomial in GF(p^m)[x] using Berlekamp algorithm, being GF(p^m) the finite field with p^m elements, with p prime and m a natural number." << std::endl;
			std::cout << "FORMAT" << std::endl;
			std:cout<< "   " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
			std:cout<< "   " << name << "((a_0, a_1, ..., a_n), p, m)" << std::endl;
		}
	};
	class CommandCantorZassenhaus: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				big_int p;
				big_int exp;
				char c;
				if (!isVector(args, v) || !coma(args) || !(args >> p) || !parenCl(args))
					throw 1; //throw new ParseError();
				if ( (args.str())[0] != ')'){
					if (!(args >> exp) || !parenCl(args))
						throw 1; //throw new ParseError();
					Fq_b f(p, exp);
					std::vector<Fqelem_b> vf;
					for (unsigned int i = 0; i < vf.size(); i++){
						vf.push_back(f.get(v[i]));
					}
					Fqxelem_b pol(vf);
					auto factors = factorizationCantorZassenhaus(pol);
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
					Fp_b f(p);
					std::vector<Fpelem_b> vf;
					for (unsigned int i = 0; i < vf.size(); i++){
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
		void help(string & name){
			std::cout << "Factorizes a polynomial in GF(p^m)[x] using Cantor-Zassenhaus algorithm, being GF(p^m) the finite field with p^m elements, with p prime and m a natural number." << std::endl;
			std::cout << "FORMAT" << std::endl;
			std:cout<< "   " << name << "((a_0, a_1, ..., a_n), p)" << std::endl;
			std:cout<< "   " << name << "((a_0, a_1, ..., a_n), p, m)"<< std::endl;
		}
	};
	class CommandHensel: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				char c;
				if ( !isVector(args, v) || !parenCl(args) )
					throw 1; //throw new ParseError();
				Zxelem_b pol(v);
				auto factors = factorizationHensel(pol);
				std::cout << "Factors:" << std::endl;
				for (auto &pair: factors){
					if (pair.second != 1)
						std::cout << "( ";
					std::cout << pair.first;
					if (pair.second != 1)
						std::cout << ")^" << pair.second;
					std::cout << std::endl;
				}
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
			std::cout << "Factorizes a polynomial in Z[x] using Hensel algorithm" << std::endl;
			std::cout << "FORMAT" << std::endl;
std:cout<< "   " << name << "((a_0, a_1, ..., a_n))" << std::endl;
		}
	};
	class CommandModularGCD: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v1, v2;
				char c;
				if (!isVector(args, v1) || !coma() ||
						|| !isVector(args, v2)	|| parenCl(args))
					throw 1; //throw new ParseError();
				Zxelem_b pol1(v1), pol2(v2);
				std::cout << modularGCD(pol1, pol2) << std::endl;
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
			std::cout << "Computes the greatest common divisor of two polynomials in Z[x] using the modular gcd algorithm." << std::endl;
			std::cout << "FORMAT" << std::endl;
std:cout<< "   " << name << "((a_0, a_1, ..., a_n), (b_0, b_1, ..., b_m))" << std::endl;
		}
	};
	class CommandCRA: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				big_int aux;
				std::vector<big_int> m, u;
				if ( !(args >> aux) || !coma(args))
					throw 1;
				m.push_back(aux);
				if ( !(args >> aux) )
					throw 1;
				u.push_back(aux);
				while (coma(args) && args >> aux ){
					m.push_back(aux);
					if (!coma(args) || !(args >> aux))
						throw 1; //throw new ParseError();
					u.push_back(aux);
				}
				if (!parenCl(args))
					throw 1;
				std::cout << integerCRA(m, u) << std::endl;
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
			std::cout << "Given positive moduli m_i \in Z (0 \leq i \leq n) which are relatively prime and given corresponding residues u_i \in Z_{m_i} t computes the unique integer u \in Z_m (where m = \prod m_i) such that u = u_i (mod m_i) i = 0,...,n. The behavior is not specified if m_i are not relatively prime." << std::endl;
			std::cout << "FORMAT" << std::endl;
			std:cout<< "   " << name << "((m_0, u_0, m_1, u_1, ..., m_n, u_n)" << std::endl;
		}
	};
	class CommandEEA_ED: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
	//TODO
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
		}
	};
	class CommandPolardFactor: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				big_int aux;
				char c;
				if (!(args >> aux) || !parenCl(args))
					throw 1; //throw new ParseError();
				//TODO	//std::cout << factorizationPollardRhoBrent(aux);
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
			std::cout << "Given an integer, it returns its factors" << std::endl;
			std::cout << "FORMAT" << std::endl;
			std:cout<< "   " << name << "(a)" << std::endl;
		}
	};
	class CommandPolardLog: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				//TODO
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
		}
	};
	class CommandMillerRabin: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				big_int num;
				char c;
				if (!(args >> num) || !parenCl(args))
					throw 1; //throw new ParseError();
				std::cout << num << "is ";
				if(!millerRabin(num))
					std::cout << "not ";
				std::cout << "prime" << std::endl;
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
		void help(string & name){
			std::cout << "Given an integer, the algorithm determines whether a number is prime. The output is correct with a very high probability.
			std::cout << "FORMAT" << std::endl;
			std:cout<< "   " << name << "(a)" << std::endl;
		}
	};
	class CommandIrrGFp: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				big_int p;
				if (!isVector(args, v) || !coma(args) || !(args >> p) || !parenCl(args))
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
		void help(string & name){
			std::cout << "Given a polynomial in GF(p)[x], the algorithm determines whether  it is irreducible. GF(p) refers to the finite field with p elements where p is prime" << std::endl;
			std::cout << "FORMAT" << std::endl;
			std:cout<< "   " << name << "((a_0, a_1, ... a_n), p)" << std::endl;
		}
	};

	UserInterface::UserInterface(){
		cmds = {
			std::make_pair<std::string, Command>(
					"help", CommandHelp());
			std::make_pair<std::string, Command>(
					"factorBerlekamp", CommandBerlekamp());
			std::make_pair<std::string, Command>(
					"factorCZ", CommandCantorZassenhaus());
			std::make_pair<std::string, Command>(
					"factorHensel", CommandHensel());
			std::make_pair<std::string, Command>(
					"modularGCD", CommandModularGCD());
			std::make_pair<std::string, Command>(
					"chineseRA", CommandCRA());
			std::make_pair<std::string, Command>(
					"eea_ed", CommandEEA_ED());
			std::make_pair<std::string, Command>(
					"factorInteger", CommandPolardFactor());
			std::make_pair<std::string, Command>(
					"discreteLog", CommandPolardLog());
			std::make_pair<std::string, Command>(
					"isPrime", CommandMillerRabin());
			std::make_pair<std::string, Command>(
					"isIrreducibleGFp", CommandIrrGFp());
		};
		//cmds[""] = new CommandMillerRabin();
	}
	static void UserInterface::run(){
		std::string cmdline;
		std::cout << "Computer algebra system by Mario Lezcano and David Martinez" << std::endl;
		std::cout << "Write \"help\" for help" << std::endl;
		std::cout << ">> ";

		std::stringstream ss;
		getline(ss, cmdline);
		while (ss.str() != "quit"){
			std::string cmd;
			getline(ss, cmd, '(');
			auto it = cmds.find(cmd);
			if ( it == cmds.end()){
				std::cout << "Unrecognized option";
				this->help();
			}
			else{
				(it->second).parseAndRun(ss);
			}
		}
	}
	static void UserInterface::help(){
		for (auto &pair: cmds){
			std::cout << pair.first << std::endl;
			std::cout << "    " << (pair.second).help(pair.first);
		}
	}
	static bool isCommand (const string & s){
		return cmds.find(s) != cmds.end();
	}
	static Command && getCommand(const string & s){
		return cmds[s];
	}

	bool coma(std::stringstream & args){
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
			if ( args.eof() )
				return false;
			std::stringstream aux(s);
			big_int num;
			if (!args >> num)
				return false;
			v.push_back(num);
			while (aux >> c >> num){
				if (c != ',')
					return false;
				v.push_back(num);
			}
			return true;
		}
	bool parenCl(std::stringstream & args){
		char c;
		if (!parenCl(args))
			return false;
		return true;
	}
}
