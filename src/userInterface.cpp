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

	bool isVector(std::stringstream & args, std::vector<big_int> & v){
			std::string s;
			getline(args, s, ')');
			if ( args.eof() )
				return false;
			std::stringstream aux(s);
			big_int num;
			while (aux >> num)
				v.push_back(num);
			return true;
		}
	class CommandBerlekamp: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				big_int p;
				big_int exp;
				char c;
				if (!isVector(args, v) || !(args >> p))
					throw 1; //throw new ParseError();
				if ( (args.str())[0] != ')'){
					if ((args >> exp) && !(args >> c && c ==')'))
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
	};
	class CommandCantorZassenhaus: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				big_int p;
				big_int exp;
				char c;
				if (!isVector(args, v) || !(args >> p))
					throw 1; //throw new ParseError();
				if ( (args.str())[0] != ')'){
					if ((args >> exp) && !(args >> c && c ==')'))
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
	};
	class CommandHensel: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				char c;
				if ( !isVector(args, v) || !((args >> c) && c == ')') )
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
	};
	class CommandModularGCD: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v1, v2;
				char c;
				if (!isVector(args, v1) || !isVector(args, v2)
						|| !((args >> c) && c == ')'))
					throw 1; //throw new ParseError();
				Zxelem_b pol1(v1), pol2(v2);
				std::cout << modularGCD(pol1, pol2) << std::endl;
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
	};
	class CommandCRA: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				big_int aux;
				std::vector<big_int> m, u;
				while (args >> aux ){
					m.push_back(aux);
					if (!(args >> aux))
						throw 1; //throw new ParseError();
					u.push_back(aux);
				}
				std::cout << integerCRA(m, u) << std::endl;
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
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
	};
	class CommandPolardFactor: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				big_int aux;
				char c;
				if (!(args >> aux) || !(args >> c) || c!= ')')
					throw 1; //throw new ParseError();
				//TODO	//std::cout << factorizationPollardRhoBrent(aux);
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
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
	};
	class CommandMillerRabin: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				big_int num;
				char c;
				if (!(args >> num) || !(args >> c) || c != ')')
					throw 1; //throw new ParseError();
				std::cout << num << "is ";
				if(!millerRabin(num))
					std::cout << "not ";
				std::cout << "prime" << std::endl;
			}catch(...){
				std::cout << "Parse error" << std::endl;
			}
		}
	};
	class CommandIrrGFp: public Command {
	public:
		void parseAndRun(std::stringstream & args){
			try{
				std::vector<big_int> v;
				big_int p;
				if (!isVector(args, v) || !(args >> p))
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
	};

	UserInterface::UserInterface(){
		cmds["factorBerlekamp"] = CommandBerlekamp();
		cmds["factorCZ"] = CommandCantorZassenhaus();
		cmds["factorHensel"] = CommandHensel();
		cmds["modularGCD"] = CommandModularGCD();
		cmds["chineseRA"] = CommandCRA();
		cmds["eea_ed"] = CommandEEA_ED();
		cmds["factorInteger"] = CommandPolardFactor();
		cmds["discreteLog"] = CommandPolardLog();
		cmds["isPrime"] = CommandMillerRabin();
		cmds["isIrreducibleGFp"] = CommandIrrGFp();
		//TODO: Preguntar a Mario por el Russian peasant
		//cmds["O"] = new CommandMillerRabin();
	}
	void UserInterface::run(){
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
			if ( it == cmds.end() || cmdline[cmdline.size()-1] != ')'){
				std::cout << "Unrecognized option";
				this->help();
			}
			else{
				(it->second).parseAndRun(ss);
			}
		}
	}
	void UserInterface::help(){
		std::cout << "esto es una ayuda" << std::endl;

	}
}
