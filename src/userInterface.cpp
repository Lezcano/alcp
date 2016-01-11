#include "userInterface.hpp"
#include <map>
#include <string>
#include <sstream>

UserInterface::UserInterface(){
	cmds["factorBerlekamp"] = new ComandBerlekamp();
	cmds["factorCZ"] = new ComandCantorZassenhaus();
	cmds["factorHensel"] = new ComandHensel();
	cmds["modularGCD"] = new ComandModularGCD();
	cmds["chineseRA"] = new ComandCRA();
	cmds["eea_ed"] = new ComandEEA_ED();
	cmds["factorInteger"] = new ComandPolardFactor();
	cmds["discreteLog"] = new ComandPolardLog();
	cmds["isPrime"] = new ComandMillerRabin();
	cmds["isIrreducibleGFp"] = new ComandIrrGFp();
	//TODO: Preguntar a Mario por el Russian peasant
	//cmds["O"] = new ComandMillerRabin();
}
void run(){
	string cmdline;
	cout << "Computer algebra system by Mario Lezcano and David Martinez\n"
	cout << "Write \"help\" for help\n";
	cout << ">> ";
	getline (cin, cmdline);
	while (cmd != "quit"){
		string cmd = cmdline.substr(0, cmdline.find('('));
		auto it = cmds.find(cmd); 
		if ( it == cmds.end() || cmdline[cmdline.size()-1] != ')'){
			cout << "Unrecognized option";
			this.help();
		}
		else{
			it	
		}
	}
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
	std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

class Comand {
	public:
	private: 
}


class ComandBerlekamp: public Comand {
public:
	run(vector<string> args){
		try{
			if (args.size() != 2  || args.size() != 3)
				throw new ErrIncorrectNumOfArgs();
			vector<big_int> v;
			if (!isVector(args[0], v))
				throw new ParseError();
			std::vector<std::pair<Fxelem, unsigned int> > factors;
			unsigned int exp;
			if (args.size() == 2 || (exp = stoi(args[2])) == 1){
				Fp f(stoi(args[1]));
				vector<Fpelem> vf;
				for (unsigned int i = 0; i < v.size(); i++){
					vf.push_back(f.get(v[i]));	
				}
				Fpxelem pol(vf);
				factors = factorizationBerlekamp(pol);
			}
			else{
				Fq f(stoi(args[1]));
				vector<Fqelem> vf;
				for (unsigned int i = 0; i < v.size(); i++){
					vf.push_back(f.get(v[i]));	
				}
				Fqxelem pol(vf);
				factors = factorizationBerlekamp(pol);
			}
 			cout << "Factors:" << endl;
			for (auto &pair: factors){
				if (pair.second != 1)
					cout << "( ";
				cout << pair.first;
				if (pair.second != 1) 
					cout << ")^" << pair.second;
				cout << endl << endl;
			}
		}catch(Exception e){
			cout << e;
		}
	}
}
class ComandCantorZassenhaus: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandHensel: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandModularGCD: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandCRA: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandEEA_ED: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandPolardFactor: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandPolardLog: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandMillerRabin: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class ComandIrrGFp: public Comand {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
