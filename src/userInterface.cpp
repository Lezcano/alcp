#include "userInterface.hpp"
#include <map>
#include <string>
#include <sstream>

UserInterface::UserInterface(){
	cmds["factorBerlekamp"] = new CommandBerlekamp();
	cmds["factorCZ"] = new CommandCantorZassenhaus();
	cmds["factorHensel"] = new CommandHensel();
	cmds["modularGCD"] = new CommandModularGCD();
	cmds["chineseRA"] = new CommandCRA();
	cmds["eea_ed"] = new CommandEEA_ED();
	cmds["factorInteger"] = new CommandPolardFactor();
	cmds["discreteLog"] = new CommandPolardLog();
	cmds["isPrime"] = new CommandMillerRabin();
	cmds["isIrreducibleGFp"] = new CommandIrrGFp();
	//TODO: Preguntar a Mario por el Russian peasant
	//cmds["O"] = new CommandMillerRabin();
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

class Command {
	public:
	private: 
}


class CommandBerlekamp: public Command {
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
class CommandCantorZassenhaus: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandHensel: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandModularGCD: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandCRA: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandEEA_ED: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandPolardFactor: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandPolardLog: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandMillerRabin: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
class CommandIrrGFp: public Command {
public:
	run(vector<string> args){
		if (args.size() != )
			errIncorrectNumOfArgs();
		else{
		
		}

	}
}
