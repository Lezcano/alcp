#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace alcp {

	bool coma(std::stringstream & args);
	bool isVector(std::stringstream & args, std::vector<big_int> & v);
	bool parenCl(std::stringstream & args);



	UserInterface::UserInterface() = default;
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
