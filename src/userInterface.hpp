#ifndef __USER_INTERFACE
#define __USER_INTERFACE

#include "command.hpp"
#include <iostream>
#include <map>
#include <string>
#include <sstream>

namespace alcp {

//TODO: This should be a singleton,
	class UserInterface {
	public:
		UserInterface();
		static void run();
		static void help();
		static bool isCommand (const std::string & s);
		static Command && getCommand(const std::string & s);
	private:
		static std::map< std::string, Command> cmds = {
				std::make_pair<std::string, Command>(
						"help", CommandHelp()),
				std::make_pair<std::string, Command>(
						"factorBerlekamp", CommandBerlekamp()),
				std::make_pair<std::string, Command>(
						"factorCZ", CommandCantorZassenhaus()),
				std::make_pair<std::string, Command>(
						"factorHensel", CommandHensel()),
				std::make_pair<std::string, Command>(
						"modularGCD", CommandModularGCD()),
				std::make_pair<std::string, Command>(
						"chineseRA", CommandCRA()),
				std::make_pair<std::string, Command>(
						"eea_ed", CommandEEA_ED()),
				std::make_pair<std::string, Command>(
						"factorInteger", CommandPolardFactor()),
				std::make_pair<std::string, Command>(
						"discreteLog", CommandPolardLog()),
				std::make_pair<std::string, Command>(
						"isPrime", CommandMillerRabin()),
				std::make_pair<std::string, Command>(
						"isIrreducibleGFp", CommandIrrGFp())
			};
			//cmds[""] = new CommandMillerRabin();
};
#endif // __FACTORIZATION_FQ

