#ifndef __USER_INTERFACE
#define __USER_INTERFACE

#include <iostream>
#include <map>
#include <string>
#include <sstream>

namespace alcp {
	class Command {
	public:
		virtual void parseAndRun(std::stringstream & args){}
		virtual void help (string & name){}
	private:
	};

//TODO: This should be a singleton,
	class UserInterface {
	public:
		UserInterface();
		static void run();
		static void help();
		static bool isCommand (const string & s);
		static Command && getCommand(const string & s);
	private:
		static std::map< std::string, Command> cmds;
}
#endif // __FACTORIZATION_FQ

