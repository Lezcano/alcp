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
	private:
	};

	class UserInterface {
	public:
		UserInterface();
		void run();
		void help();
	private:
		static std::map< std::string, Command> cmds;
	};
}
#endif // __FACTORIZATION_FQ

