#ifndef __USER_INTERFACE
#define __USER_INTERFACE

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

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

namespace alcp {

	class UserInterface {//This is a singleton
	private: class Command;
	public:

		UserInterface();
		static UserInterface& instance();
		bool isCommand (const std::string & s);
		void run();
		void help();
		void callHelp(const std::string & s);
	private:
		std::map< std::string, std::unique_ptr<Command> > cmds;

///////////////////////////////////////////////////////////////////////////////
// Commands
///////////////////////////////////////////////////////////////////////////////
		class Command {
			public:
				virtual void parseAndRun(std::istringstream & args) = 0;
				virtual void help (const std::string & name) = 0;
		};

		class CommandHelp: public Command{
			public:
				CommandHelp() = default;
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandBerlekamp: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandCantorZassenhaus: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandHensel: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandModularGCD: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandCRA: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandEEA_ED: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandPollardFactor: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandPollardLog: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandMillerRabin: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};

		class CommandIrrGFp: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};
		class CommandBCH: public Command{
			public:
				void parseAndRun(std::istringstream & args) override;
				void help (const std::string & name) override;
		};
	};
}
#endif // __USER_INTERFACE

