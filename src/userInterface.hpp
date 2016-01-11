#ifndef __USER_INTERFACE
#define __USER_INTERFACE

#include <iostream>
#include <map>


class UserInterface {
	public:
		void run();	
	private:
		static map<string, Comands> cmds;

}

#endif // __FACTORIZATION_FQ

