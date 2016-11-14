#include "userInterface.hpp"

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include <cstdarg>
#include <algorithm>
#include <utility>
#include <cstring>

#include "types.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "fqelem.hpp"
#include "zxelem.hpp"
#include "exceptions.hpp"
#include "factorizationFq.hpp"
#include "generalPurpose.hpp"
#include "bchCodes.hpp"

namespace alcp {
	void asdf(){
		while (true){
			bool salir = false;
			std::cout << "Choose an example (write its number):" << std::endl;
			std::cout << "1) p = 2, prim. pol. of degree 4, l = 15, c = 1, d = 5 " << std::endl;
			std::cout << "2) p = 3, prim. pol. of degree 2, l = 8, c = 3, d = 3 " << std::endl;
			std::cout << "3) p = 7, prim. pol. of degree 2, l = 48, c = 3, d = 9 " << std::endl;
			std::cout << "4) p = 7, prim. pol. of degree 3, l = 342, c = 1, d = 11 " << std::endl;
			std::cout << "5) p = 3, prim. pol. of degree 6, l = 728, c = 146, d = 15 " << std::endl;
			std::cout << "6) p = 7, prim. pol. of degree 4, l = 2400, c = 1, d = 31 " << std::endl;
			int option;
			std::cin >> option;
			char asdf;
			std::cin.get(asdf);
			size_t l, c, d;
			std::vector<big_int> v;
			big_int p;
			switch(option){
			case 1:
				l=15;
				c=1;
				d=5;
				v={1, 1, 0, 0, 1};
				p = 2;
				break;
			case 2:
				l=8;
				c=3;
				d=3;
				v={1, 1, 2};
				p = 3;
							break;
			case 3:
				l=48;
				c=3;
				d=9;
				v={5, 4, 1};
				p = 7;
							break;
			case 4:
				l=342;
				c=1;
				d=11;
				v={4, 4, 6, 1};
				p = 7;
							break;
			case 5:
				l=728;
				c=146;
				d=15;
				v={2, 1, 0, 0, 0, 0, 1 };
				p = 3;
							break;
			case 6:
				l=2400;
				c=1;
				d=31;
				v={3, 0, 4, 2, 1};
				p = 7;
							break;
			default: salir = true;
			}
			if (salir)
				break;
		//You can get primitive polynomials here: http://fchabaud.free.fr/English/default.php?COUNT=3&FILE0=Poly&FILE1=GF(7)&FILE2=Primitive

		Fpxelem_b prim_pol(v, p);
		BCH bch = BCH(prim_pol, 1, l, c, d);
		//std::cout << "Give me the message(Intro for random message): " << std::endl;
		//std::string cmdline;
		//std::cout << ">> ";
		//getline(std::cin, cmdline);
			while(true){
				Fpxelem_b message;
				message = randomPol(Fp_b(p), bch.getDimension()-1);
				std::cout << "The random message chosen is" << std::endl;
				std::cout << "MESSAGE := " << message << std::endl;
				std::cout << std::endl <<"The codification of your message is:" << std::endl;
				Fpxelem_b sent = bch.encode(message);
				std::cout << "CODE SENT:= " << sent << std::endl;
				std::cout << "Sending your message through an error-prone channel" << std::endl;
				std::set<int> error_indices;
				Fpxelem_b received;
				std::tie(error_indices, received) = bch.randomErrors(sent);
				std::cout << "An oracle tells us that " << error_indices.size() <<
						" error/s ocurred during the transmission at index/indices: ";
				for(auto & elem: error_indices){
					std::cout << elem << ", ";
				}
				std::cout << "The code received is" << std::endl;
				std::cout << "CODE RECEIVED:= " << received << std::endl;


				Fpxelem_b decoded = bch.decode(received);
				std::cout << "The corrected code is " << std::endl;
				std::cout << "CODE CORRECTED := " << decoded << std::endl << std::endl;
				if (sent == decoded){
					std::cout << "I have checked it is the same code that was sent" << std:: endl;
					std::cout << "The original message is " << decoded/bch.getG()<< std:: endl;
				}
				else{
					std::cout << "OH! YOUR HAVE A MISTAKE. YOUR CODE HAS NOT WORKED." << std::endl;
					while (true)
						std::cout << "asdf";
				}
				char c;
				std::cin.get(c);
				std::cout << std:: endl << std:: endl << std:: endl << std:: endl << std:: endl;
				if (c=='o')
					break;
			}
		}
	}
}
int main () {
	alcp::asdf();
    //alcp::UserInterface::instance().run();

	return 0;
}

