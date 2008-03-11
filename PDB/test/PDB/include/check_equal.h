/*
 *  check_equal.h
 *  PDB-xcode
 *
 *  Created by Daniel Russel on 2/7/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <cstring> // for std::strcmp

int return_code__=EXIT_SUCCESS;

void check_equal(std::istream &a, std::istream &b) {
	do {
		char bufa[1000];
		char bufb[1000];
		a.getline(bufa, 1000);
		b.getline(bufb, 1000);
		if (!a && !b) {
			break;
		} else if (!a) {
			std::cerr << "Source missing: " << bufb << std::endl;
			return_code__=EXIT_FAILURE;
		} else if (!b) {
			std::cerr << "Target missing: " << bufa << std::endl;
			return_code__=EXIT_FAILURE;		
		} else {
			if (std::strcmp(bufa, bufb) != 0) {
				std::cerr << "Lines are not equal. They are \n";
				std::cerr << bufa << std::endl;
				std::cerr << bufb << std::endl;
				return_code__=EXIT_FAILURE;
			}
		}
	} while (true);
}

