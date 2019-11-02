/* 
 * File:   main.cpp
 * Author: alex
 *
 * Created on October 31, 2019, 7:37 PM
 */

#include <iostream>

#include "QCC.h"

int main(int argc, char** argv) {

	QCC* comp = new QCC();

	std::cout << "Test " << comp->compile() << std::flush;

	return 0;
}

