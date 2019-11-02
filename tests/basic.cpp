/**
 * File:   newsimpletest.cpp
 * Author: alex
 *
 * Created on November 1, 2019, 11:37 AM
 */

#include <stdlib.h>
#include <iostream>
#include <cassert>

#include "QCC.h"

void test1() {
	QCC* comp = new QCC();
	assert(comp->compile() == 1234);
}

int main(int argc, char** argv) {
    test1();

    return (EXIT_SUCCESS);
}


