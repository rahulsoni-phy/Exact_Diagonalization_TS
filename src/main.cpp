#include <iostream> //for basic I/O services! (specifically: cin, cout etc.)
#include <iomanip>  //for setting precision of the results!
#include <fstream>  //to read and write from/to files!
#include <stdio.h>  //prints the strings, integer, character etc. on the output screen
#include <vector>   //for arrays that can change size!
#include <cstdlib>  //defines standard C++ functions and macros! (like: atof, atoi, abs, rand etc.)
#include <string>   //for objects that can be represented as stream of characters!
#include <stdexcept>//defines a set of standard exceptions and report common errors!
#include <random>   //for random numbers!
#include <complex>  //for complex numbers!
#include <cmath>    //for common mathematical functions!
#include <cassert>  //error handling library that aborts the code if defined condition are not met!

extern "C" {
    #include <lapacke.h>
    #include <cblas.h>
}

#include "Tensor.hpp"
#include "Parameters_BHZ.hpp"
#include "Connection_BHZ.hpp"
#include "Hamiltonian_BHZ.hpp"

using namespace std;

int main(int argc, char *argv[]){

    string Model_type = argv[1];
    string input_file = argv[2];

    if(Model_type=="BHZ"){
        Parameters_BHZ Parameters_BHZ_;
        Parameters_BHZ_.Initialize(input_file);

        Connection_BHZ Connection_BHZ_(Parameters_BHZ_);
        Connection_BHZ_.ConnectionMatrix();
    }

    return 0;
}
