#include "Tensor.hpp"
#include "Parameters_BHZ.hpp"
#include "Connection_BHZ.hpp"

#ifndef Hamiltonian_BHZ_HPP
#define Hamiltonian_BHZ_HPP

class Hamiltonian_BHZ
{
public:
    Hamiltonian_BHZ(Parameters_BHZ& Parameters_BHZ__, Connection_BHZ& Connection_BHZ__):
     Parameters_BHZ_(Parameters_BHZ__),Connection_BHZ_(Connection_BHZ__){
        Initialize();
    }

    void Initialize();
    void Diagonalizer(char option);
    double FermiFunction(double en_, double mu_);
    double ChemicalPotential(double muin_, double particles_);
    
    

    Mat_2_Complex_doub Ham_,Evecs_;
    Mat_1_doub Evals_;

    Parameters_BHZ &Parameters_BHZ_;
    Connection_BHZ &Connection_BHZ_;

};


#endif