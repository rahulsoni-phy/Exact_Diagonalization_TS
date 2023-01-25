#include <Eigen/Dense>
#include "Parameters_BHZ.hpp"
#include "Connection_BHZ.hpp"
#include "Hamiltonian_BHZ.hpp"

#ifndef Observables_BHZ_HPP
#define Observables_BHZ_HPP

using namespace Eigen;

class Observables_BHZ
{
public:
    Observables_BHZ(Parameters_BHZ& Parameters_BHZ__, Connection_BHZ& Connection_BHZ__, Hamiltonian_BHZ& Hamiltonian_BHZ__):
    Parameters_BHZ_(Parameters_BHZ__),Connection_BHZ_(Connection_BHZ__),Hamiltonian_BHZ_(Hamiltonian_BHZ__){
        Initialize();
    }
    int lx_, ly_, orbs_, spin_;
    int size_,cells_, mhs_;

    double mu_val;
    int w_size;
    double dw,eta;
    
    VectorXd evals_;
    MatrixXcd evecs_;

    complex<double> Zero_Complex, One_Complex, Iota_Complex;

    void Initialize();
    void Calculate_Local_Density_of_Electrons();
    void Calculate_MomSpace_Occupation_Number();
    void Calculate_Energy_Bands_on_Path();
    void Calculate_Spectral_Function();
    void Calculate_Spin_Chern_Number();

    Parameters_BHZ &Parameters_BHZ_;
    Connection_BHZ &Connection_BHZ_;
    Hamiltonian_BHZ &Hamiltonian_BHZ_;
    
};


#endif