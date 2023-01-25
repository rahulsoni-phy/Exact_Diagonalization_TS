#include <complex>
#include <Eigen/Dense>
#include "Parameters_BHZ.hpp"

#ifndef Connection_BHZ_HPP
#define Connection_BHZ_HPP

using namespace std;
using namespace Eigen;

class Connection_BHZ
{
public:
    Connection_BHZ(Parameters_BHZ& Parameters_BHZ__):Parameters_BHZ_(Parameters_BHZ__){
        Initialize();
        ConnectionMatrix();
    }
    int lx_, ly_, orbs_, spin_,size_;
    double A_, B_, M_;
    complex<double> Zero_Complex, One_Complex, Iota_Complex;

    MatrixXcd C_mat;

    bool Periodic_X,Periodic_Y;

    Parameters_BHZ &Parameters_BHZ_;
    
    void Initialize();
    void ConnectionMatrix();
    void PrintConnection();
    
};


#endif