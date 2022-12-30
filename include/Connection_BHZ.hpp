
#include "Tensor.hpp"
#include "Parameters_BHZ.hpp"

#ifndef Connection_BHZ_HPP
#define Connection_BHZ_HPP

class Connection_BHZ
{
public:
    Connection_BHZ(Parameters_BHZ& Parameters_BHZ__):Parameters_BHZ_(Parameters_BHZ__){
        Initialize();
        ConnectionMatrix();
    }
    int lx_, ly_, orbs_, spin_,size_;
    double A_, B_, M_;
    Mat_2_Complex_doub C_mat;

    bool Periodic_X,Periodic_Y;

    Parameters_BHZ &Parameters_BHZ_;
    
    void Initialize();
    void ConnectionMatrix();
};


#endif