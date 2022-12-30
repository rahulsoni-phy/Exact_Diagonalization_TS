#include "Tensor.hpp"
#include "Parameters_BHZ.hpp"
#include "Connection_BHZ.hpp"

void Connection_BHZ::Initialize(){
    lx_=Parameters_BHZ_.Lx;
    ly_=Parameters_BHZ_.Ly;

    A_=Parameters_BHZ_.A_val;
    B_=Parameters_BHZ_.B_val;
    M_=Parameters_BHZ_.M_val;

    orbs_=Parameters_BHZ_.N_Orbs;
    spin_=2;

    if(Parameters_BHZ_.PBC_X==true){
        Periodic_X=true;
    }
    else{
        Periodic_X=false;
    }

    if(Parameters_BHZ_.PBC_Y==true){
        Periodic_Y=true;
    }
    else{
        Periodic_Y=false;
    }

    size_=Parameters_BHZ_.Ham_Size;

    C_mat.resize(size_);
    for(int i=0;i<size_;i++){
        C_mat[i].resize(size_);
        for(int j=0;j<size_;j++){
            C_mat[i][j]=Zero_Complex;
        }
    }

}

void Connection_BHZ::ConnectionMatrix(){

    int r, r0, r1, r2;
    int half_size_= (int) (size_/2);

    for(int spin=0;spin<spin_;spin++){
        for(int orb=0;orb<orbs_;orb++){
            for(int rx=0;rx<lx_;rx++){
                for(int ry=0;ry<ly_;ry++){

                    r = spin*half_size_ + (orb + 2*ry + 2*ly_*rx);
                    r1= spin*half_size_ + (orb + 2*ry + 2*ly_*(rx+1));
                    r2= spin*half_size_ + (orb + 2*(ry+1) + 2*ly_*rx);

                    //Adding onsite mass term and NN hopping connections:---->
                    C_mat[r][r] = 1.0*( pow(-1.0, 1.0*orb) )*(M_-4*B_)*One_Complex;

                    if(rx!=lx_-1){
                        C_mat[r][r1] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                        C_mat[r1][r] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                    }
                    if(ry!=ly_-1){
                        C_mat[r][r2] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                        C_mat[r2][r] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                    }

                    if(Periodic_X==true){
                        r0 = spin*half_size_ + (orb + 2*ry);
                        if(rx==lx_-1){
                            C_mat[r][r0] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                            C_mat[r0][r] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                        }
                    }
                    if(Periodic_Y==true){
                        r0 = spin*half_size_ + (orb + 2*ly_*rx);
                        if(ry==ly_-1){
                            C_mat[r][r0] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                            C_mat[r0][r] = 1.0*( pow(-1.0, 1.0*orb) )*B_*One_Complex;
                        }
                    }

                }
            }
        }
    }

    //--------------------------------------------------//
    //For sigma=up:
    //H(p,R+x;s,R)=H(s,R+x;p,R)=iA/2
    //H(p,R+y;s,R)=-A/2, H(s,R+y;p,R)=A/2
    //
    //For sigma=dn:
    //H(p,R+x;s,R)=H(s,R+x;p,R)=-iA/2
    //H(p,R+y;s,R)=-A/2, H(s,R+y;p,R)=A/2
    //--------------------------------------------------//

    int r1s,r1p,r2s,r2p;

    for(int spin=0;spin<spin_;spin++){
        for(int rx=0;rx<lx_;rx++){
            for(int ry=0;ry<ly_;ry++){
                for(int orb=0;orb<orbs_;orb++){

                    r = spin*half_size_ + (orb + 2*ry + 2*ly_*rx);
                    if(orb==0){
                        r1 = spin*half_size_ + (1 + 2*ry + 2*ly_*(rx+1));
                        r2 = spin*half_size_ + (1 + 2*(ry+1) + 2*ly_*rx);
                    }
                    if(orb==1){
                        r1 = spin*half_size_ + (0 + 2*ry + 2*ly_*(rx+1));
                        r2 = spin*half_size_ + (0 + 2*(ry+1) + 2*ly_*rx);
                    }

                    //R->R+x : ---------->
                    if(rx!=lx_-1){
                        if(spin==0){
                            C_mat[r1][r] = -(1.0*A_/2.0)*Iota_Complex;
                            C_mat[r][r1] = (1.0*A_/2.0)*Iota_Complex;
                        }
                        if(spin==1){
                            C_mat[r1][r] = (1.0*A_/2.0)*Iota_Complex;
                            C_mat[r][r1] = -(1.0*A_/2.0)*Iota_Complex;
                        }
                    }

                    if(Periodic_X==true && lx_>2){ //For 2 sites system PBC doesn't make sense here
                        if(rx==lx_-1){
                            if(orb==0){ r0 = spin*half_size_ + (1 + 2*ry); }
                            if(orb==1){ r0 = spin*half_size_ + (0 + 2*ry); }

                            if(spin==0){
                                C_mat[r0][r] = -(1.0*A_/2.0)*Iota_Complex;
                                C_mat[r][r0] = (1.0*A_/2.0)*Iota_Complex;
                            }
                            if(spin==1){
                                C_mat[r0][r] = (1.0*A_/2.0)*Iota_Complex;
                                C_mat[r][r0] = -(1.0*A_/2.0)*Iota_Complex;
                            }
                        }
                    }
                    
                    //R->R+y : ---------->
                    if(ry!=ly_-1){
                        if(orb==0){
                            C_mat[r2][r] = (1.0*A_/2.0)*One_Complex;
                            C_mat[r][r2] = (1.0*A_/2.0)*One_Complex;
                        }
                        if(orb==1){
                            C_mat[r2][r] = -(1.0*A_/2.0)*One_Complex;
                            C_mat[r][r2] = -(1.0*A_/2.0)*One_Complex;
                        }
                    }

                    if(Periodic_Y==true && ly_>=2){
                        if(ry==ly_-1){
                            
                            if(orb==0){ r0 = spin*half_size_ + (1 + 2*ly_*rx); }
                            if(orb==1){ r0 = spin*half_size_ + (0 + 2*ly_*rx); }

                            if(orb==0){
                                C_mat[r0][r] = (1.0*A_/2.0)*One_Complex;
                                C_mat[r][r0] = (1.0*A_/2.0)*One_Complex;
                            }
                            if(orb==1){
                                C_mat[r0][r] = -(1.0*A_/2.0)*One_Complex;
                                C_mat[r][r0] = -(1.0*A_/2.0)*One_Complex;
                            }
                        }
                    }

                }
            }
        }
    }

}