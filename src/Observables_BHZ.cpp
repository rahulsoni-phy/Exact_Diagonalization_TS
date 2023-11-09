#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Dense>

#include "Parameters_BHZ.hpp"
#include "Connection_BHZ.hpp"
#include "Hamiltonian_BHZ.hpp"
#include "Observables_BHZ.hpp"
#define PI 3.14159265359

using namespace std;
using namespace Eigen;

void Observables_BHZ::Initialize(){
    lx_ = Parameters_BHZ_.Lx;
    ly_ = Parameters_BHZ_.Ly;
    orbs_ = Parameters_BHZ_.N_Orbs;
    spin_ = 2;
    cells_ = Parameters_BHZ_.Total_Cells;

    complex<double> One_Complex(1.0,0.0),Zero_Complex(0.0,0.0),Iota_Complex(0.0,1.0);

    evals_=Hamiltonian_BHZ_.Evals_;
    evecs_=Hamiltonian_BHZ_.Evecs_;

    mu_val = Hamiltonian_BHZ_.ChemicalPotential(Parameters_BHZ_.Total_Particles);

    dw = Parameters_BHZ_.dw_dos;
    eta = Parameters_BHZ_.eta_dos;

    w_size = (int) ( (Parameters_BHZ_.w_max - Parameters_BHZ_.w_min)/dw );

}


void Observables_BHZ::Calculate_Local_Density_of_Electrons(){

    Mat_4_Complex_doub CD_local;
    CD_local.resize(lx_);
    for(int rx=0;rx<lx_;rx++){
        CD_local[rx].resize(ly_);

        for(int ry=0;ry<ly_;ry++){
            CD_local[rx][ry].resize(orbs_);

            for(int orb=0;orb<orbs_;orb++){
                CD_local[rx][ry][orb].resize(spin_);

                for(int s=0;s<spin_;s++){
                    CD_local[rx][ry][orb][s]=Zero_Complex;
                }

            }

        }

    }

    int mp;
    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            for(int orb=0;orb<orbs_;orb++){
                for(int spin=0;spin<spin_;spin++){

                    mp = spin*mhs_ + orb + 2*ry + 2*ly_*rx;

                    for(int n=0;n<size_;n++){
                        
                        CD_local[rx][ry][orb][spin] += (conj(evecs_(n,mp)))*((evecs_(n,mp)))*( Hamiltonian_BHZ_.FermiFunction(evals_(n),mu_val) );
                    }
                }
            }
        }
    }

    string file_Avg_LDOE="Avg_LDOE_along_x.txt";
    ofstream file_Avg_CD_out(file_Avg_LDOE.c_str());

    file_Avg_CD_out<<"#rx"<<"   "<<"Orb-0,Spin-0"<<"    "<<"Orb-0,Spin-1"<<"    "<<"Orb-1,Spin-0"<<"    "<<"Orb-1,Spin-1"<<endl;

    complex<double> CD_s_up,CD_s_dn,CD_p_up,CD_p_dn;
    for(int rx=0;rx<lx_;rx++){
        CD_s_up = Zero_Complex;
        CD_s_dn = Zero_Complex;
        CD_p_up = Zero_Complex;
        CD_p_dn = Zero_Complex;
        for(int ry=0;ry<ly_;ry++){
            CD_s_up += (1.0/(1.0*ly_))*CD_local[rx][ry][0][0];
            CD_s_dn += (1.0/(1.0*ly_))*CD_local[rx][ry][0][1];
            CD_p_up += (1.0/(1.0*ly_))*CD_local[rx][ry][1][0];
            CD_p_dn += (1.0/(1.0*ly_))*CD_local[rx][ry][1][1];
        }
        file_Avg_CD_out<<rx<<"  "<<CD_s_up.real()<<"    "<<CD_s_dn.real()<<"    "<<CD_p_up.real()<<"    "<<CD_p_dn.real()<<endl;
    }
     
}


void Observables_BHZ::Calculate_Density_of_States(){

    string file_dos="Density_of_states.txt";
    ofstream file_dos_out(file_dos.c_str());

    double omega=Parameters_BHZ_.w_min;
    double dos_val;

    while(omega <= Parameters_BHZ_.w_max){
        dos_val=0.0;
        file_dos_out<< omega <<"    ";

        for(int n=0;n<size_;n++){
            dos_val += (1.0/(size_*PI))*(eta/((omega-evals_(n))*(omega-evals_(n)) + eta*eta));
        }

        file_dos_out<<dos_val<<endl;
        omega = omega + dw;
    }    
}


void Observables_BHZ::Calculate_Bmat(){
    B_mat.resize(cells_);
    for(int r1=0;r1<cells_;r1++){
        B_mat[r1].resize(cells_);

        for(int r2=0;r2<cells_;r2++){
            B_mat[r1][r2].resize(w_size);

            for(int om=0;om<w_size;om++){
                B_mat[r1][r2][om] = Zero_Complex;
            }
        }
    }

    int m1,m2;
    omega=0.0;

    for(int r1=0;r1<cells_;r1++){
        r1x = r1/ly_;
        r1y = r1%ly_;

        for(int r2=0;r2<cells_;r2++){
            r2x = r2/ly_;
            r2y = r2%ly_;

            for(int om=0;om<w_size;om++){
                omega = Parameters_BHZ_.w_min + om*dw;

                for(int orb=0;orb<orbs_;orb++){
                    for(int s=0;s<spin_;s++){
                        m1 = s*mhs_ + orb + 2*r1y + 2*ly_*r1x;
                        m2 = s*mhs_ + orb + 2*r2y + 2*ly_*r2x;

                        for(int n=0;n<size_;n++){
                            B_mat[r1][r2][om] += (1.0/PI)*( ( conj(evecs_(n,m1)) ) * (evecs_(n,m2)) )*( 
                                (eta)/((omega-evals_(n))*(omega-evals_(n))+(eta*eta)) );
                        }
                    }
                }
            }

        }

    }
}


void Observables_BHZ::Calculate_MomSpace_Occupation_Number(){

    Mat_2_Complex_doub P_r_rp;
    P_r_rp.resize(cells_);
    for(int c1=0;c1<cells_;c1++){
        P_r_rp[c1].resize(cells_);
        for(int c2=0;c2<cells_;c2++){
            P_r_rp[c1][c2]=Zero_Complex;
        }
    }

    int m,mp;
    int r1,r2;
    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){

            r1=iy + ly_*ix;

            for(int ixp=0;ixp<lx_;ixp++){
                for(int iyp=0;iyp<ly_;iyp++){

                    r2=iyp + ly_*ixp;

                    for(int orb=0;orb<orbs_;orb++){
                        for(int s=0;s<spin_;s++){

                            m = s*mhs_ + orb + 2*r1;
                            mp = s*mhs_ + orb + 2*r2;

                            for(int n=0;n<size_;n++){

                                P_r_rp[r1][r2] += (1.0/(1.0*cells_))*(conj(evecs_(n,m)))*((evecs_(n,mp)))*( 
                                    Hamiltonian_BHZ_.FermiFunction(evals_(n),mu_val) );

                            }

                        }
                    }

                }
            }

        }
    }

    complex<double> Nk;
    double kx,ky;
    string file_mom_space_occ="Mom_space_avg_occupation_distribution.txt";
    ofstream file_mom_space_occ_out(file_mom_space_occ.c_str());

    for(int kx_ind=-lx_/2;kx_ind<=lx_/2;kx_ind++){
        kx=(2.0*PI*kx_ind)/(1.0*lx_);

        for(int ky_ind=-ly_/2;ky_ind<=ly_/2;ky_ind++){
            ky=(2.0*PI*ky_ind)/(1.0*ly_);
            Nk=Zero_Complex;

            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    r1 = iy + ly_*ix;

                    for(int ixp=0;ixp<lx_;ixp++){
                        for(int iyp=0;iyp<ly_;iyp++){
                            r2 = iyp + ly_*ixp;

                            Nk += exp( Iota_Complex*(kx*1.0*(ix-ixp) + ky*1.0*(iy-iyp)) )*P_r_rp[r1][r2];

                        }
                    }

                }
            }
            file_mom_space_occ_out<<kx<<"   "<<ky<<"    "<<Nk.real()<<endl;
        }
        file_mom_space_occ_out<<endl;
    }


}


void Observables_BHZ::Calculate_Energy_Bands_on_Path(){
    
    string file_Akw="Akw_on_path.txt";
    ofstream file_Akw_out(file_Akw.c_str());

    complex<double> Akw_;
    omega=0.0;

    //Path Gamma to X--->
    ky_ind=0;
    k_ind=0;
    for(kx_ind=0;kx_ind<=lx_/2;kx_ind++){
        kx = (2.0*PI*kx_ind)/(1.0*lx_);
        ky = (2.0*PI*ky_ind)/(1.0*ly_);

        for(int om=0;om<w_size;om++){
            omega = Parameters_BHZ_.w_min + om*dw;
            Akw_ = Zero_Complex;

            for(int r1=0;r1<cells_;r1++){
                r1x = r1/ly_;       r1y = r1%ly_;
                for(int r2=0;r2<cells_;r2++){
                    r2x = r2/ly_;       r2y = r2%ly_;

                    Akw_ += (1.0/(4.0*cells_))*( exp(Iota_Complex*kx*(1.0*(r1x-r2x)))*exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*
                            B_mat[r1][r2][om] );
                }
            }
            file_Akw_out<<k_ind<<"      "<<omega<<"     "<<Akw_.real()<<endl;
        }
        file_Akw_out<<endl;
        k_ind = k_ind + 1;
    }

    //Path X to M --->
    kx_ind=lx_/2;
    for(ky_ind=0;ky_ind<=ly_/2;ky_ind++){
        kx = (2.0*PI*kx_ind)/(1.0*lx_);
        ky = (2.0*PI*ky_ind)/(1.0*ly_);

        for(int om=0;om<w_size;om++){
            omega = Parameters_BHZ_.w_min + om*dw;
            Akw_ = Zero_Complex;

            for(int r1=0;r1<cells_;r1++){
                r1x = r1/ly_;       r1y = r1%ly_;
                for(int r2=0;r2<cells_;r2++){
                    r2x = r2/ly_;       r2y = r2%ly_;

                    Akw_ += (1.0/(4.0*cells_))*( exp(Iota_Complex*kx*(1.0*(r1x-r2x)))*exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*
                            B_mat[r1][r2][om] );
                }
            }
            file_Akw_out<<k_ind<<"      "<<omega<<"     "<<Akw_.real()<<endl;
        }
        file_Akw_out<<endl;
        k_ind = k_ind + 1;
    }

    //Path M to Gamma --->
    kx_ind=lx_/2-1;
    ky_ind=ly_/2-1;
    for(kx_ind=(lx_/2)-1;kx_ind>=0;kx_ind--){
        kx = (2.0*PI*kx_ind)/(1.0*lx_);
        ky = (2.0*PI*kx_ind)/(1.0*ly_);

        for(int om=0;om<w_size;om++){
            omega = Parameters_BHZ_.w_min + om*dw;
            Akw_ = Zero_Complex;

            for(int r1=0;r1<cells_;r1++){
                r1x = r1/ly_;       r1y = r1%ly_;
                for(int r2=0;r2<cells_;r2++){
                    r2x = r2/ly_;       r2y = r2%ly_;

                    Akw_ += (1.0/(4.0*cells_))*( exp(Iota_Complex*kx*(1.0*(r1x-r2x)))*exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*
                            B_mat[r1][r2][om] );
                }
            }
            file_Akw_out<<k_ind<<"      "<<omega<<"     "<<Akw_.real()<<endl;
        }
        file_Akw_out<<endl;
        k_ind = k_ind + 1;
    }

}


void Observables_BHZ::Calculate_Spectral_Function(){

    string file_Akxw("Akxw.txt");
    ofstream file_Akxw_out(file_Akxw.c_str());

    complex<double> Akxw_;
    omega=0.0;

    for(kx_ind=-lx_/2;kx_ind<=lx_/2;kx_ind++){
        kx=((2*PI*kx_ind)/(1.0*lx_));

        for(int om=0;om<w_size;om++){
            omega=Parameters_BHZ_.w_min + om*dw;
            Akxw_ = Zero_Complex;

            for(r1x=0;r1x<lx_;r1x++){
                for(r2x=0;r2x<lx_;r2x++){

                    Akxw_ += (1.0/(4.0*lx_))*( exp(Iota_Complex*kx*(1.0*(r1x-r2x)))*B_mat[r1x][r2x][om] );

                }
            }
            file_Akxw_out<<kx<<"    "<<omega<<"     "<<Akxw_.real()<<endl;
        }
        file_Akxw_out<<endl;
    }
}


void Observables_BHZ::Calculate_Spin_Currents(){

    Mat_4_Complex_doub SpinCurrent_UP, SpinCurrent_DN;
    
    SpinCurrent_UP.resize(lx_*ly_);   SpinCurrent_DN.resize(lx_*ly_);
    for(int cell1=0;cell1<lx_*ly_;cell1++){
        SpinCurrent_UP[cell1].resize(lx_*ly_);     SpinCurrent_DN[cell1].resize(lx_*ly_);

        for(int cell2=0;cell2<lx_*ly_;cell2++){
            SpinCurrent_UP[cell1][cell2].resize(orbs_);     SpinCurrent_DN[cell1][cell2].resize(orbs_);

            for(int orb1=0;orb1<orbs_;orb1++){
                SpinCurrent_UP[cell1][cell2][orb1].resize(orbs_);     SpinCurrent_DN[cell1][cell2][orb1].resize(orbs_);

                for(int orb2=0;orb2<orbs_;orb2++){
                    SpinCurrent_UP[cell1][cell2][orb1][orb2]=Zero_Complex;
                    SpinCurrent_DN[cell1][cell2][orb1][orb2]=Zero_Complex;
                }
            }
        }
    }

    int m1, m2;
    int cell_ind1, cell_ind2;

    for(int s=0;s<2;s++){

        for(int r1x=0;r1x<lx_;r1x++){
            for(int r1y=0;r1y<ly_;r1y++){
                cell_ind1 = r1y + ly_*r1x;

                for(int orb1=0;orb1<orbs_;orb1++){
                    m1 = s*mhs_ + orb1 + 2*r1y + 2*ly_*r1x;

                    for(int r2x=0;r2x<lx_;r2x++){
                        for(int r2y=0;r2y<ly_;r2y++){
                            cell_ind2 = r2y + ly_*r2x;

                            for(int orb2=0;orb2<orbs_;orb2++){
                                m2 = s*mhs_ + orb2 + 2*r2y + 2*ly_*r2x;

                                if(s==0){
                                    if(Connection_BHZ_.C_mat(m1,m2).real()!=0 || Connection_BHZ_.C_mat(m1,m2).imag()!=0){
                                        for(int n=0;n<size_;n++){
                                            SpinCurrent_UP[cell_ind1][cell_ind2][orb1][orb2] += One_Complex*Connection_BHZ_.C_mat(m1,m2)*
                                                        ( (conj(evecs_(n,m1)) )*((evecs_(n,m2))) )*( Hamiltonian_BHZ_.FermiFunction(evals_(n),mu_val) );
                                        }
                                    }
                                }

                                if(s==1){
                                    if(Connection_BHZ_.C_mat(m1,m2).real()!=0 || Connection_BHZ_.C_mat(m1,m2).imag()!=0){
                                        for(int n=0;n<size_;n++){
                                            SpinCurrent_DN[cell_ind1][cell_ind2][orb1][orb2] += One_Complex*Connection_BHZ_.C_mat(m1,m2)*
                                                        ( (conj(evecs_(n,m1)) )*((evecs_(n,m2))) )*( Hamiltonian_BHZ_.FermiFunction(evals_(n),mu_val) );
                                        }
                                    }
                                }
                            }

                        }
                    }

                }
            }
        }

    }

    Mat_2_Complex_doub SC_UP_TOT, SC_DN_TOT;
    SC_UP_TOT.resize(lx_*ly_);      SC_DN_TOT.resize(lx_*ly_);
    for(int cell1=0;cell1<lx_*ly_;cell1++){
        SC_UP_TOT[cell1].resize(lx_*ly_);
        SC_DN_TOT[cell1].resize(lx_*ly_);

        for(int cell2=0;cell2<lx_*ly_;cell2++){
            SC_UP_TOT[cell1][cell2]=Zero_Complex;
            SC_DN_TOT[cell1][cell2]=Zero_Complex;
        }
    }

    int m1_up,m1_dn,m2_up,m2_dn;

    string file_spin_up_current="Total_spin_up_current_at_each_cell_link.txt";
    ofstream file_spin_up_current_out(file_spin_up_current.c_str());

    string file_spin_dn_current="Total_spin_dn_current_at_each_cell_link.txt";
    ofstream file_spin_dn_current_out(file_spin_dn_current.c_str());

    for(int r1x=0;r1x<lx_;r1x++){
        for(int r1y=0;r1y<ly_;r1y++){
            cell_ind1 = r1y + ly_*r1x;

            for(int r2x=0;r2x<lx_;r2x++){
                for(int r2y=0;r2y<ly_;r2y++){
                    cell_ind2 = r2y + ly_*r2x;

                    for(int orb1=0;orb1<orbs_;orb1++){
                        for(int orb2=0;orb2<orbs_;orb2++){
                            
                            SC_UP_TOT[cell_ind1][cell_ind2] += SpinCurrent_UP[cell_ind1][cell_ind2][orb1][orb2];
                            SC_DN_TOT[cell_ind1][cell_ind2] += SpinCurrent_DN[cell_ind1][cell_ind2][orb1][orb2];

                        }
                    }

                    file_spin_up_current_out<<cell_ind1<<"  "<<cell_ind2<<"     "<<(-1.0)*SC_UP_TOT[cell_ind1][cell_ind2].imag()<<endl;
                    file_spin_dn_current_out<<cell_ind1<<"  "<<cell_ind2<<"     "<<(1.0)*SC_DN_TOT[cell_ind1][cell_ind2].imag()<<endl;

                }
            }
        }
    }

    string file_avg_spin_current="Avg_spin_current_along_ry_vs_rx.txt";
    ofstream file_avg_spin_current_out(file_avg_spin_current.c_str());

    complex<double> avg_SC_val_up,avg_SC_val_dn;
    for(int r1x=0;r1x<lx_;r1x++){
        avg_SC_val_up=Zero_Complex;
        avg_SC_val_dn=Zero_Complex;

        for(int r1y=0;r1y<ly_;r1y++){
            cell_ind1 = r1y + ly_*r1x;

            for(int r2y=0;r2y<ly_;r2y++){
                cell_ind2 = r2y + ly_*r1x;

                avg_SC_val_up += SC_UP_TOT[cell_ind1][cell_ind2];
                avg_SC_val_dn += SC_DN_TOT[cell_ind1][cell_ind2];
            }
        }
        file_avg_spin_current_out<<r1x<<"   "<<(-1.0)*avg_SC_val_up.imag()<<"     "<<avg_SC_val_dn.imag()<<endl;
    }



}


void Observables_BHZ::Calculate_Spin_Chern_Number(){

}