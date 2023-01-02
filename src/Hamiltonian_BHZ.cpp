#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <complex>

extern "C" {
    #include <lapacke.h>
    #include <cblas.h>
}

#include "Tensor.hpp"
#include "Parameters_BHZ.hpp"
#include "Connection_BHZ.hpp"
#include "Hamiltonian_BHZ.hpp"

using namespace std;

void Hamiltonian_BHZ::Initialize(){
    Evals_.resize(size_);
    Evecs_.resize(size_);
    Ham_.resize(size_);
    for(int i=0;i<size_;i++){
        Evecs_[i].resize(size_);
        Ham_[i].resize(size_);
    }

}



void Hamiltonian_BHZ::Diagonalizer(Mat_2_Complex_doub Ham_){
    int LDA=size_;
    int info;
    /* Local arrays */
    double* eval = new double[size_];
    lapack_complex_double* mat = (lapack_complex_double *) calloc(size_*size_, sizeof(lapack_complex_double));

    for(int i=0;i<size_;i++){
        for(int j=0;j<size_;j++){
            mat[i*(size_)+j] = Ham_[i][j].real()+ Ham_[i][j].imag()*I;
            //mat[i*(size_)+j].real() = Ham_[i][j].real();
            //mat[i*(size_)+j].imag() = Ham_[i][j].imag();
        }
    }

    info=LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'L', size_, mat, LDA, eval);

    if(info > 0){
        cout<< "The LAPACKE_zheev failed to diagonalize."<<endl;
    }

    for(int i=0;i<size_;i++){
        Evals_[i]=eval[i];
        for(int j=0;j<size_;j++){
            Evecs_[j][i].real(lapack_complex_double_real(mat[i*(size_)+j]));
            Evecs_[j][i].imag(lapack_complex_double_imag(mat[i*(size_)+j]));
        }
    }

    string Evals_out="Eigenvalues.txt";
    ofstream Evals_file_out(Evals_out.c_str());

    for(int i=0;i<size_;i++){
        Evals_file_out<<eval[i]<<endl;
    }
}



double Hamiltonian_BHZ::ChemicalPotential(double muin_, double particles_){
    double mu_temp, eps_, dmu_by_dN, N_temp, dmu_by_dN_min, Ne_;
    eps_=1e-2;
    mu_temp=Evals_[0];
    N_temp=100000;
    Ne_=1.0*particles_;
    
    int iters=0;

    dmu_by_dN = 0.01*( Evals_[size_-1] - Evals_[0] )*( 1.0/(1.0*size_) );
    dmu_by_dN_min = 0.0001*( Evals_[size_-1] - Evals_[0] )*( 1.0/(1.0*size_) );

    while( abs(N_temp - Ne_) > eps_){
        N_temp=0;
        for(int i=0;i<size_;i++){
            N_temp += FermiFunction(Evals_[i],mu_temp);
        }

        mu_temp = mu_temp + dmu_by_dN*(Ne_-N_temp);
        iters++;

        if(iters%1000==0){
            dmu_by_dN = (1000.0/(10.0*iters))*dmu_by_dN;
            dmu_by_dN = max(dmu_by_dN, dmu_by_dN_min);
        }
    }

    cout<<"Calculated mu = "<<mu_temp<<endl;
    cout<<"Calculated # of particles = "<<N_temp<<endl;

    return mu_temp;
}



double Hamiltonian_BHZ::FermiFunction(double en_, double mu_){
    double ffn, temp;
    temp=Parameters_BHZ_.Temp;

    ffn = ( (1.0)/(1.0 + exp( (en_ - mu_)/temp )) );

    return ffn;
}