#ifndef Parameters_BHZ_HPP
#define Parameters_BHZ_HPP

#include <string>

//This class will read from the **input file**

class Parameters_BHZ{

        public:
                int Lx, Ly, N_Orbs, Wx;
                int Total_Sites, Total_Cells, Ham_Size;
                int Total_Particles;
                double A_val, B_val, M_val, Vo_val;
                double Fill;
                bool PBC_X,PBC_Y;

                double Temp,beta;

                double w_min,w_max;
                double dw_dos,eta_dos;

                void Initialize(string input_file);
                double matchstring(string file, string match);
                string matchstring2(string file, string match);
};

#endif