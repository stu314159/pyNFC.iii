#ifndef D3Q15LATTICE_H
#define D3Q15LATTICE_H

#include "Lattice.h"

class D3Q15Lattice: public Lattice
{
  public:
    D3Q15Lattice(const int Nx, const int Ny, const int Nz);
    ~D3Q15Lattice();
    void set_inlet_bc_micro(LBM_DataHandler& f);
    void set_inlet_bc_macro(LBM_DataHandler& f);
    void set_outlet_bc_micro(LBM_DataHandler& f);
    void set_outlet_bc_macro(LBM_DataHandler& f);


  private:
    static const int numSpd=15;
    float ex[numSpd];
    float ey[numSpd];
    float ez[numSpd];
    float w[numSpd];
    int bbSpd[numSpd];
    float Qflat[numSpd*9];


};

#endif
