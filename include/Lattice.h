#ifndef LATTICE_H
#define LATTICE_H

#include "LBM_DataHandler.h"
#include <omp.h>

class Lattice{
  public:
    Lattice(const int Nx, const int Ny, const int Nz);
    virtual ~Lattice();
    // "setter/getter functions for lattice variables"
    int getNx(){return Nx;}
    int getNy(){return Ny;}
    int getNz(){return Nz;}
    void setNumSpd(const int ns){numSpd=ns;}
    int getNumSpd(){return numSpd;}
    void setEx(float * x){ex = x;}
    void setEy(float * y){ey = y;}
    void setEz(float * z){ez = z;}
    void setW(float * tw){w = tw;}
    void setBBspd(int* BB){bbSpd = BB;}
    void setQflat(float * Qf){Qflat = Qf;}
   // void loadF(const int idx, const float val){f[idx]=val;} //<-- ?? consider changing this in light of LDH
    float * get_ex(){return ex;}
    float * get_ey(){return ey;}
    float * get_ez(){return ez;}
    float * get_w(){return w;}
    int * get_bb(){return bbSpd;}
    // common functions for LBM stream/collide tasks
    void computeMacroscopicData(float& rho,float& ux, float& uy, float& uz, const float * f);
    void computeMacroscopicData(LBM_DataHandler& f);
    void computeEquilibrium(float * fEq, const float ux, const float uy, const float uz, const float rho);
    void computeFout(LBM_DataHandler& f);
    void computeEquilibrium(LBM_DataHandler& f);
    void bounceBack(LBM_DataHandler& f);
    void compute_piFlat(LBM_DataHandler& f);
    void regularize(LBM_DataHandler& f);
    void relax(LBM_DataHandler& f);
    void relaxMRT(LBM_DataHandler& f);
    void set_Vz_micro(LBM_DataHandler& f);
    void computeStrainTensor(LBM_DataHandler & f);
    void applyTurbulenceModel(LBM_DataHandler &f);

    virtual void set_inlet_bc_micro(LBM_DataHandler& f) = 0;
    virtual void set_inlet_bc_macro(LBM_DataHandler& f) = 0;
    virtual void set_outlet_bc_micro(LBM_DataHandler& f) = 0;
    virtual void set_outlet_bc_macro(LBM_DataHandler& f) = 0;

    // declarations for lattice variables common to all subclasses


  private:
    const int Nx, Ny, Nz; // global lattice size
    // data members that all classes deriving from Lattice must have (and initialize)
    int numSpd;
//    float ux,uy,uz,rho; //macroscopic data variables
//    float rho_bc,u_bc; // macroscopic boundary conditions
    float * ex;
    float * ey;
    float * ez;
    float * w;
    int * bbSpd;
    float * Qflat;

};
#endif
