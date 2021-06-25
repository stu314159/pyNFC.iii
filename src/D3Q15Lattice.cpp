#include "Lattice.h"
#include "D3Q15Lattice.h"


D3Q15Lattice::D3Q15Lattice(const int Nx, const int Ny, const int Nz):
Lattice(Nx,Ny,Nz),
ex{0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1},
ey{0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1},
ez{0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1},
w{2.f/9.f,1.f/9.f,1.f/9,1.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,
    1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f,
    1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f},
bbSpd{0,2,1,4,3,6,5,14,13,12,11,10,9,8,7},
Qflat{1./3.,0.,0.,0,-1./3.,0,0,0,-1./3.,  //this is embarrassing, but better than implementing the LA
	2./3.,0,0,0,-1./3.,0,0,0,-1./3.,
	2./3.,0,0,0,-1./3.,0,0,0,-1./3.,
	-1./3.,0,0,0,2./3.,0,0,0,-1./3.,
	-1./3.,0,0,0,2./3.,0,0,0,-1./3.,
	-1./3.,0,0,0,-1./3.,0,0,0,2./3.,
	-1./3.,0,0,0,-1./3.,0,0,0,2./3.,
	2./3.,1,1,1,2./3.,1,1,1,2./3.,
	2./3.,-1,-1,-1,2./3.,1,-1,1,2./3.,
	2./3.,-1,1,-1,2./3.,-1,1,-1,2./3.,
	2./3.,1,-1,1,2./3.,-1,-1,-1,2./3.,
	2./3.,1,-1,1,2./3.,-1,-1,-1,2./3.,
	2./3.,-1,1,-1,2./3.,-1,1,-1,2./3.,
	2./3.,-1,-1,-1,2./3.,1,-1,1,2./3.,
	2./3.,1,1,1,2./3.,1,1,1,2./3.}
{
	// direct base-class pointers to lattice variables
	setNumSpd(numSpd);
	setEx(ex);
	setEy(ey);
	setEz(ez);
	setW(w);
	setBBspd(bbSpd);
	setQflat(Qflat);

}

D3Q15Lattice::~D3Q15Lattice()
{


}

void D3Q15Lattice::set_inlet_bc_micro(LBM_DataHandler& f)
{
  int sp[5]={5,7,8,9,10};
  int bbSp[5]={6,11,12,13,14};
  int numBB = 5;
  for(int s=0;s<numBB;s++)
  {
	  f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
  }
}

void D3Q15Lattice::set_inlet_bc_macro(LBM_DataHandler& f)
{
	f.uz = f.u_bc;
	f.ux = 0.; f.uy = 0.;
	f.rho = (1./(1. - f.uz))*(2.*
			(f.f[6]+f.f[11]+f.f[12]+f.f[13]+f.f[14])+
			(f.f[0]+f.f[1]+f.f[2]+f.f[3]+f.f[4]));
}

void D3Q15Lattice::set_outlet_bc_micro(LBM_DataHandler& f)
{
	int sp[5]={6,11,12,13,14};
	int bbSp[5]={5,7,8,9,10};
	int numBB = 5;
	for(int s=0;s<numBB;s++)
	{
		f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
	}
}

void D3Q15Lattice::set_outlet_bc_macro(LBM_DataHandler& f)
{
	f.rho = f.rho_bc;
	f.uz = -1. + (1./f.rho)*
			(2.*(f.f[5]+f.f[7]+f.f[8]+f.f[9]+f.f[10])+
					(f.f[0]+f.f[1]+f.f[2]+f.f[3]+f.f[4]));
}

