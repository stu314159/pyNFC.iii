/*
 * D3Q19Lattice.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: sblair
 */

#include "D3Q19Lattice.h"
#include <cstdlib>

D3Q19Lattice::D3Q19Lattice(const int Nx, const int Ny, const int Nz):
Lattice(Nx,Ny,Nz),
ex{0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0},
ey{0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1},
ez{0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
w{3.f/9.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,
    1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,
    1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f},
bbSpd{0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15},
Qflat{-1./3.,0,0,0,-1./3.,0,0,0,-1./3.,
	2./3.,0,0,0,-1./3.,0,0,0,-1./3.,
	2./3.,0,0,0,-1./3.,0,0,0,-1./3.,
	-1./3.0,0,0,-2./3.,0,0,0,-1./3.,
	-1./3.,0,0,0,2./3.,0,0,0,-1./3.,
	-1./3.,0,0,0,-1./3.,0,0,0,2./3.,
	-1./3.,0,0,0,-1./3.,0,0,0,2./3.,
	2./3.,1,0,1,2./3.,0,0,0,-1./3.,
	2./3.,-1,0,-1,2./3.,0,0,0,-1./3.,
	2./3.,-1,0,-1,2./3.,0,0,0,-1./3.,
	2./3.,1,0,1,2./3.,0,0,0,-1./3.,
	2./3.,0,1,0,-1./3.,0,1,0,2./3.,
	2./3.,0,-1,0,-1./3.,0,-1,0,2./3.,
	2./3.,0,-1,0,-1./3.,0,-1,0,2./3.,
	2./3.,0,1,0,-1./3.,0,1,0,2./3.,
	-1./3.,0,0,0,2./3.,1,0,1,2./3.,
	-1./3.,0,0,0,2./3.,-1,0,-1,2./3.,
	-1./3.,0,0,0,2./3.,-1,0,-1,2./3.,
	-1./3.,0,0,0,2./3.,1,0,1,2./3.}
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

D3Q19Lattice::~D3Q19Lattice()
{

}

void D3Q19Lattice::set_inlet_bc_micro(LBM_DataHandler& f)
{
	int sp[5]={5,11,12,15,16};
	int bbSp[5]={6,14,13,18,17};
	int numBB = 5;
	for(int s=0;s<numBB;s++)
	{
		f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
	}
}

void D3Q19Lattice::set_inlet_bc_macro(LBM_DataHandler& f)
{
	f.uz = f.u_bc;
	f.ux = 0.; f.uy = 0.;
	f.rho = (1./(1.-f.uz))*(2.*(f.f[6]+f.f[13]+f.f[14]+
			f.f[17]+f.f[18])+
			(f.f[0]+f.f[1]+f.f[2]+f.f[3]+f.f[4]+
					f.f[7]+f.f[8]+f.f[9]+f.f[10]));


}

void D3Q19Lattice::set_outlet_bc_micro(LBM_DataHandler& f)
{
	int sp[5]={6,14,13,18,17};
	int bbSp[5]={5,11,12,15,16};
	int numBB = 5;
	for(int s=0;s<numBB;s++)
	{
		f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
	}
}

void D3Q19Lattice::set_outlet_bc_macro(LBM_DataHandler& f)
{
	// check this.
	f.rho = f.rho_bc;
	f.uz = -1. + (1./f.rho)*(2.*(f.f[5]+f.f[11]+f.f[12]+
			f.f[15]+f.f[16])+(f.f[0]+f.f[1]+f.f[2]+
					f.f[3]+f.f[4]+f.f[7]+f.f[8]+f.f[9]+f.f[10]));


}


