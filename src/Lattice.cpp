#include "Lattice.h"
#include <cstdlib>
#include <iostream> //<-- used for debugging
#include <cmath> //<-- used for turbulence model sqrt calc

Lattice::Lattice(const int Nx, const int Ny, const int Nz):
Nx(Nx), Ny(Ny), Nz(Nz), numSpd(0),

ex(NULL), ey(NULL), ez(NULL),w(NULL), bbSpd(NULL), Qflat(NULL)
{

}

Lattice::~Lattice()
{


}

void Lattice::bounceBack(LBM_DataHandler& f)
{
	//float fTmp[numSpd];
	for(int spd=0;spd<numSpd;spd++)
	{
		//fTmp[spd]=f.f[bbSpd[spd]];
		f.fOut[spd] = f.f[bbSpd[spd]];
	}
//	for(int spd=0;spd<numSpd;spd++)
//	{
//		f.f[spd]=fTmp[spd];
//	}

}
void Lattice::computeMacroscopicData(float& rho, float& ux, float& uy, float& uz, const float* f)
{
	rho = 0.; ux = 0.; uy = 0.; uz = 0.;
	for(int spd = 0; spd<numSpd; spd++)
	{
		rho+=f[spd];
		ux+=ex[spd]*f[spd];
		uy+=ey[spd]*f[spd];
		uz+=ez[spd]*f[spd];
	}
	ux/=rho; uy/=rho; uz/=rho;

}
void Lattice::computeMacroscopicData(LBM_DataHandler& f)
{
	float rho, ux, uy, uz;
	rho = 0.; ux = 0.; uy = 0.; uz = 0.;
	for(int spd = 0; spd<numSpd; spd++)
	{

		rho+=f.f[spd];
		ux+=ex[spd]*f.f[spd];
		uy+=ey[spd]*f.f[spd];
		uz+=ez[spd]*f.f[spd];
	}
	ux/=rho; uy/=rho; uz/=rho;
	f.rho = rho; f.ux = ux; f.uy = uy; f.uz = uz;
	switch (f.nodeType)
	{
	case 1:
		f.ux = 0.; f.uy = 0.; f.uz = 0.;
		break;
	case 2:
		f.uz = f.u_bc;
	}

}

void Lattice::computeEquilibrium(float * fEq,const float ux, const float uy, const float uz, const float rho)
{
	//this will assume that up-to-date data exists in rho, ux, uy, and uz class variables
	//and that ex,ey,ez,and w have been updated with the lattice-specific info
	float cu;

	for(int spd = 0; spd<numSpd;spd++)
	{
		fEq[spd] = 0;
		cu = 3.f*(ex[spd]*ux+ey[spd]*uy+ez[spd]*uz);
		fEq[spd]=w[spd]*rho*(1.f+cu+0.5f*(cu*cu)-3.f/2.f*(ux*ux+uy*uy+uz*uz));

	}
}

void Lattice::computeEquilibrium(LBM_DataHandler& f)
{

	computeEquilibrium(f.fEq,f.ux,f.uy,f.uz,f.rho);
}

void Lattice::compute_piFlat(LBM_DataHandler& f)
{

	//f.piFlat = {0,0,0,0,0,0,0,0,0}; // initialized in constructor
	for(int k=0; k<9;k++)
	{
		f.piFlat[k]=0.;
	}
	float fNeq;
	for(int spd = 0; spd<numSpd; spd++)
	{
		fNeq = f.f[spd] - f.fEq[spd];
		f.piFlat[0]+=ex[spd]*ex[spd]*fNeq;
		f.piFlat[1]+=ey[spd]*ex[spd]*fNeq;
		f.piFlat[2]+=ez[spd]*ex[spd]*fNeq;
		f.piFlat[3]+=ex[spd]*ey[spd]*fNeq;
		f.piFlat[4]+=ey[spd]*ey[spd]*fNeq;
		f.piFlat[5]+=ez[spd]*ey[spd]*fNeq;
		f.piFlat[6]+=ex[spd]*ez[spd]*fNeq;
		f.piFlat[7]+=ey[spd]*ez[spd]*fNeq;
		f.piFlat[8]+=ez[spd]*ez[spd]*fNeq;

	}
}

void Lattice::regularize(LBM_DataHandler& f)
{

	float wa;

	for(int spd = 0; spd<numSpd; spd++)
	{
		// get leading constant
		wa = (9.f/2.f)*w[spd];
		f.f[spd]=f.fEq[spd];
		// load chunk of Qflat:
		for(int k=0;k<9;k++)
		{
			f.f[spd]+= wa*Qflat[spd*9+k]*f.piFlat[k];
		}
	}
}

void Lattice::relax(LBM_DataHandler& f)
{

	for(int spd=0;spd<numSpd;spd++)
	{
		f.fOut[spd]=f.f[spd]-f.omega*(f.f[spd] - f.fEq[spd]);
	}

}

void Lattice::relaxMRT(LBM_DataHandler& f)
{
	for(int spd=0;spd<numSpd;spd++)
	{
		f.fOut[spd] = 0.;
		for(int j=0;j<numSpd;j++)
		{
			f.fOut[spd]+=f.omegaMRT[spd*numSpd+j]*(f.f[j]-f.fEq[j]);
		}
		f.fOut[spd] = f.f[spd] + f.fOut[spd];//<-- remember the minus sine in omegaMRT
	}
}

void Lattice::set_Vz_micro(LBM_DataHandler & f)
{
	for(int spd = 0; spd<numSpd; spd++)
	{
		// push microscopic BCs towards desired value
		f.f[spd]+=3.*f.rho*w[spd]*(ez[spd]*(f.u_bc - f.uz) +
				ex[spd]*(0. - f.ux) + ey[spd]*(0.-f.uy));
	}
	// update macro BCs so Equilibrium will be calculated as desired.
	f.ux = 0.; f.uy = 0.; f.uz = f.u_bc;
}

void Lattice::computeStrainTensor(LBM_DataHandler & f)
{
  // here we take advantage of the fact that f.S is initialized to zero
	float e[3];
	const int nDim = 3;
	for(int spd = 0; spd<numSpd; spd++)
	{
		e[0] = ex[spd]; e[1] = ey[spd]; e[2] = ez[spd];
		for(int i = 0; i<nDim; i++)
		{
			for(int j=0; j<nDim; j++)
			{
				f.S[i*nDim+j]+=e[i]*e[j]*(f.f[spd] - f.fEq[spd]);
			}
		}
	}
}

void Lattice::applyTurbulenceModel(LBM_DataHandler & f)
{
	float nu, nu_e;
	nu = ((1./f.omega) - 0.5)/3.0;
	float P;
	P = sqrt(f.S[0]*f.S[0]+f.S[4]*f.S[4]+f.S[8]*f.S[8] +
			2.*(f.S[1]*f.S[1]+f.S[2]*f.S[2] + f.S[5]*f.S[5]));
	P*=f.Cs; P = sqrt(P+nu*nu) - nu;
	nu_e = P/6.;
	f.omega = 1./(3.*(nu+nu_e)+0.5);

}

void Lattice::computeFout(LBM_DataHandler& f)
{
	// compute macroscopic velocity and pressure
	computeMacroscopicData(f);
	if(f.nodeType==1) // solid node
	{
		f.ux = 0.; f.uy = 0.; f.uz = 0; // solid nodes, zero velocity
		bounceBack(f);
		return;
	}

	// node type 2 and 3 apply macroscopic boundary conditions
	switch (f.nodeType)
	{
	case 2:
		set_inlet_bc_macro(f);
		break;
	case 3:
		set_outlet_bc_macro(f);
		break;
	case 5:
		set_Vz_micro(f);
		break;
	}


	// compute equilibrium
	computeEquilibrium(f);

	switch(f.nodeType)
	{
	case 2:
		set_inlet_bc_micro(f);
		break;
	case 3:
		set_outlet_bc_micro(f);
		break;
	}


	// if Cs>0, apply turbulence model here to adjust f.omega
	if(f.Cs > 0)
	{
		computeStrainTensor(f);
		applyTurbulenceModel(f);
	}

	// get (flattened) second-order moment of particle density distribution

	switch(f.dynamics)
	{
	case 1:
		relax(f); break;

	case 2:
		compute_piFlat(f);
		regularize(f);
		relax(f);
		break;
	case 3:
		relaxMRT(f);

	}


}
