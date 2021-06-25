/*
 * D3Q19Lattice.h
 *
 *  Created on: Jul 10, 2017
 *      Author: sblair
 */

#ifndef D3Q19LATTICE_H_
#define D3Q19LATTICE_H_

#include "Lattice.h"


class D3Q19Lattice: public Lattice
{
public:
	D3Q19Lattice(const int Nx, const int Ny, const int Nz);
	~D3Q19Lattice();
	void set_inlet_bc_micro(LBM_DataHandler& f);
	void set_inlet_bc_macro(LBM_DataHandler& f);
	void set_outlet_bc_micro(LBM_DataHandler& f);
	void set_outlet_bc_macro(LBM_DataHandler& f);

private:
	static const int numSpd=19;
	float ex[numSpd];
	float ey[numSpd];
	float ez[numSpd];
	float w[numSpd];
	int bbSpd[numSpd];
	float Qflat[numSpd*9];

};


#endif /* D3Q19LATTICE_H_ */
