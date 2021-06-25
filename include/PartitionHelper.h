/*
 * PartitionHelper.h
 *
 *  Created on: Aug 4, 2017
 *      Author: stu
 */

#ifndef PARTITIONHELPER_H_
#define PARTITIONHELPER_H_

#include "Lattice.h"
#include "D3Q15Lattice.h"
#include "D3Q19Lattice.h"
#include "D3Q27Lattice.h"
#include <boost/python.hpp>
#include <cstdlib>


class PartitionHelper {

public:
	PartitionHelper(const int Nx, const int Ny, const int Nz, const int numSpd);
	~PartitionHelper();
	void setAdjacency(boost::python::dict& adjDict);

private:
	const int numSpd;
	const int Nx,Ny,Nz;
	Lattice * myLattice;
	float * ex;
	float * ey;
	float * ez;
	void set_latticeVelocities();


};



#endif /* PARTITIONHELPER_H_ */
