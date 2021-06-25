/*
 * PartitionHelper.cpp
 *
 *  Created on: Aug 4, 2017
 *      Author: stu
 */

#include "PartitionHelper.h"
#include <iostream>

PartitionHelper::PartitionHelper(const int nx, const int ny, const int nz, const int ns):
Nx(nx), Ny(ny), Nz(nz), numSpd(ns)
{
	// create the lattice object
		switch (numSpd)
		{
		case(15):
				myLattice = new D3Q15Lattice(Nx,Ny,Nz); break;
		case(19):
				myLattice = new D3Q19Lattice(Nx,Ny,Nz); break;
		case(27):
				myLattice = new D3Q27Lattice(Nx,Ny,Nz);

		}
		set_latticeVelocities();
}

PartitionHelper::~PartitionHelper()
{

}




void PartitionHelper::setAdjacency(boost::python::dict& adjDict)
{
	int dx,dy,dz;
	int tX,tY,tZ;
	int gId,tId;

	for(int z=0;z<Nz;z++)
	{
		for(int y=0;y<Ny; y++)
		{
			for(int x=0;x<Nx; x++)
			{
				gId = x + y*Nx + z*Nx*Ny;
				//std::vector<int> val(numSpd);
				boost::python::list l;
				for(int spd=0;spd<numSpd;spd++)
				{
					dx = (int)ex[spd]; dy = (int)ey[spd]; dz = (int)ez[spd];

					tX = (x+dx)%Nx; tY=(y+dy)%Ny; tZ = (z+dz)%Nz;
					if (tX<0)
					{
						tX = Nx-1;
					}
					if (tY<0)
					{
						tY = Ny-1;
					}
					if (tZ<0)
					{
						tZ = Nz-1;
					}
					tId = tX+tY*Nx+tZ*Nx*Ny;
					//val[spd] = tId;
					l.append(tId);

				}
				adjDict[gId] = l;

			}
		}
	}


}

void PartitionHelper::set_latticeVelocities()
{
	ex = myLattice->get_ex();
	ey = myLattice->get_ey();
	ez = myLattice->get_ez();
}

using namespace boost::python;

BOOST_PYTHON_MODULE(PartitionHelper)
{
	class_<PartitionHelper>("PartitionHelper",init<const int,const int,const int,const int>())
			.def("setAdjacency",&PartitionHelper::setAdjacency)
			;
}



