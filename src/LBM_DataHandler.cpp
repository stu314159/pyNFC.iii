/*
 * LBM_DataHandler.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: stu
 */


#include "LBM_DataHandler.h"
#include <cstdlib>



LBM_DataHandler::LBM_DataHandler(const int numSpd) :
ux(0),uy(0),uz(0),rho(0),u_bc(0),rho_bc(0),nodeType(0),dynamics(2),omega(0), Cs(0),

piFlat{0,0,0,0,0,0,0,0,0},S{0,0,0,0,0,0,0,0,0},numSpd(numSpd), omegaMRT(NULL)
{
	f = new float[numSpd];
	fEq = new float[numSpd];
	fOut = new float[numSpd];


}

LBM_DataHandler::~LBM_DataHandler()
{
	delete [] f;
	delete [] fEq;
	delete [] fOut;

}

int LBM_DataHandler::get_numSpd()
{
	return numSpd;
}


