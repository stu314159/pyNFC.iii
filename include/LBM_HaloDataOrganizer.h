/*
 * LBM_HaloDataOrganizer.h
 *
 *  Created on: Jul 21, 2017
 *      Author: stu
 */

#ifndef LBM_HALODATAORGANIZER_H_
#define LBM_HALODATAORGANIZER_H_

#include "LBM_HaloData.h"
#include <map>

class LBM_HaloDataOrganizer
{
public:
	LBM_HaloDataOrganizer();
	~LBM_HaloDataOrganizer();
	void insert_ngb(const int ngbNum,const int numData, const int numSpd);
	void initialize_ngb_pointers(const int ngbNum, int * nd_num, int * spd, float * data);
	void extractHaloData(const float * f);
	void distributeHaloData(float * f);


private:
	std::map<int,LBM_HaloData> HaloData;

};




#endif /* LBM_HALODATAORGANIZER_H_ */
