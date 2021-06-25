#include "LBM_HaloData.h"
#include <cstdlib>

LBM_HaloData::LBM_HaloData(int numData, int numSpd) :
local_nn(NULL), spd(NULL), data_buf(NULL), numData(numData),numSpd(numSpd)
{


}

LBM_HaloData::LBM_HaloData() :
			local_nn(NULL), spd(NULL), data_buf(NULL), numData(0), numSpd(0) // default constructor
{


}

LBM_HaloData::~LBM_HaloData()
{


}

LBM_HaloData::LBM_HaloData(const LBM_HaloData& hd)
{
	// shallow copy is approproate.  pyNFC partition object holds the deep copy
	numData = hd.numData;
	numSpd = hd.numSpd;
	local_nn = hd.local_nn;
	spd = hd.spd;
	data_buf = hd.data_buf;
}

void LBM_HaloData::extractHaloData(const float * f)
{
	for(int hIdx=0;hIdx<numData;hIdx++)
	{
		data_buf[hIdx] = f[local_nn[hIdx]*numSpd+spd[hIdx]];
	}


}

void LBM_HaloData::distributeHaloData(float * f)
{
	for(int hIdx=0;hIdx<numData;hIdx++)
	{
		f[local_nn[hIdx]*numSpd+spd[hIdx]] = data_buf[hIdx];
	}

}
