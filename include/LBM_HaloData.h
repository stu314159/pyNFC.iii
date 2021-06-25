#ifndef LBM_HALODATA_H_
#define LBM_HALODATA_H_

class LBM_HaloData
{
public:

	LBM_HaloData(int numData, int numSpd);
	LBM_HaloData();
	~LBM_HaloData();
	LBM_HaloData(const LBM_HaloData& hd); // copy constructor
	LBM_HaloData& operator=(const LBM_HaloData& hd)=default; //compiler generate shallow copy
	void set_numData(const int nd){numData = nd;};
	void set_numSpd(const int s){numSpd=s;};
	void set_local_nn(int* nn){local_nn = nn;};
	void set_spd(int * s){spd = s;};
	void set_data(float * db){data_buf = db;};
	int get_numData(){return numData;};
	int get_numSpd(){return numSpd;};
	void extractHaloData(const float * f);
	void distributeHaloData(float * f);

private:

	// these will hold the pointers to the python HDO object
	int * local_nn; // pointer to array of local node numbers
	int * spd; // pointer to array of speeds
	float * data_buf; // pointer to array of data

	int numData; // number of data items.
	int numSpd; // number of speeds for the associated LBM lattice (needed for dist)
};


#endif
