#ifndef OUTFILES_H
#define OUTFILES_H

#include "readfiles.h"

using namespace std;

class OutSatVeri
{
public:
	OutSatVeri()
	{
		_sys = '\0';
		_prn = 0;
		_fixflag = 0;
		_resflag = 0;
		_outflag = 0;
		_ngrid = 0;
		_el = 0.0;
		_dstec[0] = _dstec[1] = _dstec[2] = 0.0;
		_dDstec[0] = _dDstec[1] = 0.0;
	}
	~OutSatVeri(){}

	char _sys;
	int _prn;
	int _fixflag;
	int _resflag;
	int _outflag;
	int _ngrid;
	double _el;
	double _dstec[3];
	double _dDstec[2];
private:
};

class OutRovStec
{
public:
	OutRovStec()
	{
		_time = { 0 };
		_name = "";
		for (int i = 0; i < NUMSYS; i++) {
			_satNum[i] = _satRef[i] = 0;
			_rate[i] = 0.0;
		}
	}
	~OutRovStec() {}

	Gtime _time;
	string _name;
	int _satNum[NUMSYS];
	int _satRef[NUMSYS];
	double _rate[NUMSYS];
	map<int, OutSatVeri> _satveris[NUMSYS];

private:
};

bool checkRef(IN int sys, IN int ref, IN SiteAtmo& roviono);
bool checkSat(IN int sys, IN int ref, IN int prn, IN SiteAtmo& roviono, IN double minEl);
bool findStecModelSat(IN int sys, IN int prn, IN ProStecMod& stecmod, OUT ProStecModSat* satmod);
double rovRes(IN double lat, IN double lon, IN GridInfo& grid, IN ProStecModSat* satmod, IN ProStecModSat* refmod, OUT int* n);
void rovStecDiff(IN Coption& cfg, IN GridInfo& grid, IN FILE* fp, IN SiteAtmo& roviono, IN ProStecMod& stecmod, OUT OutRovStec& rovout);
void outRovStec(IN Coption& cfg, IN GridInfo& grid, IN SiteAtmos& rovaug, IN ProStecMod& stecmod, IN FileFps& rovfps, IN int type);
void createRovFile(IN Coption& cfg, OUT FileFps& fps);


#endif
