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
		_nsta = 0;
		_fixflag = 0;
		_resflag = 0;
		_outflag = 0;
		_ngrid = 0;
		_el = 0.0;
		_elm = 0.0;
		_dstec[0] = _dstec[1] = _dstec[2] = 0.0;
		_dDstec[0] = _dDstec[1] = 0.0;
		_QI[0] = _QI[1] = 0.0;
		_pergrid[0] = _pergrid[1] = _pergrid[2] = _pergrid[3] = 0;
	}
	~OutSatVeri(){}

	char _sys;
	int _prn;
	int _nsta;
	int _fixflag;
	int _resflag;
	int _outflag;
	int _ngrid;
	int _pergrid[4];
	double _el;
	double _elm;
	double _dstec[3];
	double _dDstec[2];
	double _QI[2];
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
			//_rate[i] = 0.0;
		}
	}
	~OutRovStec() {}

	Gtime _time;
	string _name;
	int _satNum[NUMSYS];
	int _satRef[NUMSYS];
	//double _rate[NUMSYS];
	map<int, OutSatVeri> _satveris[NUMSYS];

private:
};

bool checkRef(IN int sys, IN int ref, IN SiteAtmo& roviono);
bool checkSat(IN int sys, IN int ref, IN int prn, IN SiteAtmo& roviono, IN double minEl);
bool findStecModelSat(IN int sys, IN int prn, IN ProStecMod& stecmod, OUT ProStecModSat* satmod);
double rovRes(IN double lat, IN double lon, IN GridInfo& grid, IN ProStecModSat* satmod, IN ProStecModSat* refmod, OUT int* n, OUT int* stas);
double rovResMSF(IN double lat, IN double lon, IN GridInfo& grid, IN ProStecModSat* satmod, IN ProStecModSat* refmod, OUT int* n, OUT int* stas);
void rovStecDiff(IN Coption& cfg, IN GridInfo& grid, IN FILE* fp, IN SiteAtmo& roviono, IN ProStecMod& stecmod, OUT OutRovStec& rovout);
void outRovStec(IN Coption& cfg, IN GridInfo& grid, IN SiteAtmos& rovaug, IN ProStecMod& stecmod, IN FileFps& rovfps, IN int type);
void writeModelHead(IN  Coption& cfg, IN FILE* fp);
void printSatMod(IN int useres, IN ProStecModSat& dat, IN FILE* fp);
void printEpochSTEC(IN Gtime tnow, IN ProStecMod& stecmod, IN int ngrid, IN FILE* fp);
void outStecModel(IN Gtime tnow, IN  Coption& cfg, IN GridInfo& grid, IN ProStecMod& stecmod, IN FILE* fp);
void printEpochZTD(IN Gtime tnow, IN ProZtdMod& ztdmod, IN FILE* fp);
void printEpoch(IN Gtime tnow, IN bool* stat, IN ProStecMod& stecmod, IN ProZtdMod& ztdmod, IN int ngrid, IN FILE* fp);
void outZtdModel(IN Gtime tnow, IN GridInfo& grid, IN ProZtdMod& ztdmod, IN FILE* fp);
void outAtmoModel(IN Gtime tnow, IN  Coption& cfg, IN GridInfo& grid, IN bool* stat,
	IN ProStecMod& stecmod, IN ProZtdMod& ztdmod, IN FILE* fp);
void createRovFile(IN Coption& cfg, OUT FileFps& fps);

#endif
