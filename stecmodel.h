#ifndef STEC_MODEL_H
#define STEC_MODEL_H

#include "comminterface.h"
#include "localmodclass.h"
#include "algmatrix.h"
#include "stecroti.h"

/* pos */
class Pos
{
public:
	Pos()
	{
		_lat = _lon = 0.0;
	};
	~Pos() {};
	Pos(double lat, double lon) :_lon(lon), _lat(lat) {};

	bool operator == (const Pos& src) const {
		return (fabs(src._lat - _lat) < DBL_EPSILON) && (fabs(src._lon - _lon) < DBL_EPSILON);
	}

	double _lat;	// latitude  (rad)
	double _lon;	// longitude (rad)
	
private:
};
/* stec grid */
class StecGrid
{
public:
	StecGrid()
	{
		_gridid = 0;
		_lon = _lat = 0.0;
		_stec = _rms = 0.0;		
	};
	~StecGrid() {};
	bool isVaild(IN double lat, IN double lon, IN double dlat, IN double dlon) const;

	int _gridid;	// grid ID
	double _lon;	// longitude (rad)
	double _lat;	// latitude  (rad)
	double _stec;	// stec (TECu)
	double _rms;	// rms (TECu)
};
/* sat stec model info */
class StecModSat
{
public:
	StecModSat()
	{
		_system = '\0';
		_sat = 0;
		for (int i = 0; i < 4; i++) {
			_coeff[i] = _coeff_rms[i] = 0.0;
		}
		for (int i = 0; i < 2; i++) {
			_QI[i] = 0.0;
		}
		_gridNum = 0;
		_satreslevel = 0;
		_ele = 0.0;
		_nsta = 0;
	}
	~StecModSat() {};
	void reset();

	/* sat info */
	char   _system;						// system captial GREC
	int    _sat;						// sat number
	double _ele;						// elevation (rad)
	int    _nsta;
	/* coffeicients & residual */
	double _coeff[4];					// coeffeicient C00 C01 C10 C11
	double _coeff_rms[4];				// coeffeicient rms
	double _QI[2];						// [0]without res [1]with res
	int    _gridNum;					// grid number
	map<int, StecGrid> _stecpergrid;	// stec grids
	unsigned char _satreslevel;			// 2.0.0 add 0~7

private:
};
/* stec model epoch data */
class StecModEpoch 
{
public:
	StecModEpoch()
	{
		_time = { 0 };
		_modetype = 0;
		for (int i = 0; i < NUMSYS; i++) {
			_satNum[i] = 0;
			_refsat[i] = 0;
		}
	};
	~StecModEpoch() {};

	Gtime _time;								// current time
	int	_modetype;								// [0]C00 [1]C00,C01,C10 [2]C00,C01,C10,C11
	int _satNum[NUMSYS];						// sys-sat number
	int _refsat[NUMSYS];						// sys ref sat
	map<int, StecModSat> _stecmodgnss[NUMSYS];	// sys-sat model infos

private:
};

class stecOBS
{
public:
	stecOBS()
	{
		_id = 0;
		_lat = 0.0;
		_lon = 0.0;
		_stec = 0.0;
		_wgt = 0.0;
	}
	~stecOBS() {}
	bool operator == (const stecOBS& src) const {
		return (fabs(src._lat - _lat) < DBL_EPSILON) &&
			   (fabs(src._lon - _lon) < DBL_EPSILON);
	}

	int    _id;
	double _lat;
	double _lon;
	double _stec;
	double _wgt;
private:
};

class StecModel
{
public:
	StecModel()
	{
		_tnow = { 0 };
		_useres = 1;
		_cursys = SYS_NONE;
		for (int i = 0; i < NUMSYS; i++) {
			_fixsys[i] = 0;
			for (int j = 0; j < MAXSAT; j++) {
				_X[i][j] = MatXd::Zero(3, 1);
				_Q[i][j] = MatXd::Zero(3, 3);
			}
		}
		_sysidx[0] = SYS_GPS;
		_sysidx[1] = SYS_GLO;
		_sysidx[2] = SYS_GAL;
		_sysidx[3] = SYS_BDS2;
		_sysidx[4] = SYS_BDS3;
		_refsatsmooth = 0;
		_minel = 0.0;
		for (int i = 0; i < 5; i++) {
			_roti[i] = 0.0;
		}
		_qi_multi = 0.68;
		_qi_base  = 0.0;
		_qi_coeff = 1.0;
	}
	~StecModel() {}

	void settime(IN const Gtime t);
	void setBasicOption(IN ProOption& opt, IN int res);
	bool setCurSys(IN int* usesys, IN int symbol);
	int findSatStec(IN Gtime tnow, IN string site, IN int prn, IN int ref, IN int idx, 
		IN AtmoEpochs& group, OUT map<int, double>& SD);
	int getAtmo(IN Gtime tnow, IN AtmoEpochs& group, OUT AtmoEpoch& atmo);
	int uniformRef(IO AtmoEpoch& atmo);
	bool preCheckSatModel(IN Gtime tnow, IN AtmoEpochs& group, OUT AtmoEpoch& atmo);
	const StecModEpoch* StecModInLastEpoch();
	void initStecMod(IN AtmoEpoch& atmo);
	bool checkSatStatus(IN AtmoEpoch& atmo, IN int sys, IN int prn);
	bool getStecObs(IN AtmoEpoch& atmo, IN int sys, IN int prn, OUT vector<stecOBS>& obss);
	int markUnhalthSites(IO vector<stecOBS>& obss);
	int markUnhalthSitesRes(IN int sys, IN int prn, OUT int& nsta);
	bool estCoeff(IN vector<stecOBS>& obss, IN GridInfo& grid, IN int sys, IN int prn, OUT StecModSat& dat);
	bool oneSatModelEst(IN AtmoEpoch& atmo, IN GridInfo& grid, IN int sys, IN int prn, OUT StecModSat& dat);
	double calcGridTecRes(IN int sys, IN int prn, IN GridEach& grid);
	double calcRovTecRes(IN string site, IN const double* blh, IN GridInfo& grid, IN StecModSat& dat);
	bool oneSatStecRes(IN AtmoEpoch& atmo, IN GridInfo& grid, IN int sys, IN int prn, IO StecModSat& dat);
	bool recalculateQI(IN AtmoEpoch& atmo, IN int sys, IN int prn, IN GridInfo& grid, OUT StecModSat& dat);
	void satModAndRes(IN int id, IN double dlat, IN double dlon, IN const StecModSat& satdat, OUT double& stec, OUT double& res);
	bool checkSatContinuous(IN Gtime tnow, IN int sys, IN int prn, IN int ref, IN GridInfo& grid, IN StecModSat& dat);
	void addStecMod(IN Gtime tnow);
	void refSatSmoothing(IN Gtime tnow);
	void satModEst(IN AtmoEpoch& atmo, IN GridInfo& grid);

	/* basic option */
	Gtime _tnow;							// current time
	int _useres;							// use grid resudial or not [0]no use [1]use
	int _cursys;							// current system *default: SYS_NONE
	int _fixsys[NUMSYS];					// fix system or not *default: [0,0,0,0,0]
	int _sysidx[NUMSYS];					// system index [1,2,4,8,16]
	int _refsatsmooth;						// reference satellite smooth [0]no use [1]use *default: 0
	double _minel;							// minimum elevation angle (deg)
	double _roti[5];						// max roti
	double _qi_multi;						// QI for muti-system
	double _qi_base;						// QI for base
	double _qi_coeff;						// QI for coeffieients
	/* model */
	StecModEpoch _stecModCur;				// current sat stec model
	AtmosRes     _stecModRes;				// sat stec grid residual
	SiteGridDist _siteGridDist;				// distance from each sites to each grids
	map<Gtime, StecModEpoch> _stecModList;	// previous sat stec model (10 epoch)
	set<int> _satList[NUMSYS];				// sat list
	vector<Pos> _badsta;					// bad sites
	StecRoti _stecRoti;						// stec roti
	/* math */
	MatXd _X[NUMSYS][MAXSAT];				// coffeicient matrix
	MatXd _Q[NUMSYS][MAXSAT];				// covariance  matrix
	/* debug */

private:
};

#endif // !STEC_MODEL_H