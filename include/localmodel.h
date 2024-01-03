#ifndef LOCAL_MODEL_H
#define LOCAL_MODEL_H

#include "comminterface.h"
#include "localmodclass.h"
#include "stecmodel.h"
#include "ztdmodel.h"

class EleStation
{
public:
	EleStation()
	{
		_prn = _num = 0;
		_el = _rot = 0.0;
	}
	~EleStation() {}

	void operator = (const EleStation& src) {
		this->_prn = src._prn;
		this->_num = src._num;
		this->_el  = src._el;
		this->_rot = src._rot;
	}

	bool operator < (const EleStation& src) const {
		return this->_el < src._el;
	}

	int _prn;		// prn
	int _num;		// site numbers
	double _el;		// elevation (rad)
	double _rot;	// rot
private:
};

class LocalAtmoModel
{
public:
	LocalAtmoModel() 
	{
		_stanumGEC = _stanumR = 0;
		_nbadroti = 0;
	}
	~LocalAtmoModel() {}

	void setBasic(IN int res);
	void setReigon(IN GridInfo& grid);
	void setOption(IN ProOption& opt);
	bool setRefSites(IN SiteAtmos& stas);

	bool getAtmoEpoch(IN Gtime tnow, IN int symbol, OUT AtmoEpochs::iterator& it);

	int stationSys(IN SiteSol& sol);
	void cntFixSat(IN SatSols& sol, OUT map<int, int>& track);
	int fixDecision(IN AtmoInfo& stecinf, OUT map<int, int>* track);
	void sortSiteNumber(IN int isys, IN map<int, EleStation>* numstas, OUT map<int, set<EleStation>>& vprn);
	void findRefSat(IN int isys, IN map<int, EleStation>* numstas, IN map<int, set<EleStation>>& vprn, OUT int* refsat);
	void uniformDatum(IN Gtime tnow, IN int* refsat, IN int symbol);
	bool selectOneRefSat(IN map<int, EleStation>* numstas, OUT int* refsat);
	void inputAtmoSat(IN int isys, IN SatSols* src, OUT SatInfos* dst, OUT map<int, EleStation>* numstas);
	bool inputAtmoSite(IN Gtime tnow, IN AtmoInfo& stecinf, OUT map<int, int>* track, OUT map<int, EleStation>* numstas);
	bool inputAtmoEpoch(IN Gtime tnow, IN AtmoInfo& stecinf, IN bool selectRef);
	bool doStecModSys(IN int symbol);
	void resetRefSatBD2(IO StecModEpoch& mod);
	void setSatResLevel(IO StecModEpoch& mod);
	void copyStecSat(IN int sys, IN int ref, IN map<int, StecModSat>& src, OUT map<int, ProStecModSat>& dst);
	void copyStecMod(IN StecModEpoch& src, OUT ProStecMod& dst);
	bool doStecMod(IN Gtime tnow, IN AtmoInfo& stecinf, OUT ProStecMod& stecmod);

	int inputZtdEpoch(IN Gtime tnow, IN AtmoInfo& ztdinf);
	int sparseSites();
	void copyZtdMod(IN int npara, IN ZtdModEpoch& src, OUT ProZtdMod& dst);
	bool doZtdMod(IN Gtime tnow, IN AtmoInfo& ztdinf, OUT ProZtdMod& ztdmod);
	
	int _stanumGEC;							// station number GEC
	int _stanumR;							// station number R

	ProOption _proOption;					// process option
	StecModel _stecPro;						// stec process model
	ZtdModel  _ztdPro;						// ztd process model

	GridInfo              _gridinfo;		// basic grid info
	map<int, SiteInfo>    _allsites;		// all ref sites info

	AtmoEpochs _stecGroupGEC;	// stec group GEC
	AtmoEpochs _stecGroupR;		// stec group R

	ZtdInfos _ztds;
	/* debug */
	int _nbadroti;

private:
};


#endif // !LOCAL_MODEL_H

