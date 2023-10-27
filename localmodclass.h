#ifndef LOCAL_MOD_CLASS_H
#define LOCAL_MOD_CLASS_H

#include "comminterface.h"

/* satellite-atmosphere */
class SatAtmo
{
public:
	SatAtmo()
	{
		_sys = "";
		_prn = 0;
		_fixflag = 0;
		_iono = _covIono = 0.0;
		_trop = _covTrop = 0.0;
		_quality = _rot = 0.0;
		for (int i = 0; i < 2; i++) {
			_azel[i] = 0.0;
			_nflo[i] = _nfix[i] = 0.0;
		}
		for (int i = 0; i < 3; i++) {
			_xyz[i] = 0.0;
		}
	};
	~SatAtmo() {};

	/* sat info */
	string _sys;		// sat system (GPS/GLO/GAL/BDS2/BDS3)
	int    _prn;		// PRN number
	int    _fixflag;	// [0]flo [1]fix
	double _azel[2];	// [0]az [1]el (rad)
	double _xyz[3];		// sat ecef pos (m)

	/* sat atmo */
	double _iono;		// STEC (TECu)
	double _covIono;	// sigma of STEC (TECu)
	double _trop;		// Slant Trop Delay (m)
	double _covTrop;	// sigma of STD (m)
	double _quality;	// roti
	double _rot;		// rot

	/* sat ambiguity */
	double _nflo[2];	// [0]N1 [1]N2 (cyc)
	double _nfix[2];	// [0]N1 [1]N2 (cyc)

private:
};
typedef map<int, SatAtmo> SatAtmos;
/* station-atmosphere */
class SiteAtmo
{
public:
	SiteAtmo()
	{
		_name = "";
		_ID = 0;
		for (int i = 0; i < 3; i++) {
			_xyz[i] = _blh[i] = 0.0;
		}
		_zhd = _zwd = _std_zwd = 0.0;
	};
	~SiteAtmo() {};

	/* site info*/
	string _name;	// site name (4-char)
	int    _ID;		// site ID
	double _xyz[3];	// site xyz (m)
	double _blh[3];	// site blh (rad,rad,m)

	/* site atmo */
	double _zhd;	// ZHD (m)
	double _zwd;	// ZWD (m)
	double _std_zwd;// sigma of ZWD (m)
	SatAtmos _satIon[NUMSYS]; // sys->prn->data

private:
};
typedef std::map<int, SiteAtmo> SiteAtmos;
/* station information */
class Cstation
{
public:
	Cstation()
	{
		for (int i = 0; i < 3; i++) {
			_xyz[i] = 0.0;
			_blh[i] = 0.0;
		}
		_fp = NULL;
		_ID = 0;

		_nbad = 0;
		_nepo = 0;
	};
	~Cstation() {};

	std::string _name;
	int _ID;
	double _xyz[3];
	double _blh[3];
	std::string _path;
	std::string _line;
	FILE* _fp;

	/* debug */
	int _nbad;
	int _nepo;

private:
};
/* basic grid */
struct GridEach
{
	int    _id;		// grid ID
	double _lat;	// latitude  (rad)
	double _lon;	// longitude (rad)
};
/* grid information */
class GridInfo
{
public:
	GridInfo()
	{
		_gridNum = 0;
		for (int i = 0; i < 2; i++) {
			_latgrid[i] = _longrid[i] = 0.0;
			_latcell[i] = _loncell[i] = 0.0;
			_step[i] = 0.0;
			_center[i] = _length[i] = 0.0;
		}
		for (int i = 0; i < MAX_GRID; i++) {
			_grids[i] = { 0 };
		}
	}
	~GridInfo() {}

	void deepcopy(const GridInfo& src);

	int    _gridNum;			// number of grids
	double _latgrid[2];			// lat grid [0]min [1]max (deg)
	double _latcell[2];			// lat cell [0]min [1]max (deg)
	double _longrid[2];			// lon grid [0]min [1]max (deg)
	double _loncell[2];			// lon cell [0]min [1]max (deg)
	double _step[2];			// step     [0]lat [1]lon (deg)
	double _center[2];			// center   [0]lat [1]lon (deg)
	double _length[2];			// length   [0]lat [1]lon (deg)
	GridEach _grids[MAX_GRID];	// grid datas

private:
};
/* sat solution */
class SatSol
{
public:
	SatSol()
	{
		_sys = "";
		_prn = 0;
		_fixflag = _locktime = 0;
		_iono = _covIono = 0.0;
		_trop = _covTrop = 0.0;
		_quality = _rot = 0.0;
		for (int i = 0; i < 2; i++) {
			_azel[i] = 0.0;
		}
		for (int i = 0; i < 3; i++) {
			_xyz[i] = 0.0;
		}
	}
	~SatSol() {}

	/* sat info */
	string _sys;		// sat system (GPS/GLO/GAL/BDS2/BDS3)
	int    _prn;		// PRN number
	int    _fixflag;	// [0]flo [1]fix
	int    _locktime;	// locktime
	double _azel[2];	// [0]az [1]el (rad)
	double _xyz[3];		// sat ecef pos (m)

	/* sat atmo */
	double _iono;		// STEC (TECu)
	double _covIono;	// sigma of STEC (TECu)
	double _trop;		// Slant Trop Delay (m)
	double _covTrop;	// sigma of STD (m)
	double _quality;	// roti
	double _rot;		// rot

	/* sat ambiguity */
	//double _nflo[2];	// [0]N1 [1]N2 (cyc)
	//double _nfix[2];	// [0]N1 [1]N2 (cyc)
private:
};
typedef map<int, SatSol> SatSols;
/* site solution */
class SiteSol
{
public:
	SiteSol()
	{
		_time = { 0 };
		_name = "";
		_ID = 0;
		for (int i = 0; i < 3; i++) {
			_xyz[i] = _blh[i] = 0.0;
		}
		_zhd = _zwd = _std_zwd = 0.0;

		for (int i = 0; i < NUMSYS; i++) {
			_satNum[i] = 0;
		}
	}
	~SiteSol() {}

	/* site info*/
	Gtime  _time;			  // solution time
	string _name;			  // site name (4-char)
	int    _ID;				  // site ID
	double _xyz[3];			  // site xyz (m)
	double _blh[3];			  // site blh (rad,rad,m)
	double _zhd;			  // ZHD (m)
	double _zwd;			  // ZWD (m)
	double _std_zwd;		  // sigma of ZWD (m)

	/* site atmo */
	int _satNum[NUMSYS];	  // system-sat numbers
	SatSols _satsols[NUMSYS]; // prn->data

private:
};
typedef map<int, SiteSol> SiteSols;
/* sat info */
class SatInfo
{
public:
	SatInfo()
	{
		_sys = "";
		_prn = 0;
		_fixflag = 0;
		_iono = _covIono = 0.0;
		_trop = _covTrop = 0.0;
		_quality = _rot = 0.0;
		for (int i = 0; i < 2; i++) {
			_azel[i] = 0.0;
			//_nflo[i] = _nfix[i] = 0.0;
		}
		for (int i = 0; i < 3; i++) {
			_xyz[i] = 0.0;
		}
		_wgt = 0.0;
	}
	SatInfo(IN SatSol src);
	~SatInfo() {}

	/* sat info */
	string _sys;		// sat system (GPS/GLO/GAL/BDS2/BDS3)
	int    _prn;		// PRN number
	int    _fixflag;	// [0]flo [1]fix
	double _azel[2];	// [0]az [1]el (rad)
	double _xyz[3];		// sat ecef pos (m)

	/* sat atmo */
	double _iono;		// STEC (TECu)
	double _covIono;	// sigma of STEC (TECu)
	double _trop;		// Slant Trop Delay (m)
	double _covTrop;	// sigma of STD (m)
	double _quality;	// roti
	double _rot;		// rot

	/* sat ambiguity */
	//double _nflo[2];	// [0]N1 [1]N2 (cyc)
	//double _nfix[2];	// [0]N1 [1]N2 (cyc)

	/* math */
	double _wgt;			// weight

private:
};
typedef unordered_map<int, SatInfo> SatInfos;
/* site info */
class SiteInfo
{
public:
	SiteInfo()
	{
		_name.clear();
		_ID = 0;
		for (int i = 0; i < 3; i++) {
			_xyz[i] = _std_xyz[i] = _blh[i] = 0.0;
		}
		for (int i = 0; i < NUMSYS; i++) {
			_satNum[i] = _supSys[i] = 0;
		}
	}
	~SiteInfo() {};
	void operator = (const SiteInfo& src)
	{
		this->_name = src._name;
		this->_ID = src._ID;
		for (int i = 0; i < 3; i++) {
			this->_xyz[i] = src._xyz[i];
			this->_blh[i] = src._blh[i];
			this->_std_xyz[i] = src._std_xyz[i];
		}
		for (int i = 0; i < NUMSYS; i++) {
			this->_satNum[i] = src._satNum[i];
			this->_supSys[i] = src._supSys[i];
		}
	}
	void coypSiteSol(IN SiteSol src);

	/* site info*/
	string _name;			// site name (4-char)
	int    _ID;				// site ID
	double _xyz[3];			// site xyz (m)
	double _std_xyz[3];		// site std xyz (m)
	double _blh[3];			// site blh (rad,rad,m)
	int _satNum[NUMSYS];	// sat number of each system
	int _supSys[NUMSYS];	// [0]unsupport [1]support

private:
};
/* site atmo info */
class AtmoSite
{
public:
	AtmoSite()
	{
		for (int i = 0; i < NUMSYS; i++) {
			_refSat[i] = 0;
		}
		_pppflag = 6;
	}
	AtmoSite(IN SiteSol src);
	~AtmoSite() {}

	int _refSat[NUMSYS];
	int _pppflag;
	SiteInfo _staInfo;
	unordered_map<int, SatInfo> _satInfos[NUMSYS];  // sys->prn->satinfo
};
/* all sites epoch atmo data */
class AtmoEpoch
{
public:
	AtmoEpoch()
	{
		_time = { 0 };
		for (int i = 0; i < 6; i++) {
			_ep[i] = 0;
		}
		_stanum = 0;
		memset(_refSat, 0, sizeof(_refSat));
	}
	AtmoEpoch(IN Gtime t, IN int nsta);
	~AtmoEpoch() {}
	

	Gtime _time;
	int _ep[6];
	int _stanum;
	int _refSat[NUMSYS];

	unordered_map<string, AtmoSite> _staAtmos; // all stations
};
typedef map<Gtime, AtmoEpoch> AtmoEpochs;
/* all sites grid residual data */
class AtmosRes
{
public:
	AtmosRes()
	{
		_time = { 0 };
		for (int i = 0; i < 6; i++) {
			_ep[i] = 0;
		}
		_stanum = 0;
	}
	~AtmosRes() {}

	Gtime _time;
	int _ep[6];
	int _stanum;
	unordered_map<string, AtmoSite> _staRes; // all stations map
};
/* atmo information */
class AtmoInfo
{
public:
	AtmoInfo()
	{
		_time = { 0 };
		_sitenum = 0;
	}
	~AtmoInfo() {}
	void reset();

	Gtime _time;		// solution time
	int _sitenum;		// site number
	SiteSols _sitesols; // site solutions

private:
};
/* process stec residual grid */
class ProStecGrid
{
public:
	ProStecGrid()
	{
		_gridid = 0;
		_nsta = 0;
		_lat = _lon = _stec = 0.0;
		_gridreslevel = 0;
	}
	void deepcopy(IN ProStecGrid& src);
	void reset();
	~ProStecGrid() {}

	int    _gridid;				 // grid ID
	int    _nsta;				 // number of sta res
	double _lat;				 // latitude  (rad)
	double _lon;				 // longitude (rad)
	double _stec;				 // stecc (TECu)
	unsigned char _gridreslevel; // 1.9.0 add 0~7

private:
};
/* process sat stec model */
class ProStecModSat
{
public:
	ProStecModSat()
	{
		_sys = '\0';
		_sat = 0;
		_nsta = 0;
		_ele = 0.0;
		for (int i = 0; i < 4; i++) {
			_coff[i] = _coff_res[i] = 0.0;
		}
		for (int i = 0; i < 2; i++) {
			_QI[i] = 0;
		}
		_gridNum = 0;
		_satreslevel = 0;
	}
	void deepcopy(IN ProStecModSat& src);
	void reset();
	~ProStecModSat() {}

	char   _sys;						// system captial GREC
	int    _sat;						// sat number
	int    _nsta;						// number of stations
	double _ele;						// elevation
	double _coff[4];					// coffeicient C00 C01 C10 C11
	double _coff_res[4];				// coffeicient residual
	double _QI[2];						// [0]withput res [1]with res
	int    _gridNum;					// grid number
	ProStecGrid _stecpergrid[MAX_GRID]; // stec grids
	unsigned char _satreslevel;			// 2.0.0 add 0~7
private:
};
/* process stec moddel */
class ProStecMod
{
public:
	ProStecMod()
	{
		_time = { 0 };
		_modetype = 1;
		for (int i = 0; i < NUMSYS; i++) {
			_satNum[i] = _satRef[i] = 0;
		}
	}
	~ProStecMod() {}
	void reset();

	Gtime _time;								// current time
	int _modetype;								// [0]C00 [1]C00,C01,C10 [2]C00,C01,C10,C11
	int _satNum[NUMSYS];						// sys-sat number
	int _satRef[NUMSYS];						// sys ref sat
	map<int, ProStecModSat> _stecmod[NUMSYS];	// sys-sat model infos

private:
};
/* process option */
class ProOption
{
public:
	ProOption()
	{
		_ti = 0;
		for (int i = 0; i < NUMSYS; i++) {
			_usesys[i] = _fixsys[i] = _ressys[i] = 0;
		}
		_maxsatres = 20;
		_algotype = 1;
		_fittype = 0;
		_diffmode = 1;
		_refsatsmooth = 0;
		_minel = 15.0;
		for (int i = 0; i < 5; i++) {
			_maxroti[i] = 0.0;
		}
		_qimulti = 0.68;
		_qibase = 0.0;
		_qicoeff = 1.0;
	}
	~ProOption() {}

	void deepcopy(IN ProOption& src);

	int _ti;					// time interval
	int _usesys[NUMSYS];		// use system or not [0]no use [1]use
	int _fixsys[NUMSYS];		// fix system or not [0]no use [1]use
	int _ressys[NUMSYS];		// use res or not    [0]no use [1]use
	int _maxsatres;				// maximum satellite number use residuals *default: 20
	int _algotype;				// innput data type [0]mixed [1]fixed
	int _fittype;				// fitting model [0]IDW [1]MSF
	int _diffmode;				// difference mode [1]SD [2]SD+UD
	int _refsatsmooth;			// reference satellite smooth [0]no use [1]use *default: 0
	double _minel;				// minimum elevation angle (deg)
	double _maxroti[5];			// roti of stec *default: {0.5,0,0,0,0}
	double _qimulti;			// QI for muti-system
	double _qibase;				// QI for base
	double _qicoeff;			// QI for coeffieientss

private:
};
/* distance from site to grid */
using SiteDist = unordered_map<string, double>;
class SiteGridDist
{
public:
	void calcAllDist(IN AtmoEpoch& atmo, IN GridInfo& grid);
	double getDist(IN int id, IN string site);

	unordered_map<int, SiteDist> _dist;	// gridid->site->dist
private:
};

class SatRoti
{
public:
	SatRoti()
	{
		_prn = 0;
		_roti = 0.0;
		_tlast = { 0 };
	}
	SatRoti(IN Gtime tnow, IN int prn)
	{
		_prn = prn;
		_roti = 0.0;
		_tlast = tnow;
	}
	~SatRoti() {}
	void reinit(IN Gtime tnow, IN double stec);
	double proRoti(IN Gtime tnow, IN SatAtmo& satinf);

	int _prn;
	double _roti;
	Gtime _tlast;
	map<Gtime, double> _stec;
	map<Gtime, double> _rot;
private:
};

class StaRoti
{
public:
	StaRoti()
	{
		_site = "";
	}
	StaRoti(IN string site)
	{
		_site = site;
	}
	~StaRoti() {}
	void proRoti(IN Gtime tnow, IN SiteAtmo& siteinf);
	
	string _site;
	map<int, SatRoti> _rotis[NUMSYS];
private:
};
typedef map<string, StaRoti> StaRotiMap;

#endif // !LOCAL_MOD_CLASS_H

