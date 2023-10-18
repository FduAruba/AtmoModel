#ifndef READFILES_H
#define READFILES_H

#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "io.h"
#include "direct.h"
#include "dirent.h"

#include "comminterface.h"
#include "localmodel.h"

using namespace std;

/* option */
class Coption
{
public:
	Coption()
	{
		for (int i = 0; i < 6; i++) {
			_ts[i] = _te[i] = 0;
		}
		_ti = 0;

		_useres = 0;
		_usesys = 1;
		_fixsys = 0;
		_ressys = 1;
		_maxsatres = 20;
		_refsatsmooth = 0;

		_algotype = 1;
		_modeltype = 1;
		_minel = 15.0;
		_qimulti = 0.68;
		_qibase = 0.0;
		_qicoeff = 1.0;

		for (int i = 0; i < 2; i++) {
			_latgrid[i] = _longrid[i] = 0.0;
			_latcell[i] = _loncell[i] = 0.0;
			_step[i] = 0.0;
			_center[i] = _length[i] = 0.0;
		}

		_nlack = 0;
		_nvali = _ngoodres = _nbadres = _noutres = 0;
	};
	~Coption() {};

	/* file setting */
	std::string _pathin;		// input folder
	std::string _pathou;		// output folder
	std::vector<Cstation> _sta; // stations
	std::vector<Cstation> _rov; // rovers
	/* time setting */
	int _ts[6];					// time begin
	int _te[6];					// time end
	int _ti;					// time interval
	/* system setting */
	int _useres;				// use grid resudial or not [0]no use [1]use
	int _usesys;				// use system or not   [1]GPS [2]GLO [4]GAL [8]BDS2 [16]BDS3
	int _fixsys;				// fix system or not   [1]GPS [2]GLO [4]GAL [8]BDS2 [16]BDS3
	int _ressys;				// use residual or not [1]GPS [2]GLO [4]GAL [8]BDS2 [16]BDS3
	int _maxsatres;				// maximum satellite number use residuals *default: 20
	int _refsatsmooth;			// reference satellite smooth   [0]no use [1]use *default: 0
	/* model setting */
	int _algotype;				// [0]mixed [1]fixed
	int _modeltype;				// [1]STEC [2]ZTD [4]STD [8]VTEC
	double _minel;				// minimum elevation angle (deg)
	double _qimulti;			// QI for muti-system
	double _qibase;				// QI for base
	double _qicoeff;			// QI for coeffieients
	/*  region setting */
	double _latgrid[2];			// lat grid [0]min [1]max (deg)
	double _latcell[2];			// lat cell [0]min [1]max (deg)
	double _longrid[2];			// lon grid [0]min [1]max (deg)
	double _loncell[2];			// lon cell [0]min [1]max (deg)
	double _step[2];			// step     [0]lat [1]lon (deg)
	double _center[2];			// center   [0]lat [1]lon (deg)
	double _length[2];			// length   [0]lat [1]lon (deg)
	/* debug*/
	int _nlack;
	int _nvali;
	int _ngoodres;
	int _nbadres;
	int _noutres;
	map<string, vector<double>> _rovstatic;

private:

};

/* @brief get stations & rovers */
static int getStations(IN set<string> rovset, OUT Coption& config);
/* @brief get file option */
static bool configFile(IN FILE* fp, OUT Coption& config);
/* @brief get time option */
static bool configTime(IN FILE* fp, OUT Coption& config);
/* @brief get system option */
static bool configSystem(IN FILE* fp, OUT Coption& config);
/* @brief get model option */
static bool configModel(IN FILE* fp, OUT Coption& config);
/* @brief get region option */
static bool configRegion(IN FILE* fp, OUT Coption& config);
/* @brief read config file and set option */
extern bool readConfigFile(IN const char* fname, OUT Coption& config);

/* @brief decode epoch data */
static int decodeEpoch(IN Gtime tnow, IN char* line, OUT SiteAtmo& sta, OUT Gtime& tt);
/* @brief sat data */
static int decodeSatData(IN Gtime tnow, IN FILE* fp, IN int ns, OUT SatAtmos* satinfo);
/* @brief entrance of decode data */
static int decodeData(IN Gtime tnow, OUT vector<Cstation>& sta, OUT SiteAtmos& stas);
/* @brief read augment data */
extern bool readAugmentData(IN Gtime tnow, IN Coption& config, OUT SiteAtmos& stas, OUT SiteAtmos& rovs);

/* @brief set option from config */
extern void movOption(IN Coption& config, OUT ProOption& opt);
/* @brief set grids from config */
extern bool movGrids(IN Coption& config, OUT GridInfo& grid);
/* @brief set atmo epoch-wise */
extern bool movAtmos(IN Gtime tnow, SiteAtmos& stas, OUT AtmoInfo& stecinf);

#endif