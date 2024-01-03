#ifndef COMM_FUNCTION_H
#define COMM_FUNCTION_H

#include "const.h"
#include "commtime.h"
#include "comminterface.h"

typedef Eigen::MatrixXd MatXd;
typedef Eigen::VectorXd VecXd;
typedef map<string, map<int, FILE*>> FileFps;
typedef std::vector<double> Dvec;

struct StaDistIon
{
	StaDistIon() = default;
	StaDistIon(double d, double i) :dist(d), ion(i) {}
	bool operator < (const StaDistIon& other) const {
		return dist < other.dist;
	}

	double dist = 0.0;
	double ion = 0.0;
};
typedef vector<StaDistIon> StaDistIonArr;

class StaDistMSF
{
public:
	StaDistMSF() {
		_site = "";
		_lat = _lon = 0.0;
		_g0i = _ion = 0.0;
		_gij = Eigen::VectorXd::Zero(6);
	}
	StaDistMSF(string site, double lat, double lon, double dist, double ion) {
		_site = site;
		_lat = lat;
		_lon = lon;
		_g0i = dist;
		_ion = ion;
		_gij = Eigen::VectorXd::Zero(6);
	}
	~StaDistMSF() {}
	bool operator < (const StaDistMSF& src) const {
		return _g0i < src._g0i;
	}

	double _g0i;
	string _site;
	double _lat;
	double _lon;
	double _ion;
	Eigen::VectorXd _gij;
private:
};
typedef vector<StaDistMSF> StaDistMSFArr;

/* string */
extern void modifyPath(IO char* src);
extern void trimSpace(IO char* src);
extern void cutSep(IO char* str);
extern double str2num(IN char* s, IN int i, IN int n);

/* position */
extern void pos2ecef(IN double* pos, OUT double* r);
extern double sphereDist(IN double latG, IN double lonG, IN double latB, IN double lonB);
extern bool gridVaild(IN double latG, IN double lonG, IN double lat, IN double lon, IN double dlat, IN double dlon);

/* sat system */
extern bool checksys(IN char s);
extern string sysstr(IN int sys);
extern char syschar(IN string system);
extern char idx2sys(IN int idx);
extern int satsys(IN int sat, OUT int* prn);
extern int satno(IN int sys, IN int prn);
extern int satid2no(IN char* id);
extern int satid2idx(IN char code, IN int num);
extern double satfreq(IN int sys, IN int f);

/* gps time */
extern string strtime(IN Gtime t, IN int opt);
extern void time2str(IN Gtime t, IN char* s, IN int opt);

/* math */
extern void calcMeanStd(IN vector<double> data, OUT double& vmean, OUT double& vstd);
extern void calcMeanStd(IN VecXd data, OUT double& vmean, OUT double& vstd);
extern double robust(IN double V, OUT double rms);
extern double modelIDW(IN StaDistIonArr& list, IN int nused, IN double maxdist, IN int k, OUT int* n);
extern double modelMSF(IN StaDistMSFArr& list, IN int sz);
extern void polynomial2(IN double x, IN double y, IN int order, IN int idmax, IO int& id, OUT VecXd& vec);
extern void polynomial3(IN double x, IN double y, IN double z, IN int order, IN int idmax, IO int& id, OUT VecXd& vec);

#endif
