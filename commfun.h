#ifndef COMM_FUNCTION_H
#define COMM_FUNCTION_H

#include "const.h"
#include "commtime.h"

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

/* string */
extern void modifyPath(IO char* src);
extern void trimSpace(IO char* src);
extern void cutSep(IO char* str);
extern double str2num(IN char* s, IN int i, IN int n);

/* position */
extern void pos2ecef(IN double* pos, OUT double* r);
extern double sphereDist(IN double latG, IN double lonG, IN double latB, IN double lonB);

/* sat system */
extern bool checksys(IN char s);
extern string sysstr(IN int sys);
extern char syschar(IN string system);
extern char idx2sys(IN int idx);
extern int satsys(IN int sat, OUT int* prn);
extern int satno(IN int sys, IN int prn);
extern int satid2no(IN char* id);
extern int satid2idx(IN char code, IN int num);

/* gps time */
extern string strtime(IN Gtime t, IN int opt);
extern void time2str(IN Gtime t, IN char* s, IN int opt);

/* math */
extern void calcMeanStd(IN vector<double> data, OUT double& vmean, OUT double& vstd);
extern double robust(IN double V, OUT double rms);
extern double modelIDW(IN StaDistIonArr& list, IN int nused, IN double maxdist, IN int k);

#endif
