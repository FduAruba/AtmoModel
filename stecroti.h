#ifndef STEC_ROTI_H
#define STEC_ROTI_H

#include <string>
#include <map>

#include "comminterface.h"
#include "localmodclass.h"

typedef map<Gtime, double> StecEpoch;

class RotKey
{
public:
	
	RotKey()
	{
		_site = "";
		_system = '\0';
		_prn = 0;
	}
	RotKey(IN string site, IN char system, IN int prn);
	~RotKey() {}
	
	bool operator < (const RotKey& src) const {
		if (src._site == _site) {
			if (src._system == _system) {
				return _prn < src._prn;
			}
			else {
				return _system < src._system;
			}
		}
		else {
			return _site < src._site;
		}
	}

	string _site;
	char _system;
	int _prn;
private:
};

class SatRot
{
public:
	SatRot()
	{
		_roti = _sum = 0.0;
	}
	~SatRot() {}

	double _roti;
	double _sum;
	vector<double> _rotArray;
private:
};
typedef map<int, SatRot> SatRots;

class StecRoti
{
public:
	StecRoti() {}
	~StecRoti() {}

	void reset();
	int pushTEC(IN AtmoEpochs& group, OUT set<string>* list);
	int calcROT();
	bool procRoti(IN Gtime tnow, IN AtmoEpochs& group);

	map<RotKey, StecEpoch> _siteSatStec;
	SatRots _cellSatsRot[NUMSYS];
private:
};
#endif // !STEC_ROTI_H

