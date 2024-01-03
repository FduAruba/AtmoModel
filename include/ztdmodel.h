#ifndef ZTD_MODEL_H
#define ZTD_MODEL_H

#include "comminterface.h"
#include "localmodclass.h"
#include "algmatrix.h"

class ZtdModEpoch
{
public:
	ZtdModEpoch() {
		_time = { 0 };
		_valid = false;
		_nsta = _zhd = 0;
		_qi = 0;
		for (int i = 0; i < 10; i++) {
			_coeff[i] = _coeff_rms[i] = 0;
		}
	}
	~ZtdModEpoch() {}

	Gtime _time;								// current time
	int _nsta;									// number of model stations
	bool _valid;								// valid flag
	double _zhd;								// region mean zenith hydrostatic delay(m)
	double _coeff[10];							// zwd model coefficients
	double _coeff_rms[10];						// zwd model coefficient resduials
	double _qi;									// zwd model rms(m)
private:
};

class ZtdModel
{
public:
	ZtdModel() {
		_tnow = { 0 };
		_useres = _fittype = 0;
		_ncoeff = 0;
	}
	~ZtdModel() {}
	void setBasicOption(IN ProOption& opt, IN int res);
	void settime(IN const Gtime t);
	bool ztdModSaas(IN double* pos, IN double* azel, IN double humi, OUT double* ztd);
	void ztd2msl(IN int bmean, IO ZtdInfos& ztds);
	bool ofcModel(IN Gtime t, IN ZtdInfos& ztds, IN GridInfo& grids, IN int bmean);
	VecXd buildRowA(IN int nx, IN int np, IN double* blh);
	int bulidMatrix(IN int nx, IN ZtdInfos& ztds, IN double* c, IN unordered_set<string> ext,
		OUT VecXd* y, OUT VecXd* P, OUT MatXd* A, OUT map<int, string>& idx);
	double estres(int nx, int ny, VecXd* x, VecXd* y, MatXd* A, VecXd* v);

	/* basic option */
	Gtime _tnow;							// current time
	int _useres;							// use grid resudial or not [0]no use [1]use
	int _fittype;							// fitting model [0]OFC [1]OFC-MSL [2]MOFC
	int _ncoeff;							// number of coefficients
	/* model */
	ZtdModEpoch _ztdModCur;					// current ztd model
	set<string> _badsta;					// bad sites

private:
};




#endif // !ZTD_MODEL_H

