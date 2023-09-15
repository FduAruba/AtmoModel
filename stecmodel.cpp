#include "stecmodel.h"

void StecModel::setBasicOption(IN ProOption& opt, IN int res)
{
	_useres   = res;
	_minel    = opt._minel;
	_qi_multi = opt._qimulti;
	_qi_base  = opt._qibase;
	_qi_coeff = opt._qicoeff;
	for (int i = 0; i < 5; i++) {
		_roti[i] = opt._maxroti[i];
	}
	for (int i = 0; i < NUMSYS; i++) {
		_fixsys[i] = opt._fixsys[i];
	}
}

void StecModel::setCurSys(IN int* usesys, IN int symbol)
{
	switch (symbol)
	{
	case 0: {
		_cursys |= usesys[IDX_GPS]  ? SYS_GPS  : 0x00;
		_cursys |= usesys[IDX_GAL]  ? SYS_GAL  : 0x00;
		_cursys |= usesys[IDX_BDS2] ? SYS_BDS2 : 0x00;
		_cursys |= usesys[IDX_BDS3] ? SYS_BDS3 : 0x00;
		break;
	}
	case 1: {
		_cursys |= usesys[IDX_GLO]  ? SYS_GLO  : 0x00;
		break;
	}
	default: { break;}
	}
}

bool StecModel::preCheckSatModel(IN Gtime tnow, IN AtmoEpochs& group, OUT AtmoEpoch& atmo)
{	
	_stecRoti.procRoti(tnow, group);

	return true;
}

const StecModEpoch* StecModel::StecModInLastEpoch()
{
	return _stecModList.size() > 0 ? &_stecModList.rbegin()->second : NULL;
}