#include "stecmodel.h"

#define		THRESHOLD_MAX_EXPAND_TIME	(120.0)

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

int StecModel::findSatStec(IN Gtime tnow, IN string site, IN int prn, IN int ref, IN int idx,
	IN AtmoEpochs& group, OUT map<int, double>& SD)
{
	for (auto pAtmo = group.rbegin(); pAtmo != group.rend(); ++pAtmo) {
		double dt = tnow - pAtmo->first;
		if (dt < 0.0) {
			continue;
		}
		if (dt > THRESHOLD_MAX_EXPAND_TIME) {
			break;
		}

		auto pSta = pAtmo->second._staAtmos.find(site);
		if (pSta == pAtmo->second._staAtmos.end()) {
			continue;
		}
		auto pSat = pSta->second._satInfos[idx].find(prn);
		auto pRef = pSta->second._satInfos[idx].find(ref);
		if (pSat == pSta->second._satInfos[idx].end()) {
			//printf("%s %s sat %c%02d not found\n", strtime(tnow, 2).c_str(), site.c_str(), idx2sys(idx), prn);
			continue;
		}
		if (pRef == pSta->second._satInfos[idx].end()) {
			//printf("%s %s ref %c%02d not found\n", strtime(tnow, 2).c_str(), site.c_str(), idx2sys(idx), ref);
			continue;
		}

		double dstec = pSat->second._iono - pRef->second._iono;

		SD.emplace((int)dt, dstec);
	}

	return (int)SD.size();
}

bool StecModel::preCheckSatModel(IN Gtime tnow, IN AtmoEpochs& group, OUT AtmoEpoch& atmo)
{	
	int cnt = 0;
	AtmoEpochs::iterator pAtmo;

	_stecRoti.procRoti(tnow, group);

	if ((pAtmo = group.find(tnow)) == group.end()) {
		return false;
	}

	atmo._stanum = 0;
	atmo._time = tnow;
	time2epoch(tnow, atmo._ep);
	for (int isys = 0; isys < NUMSYS; isys++) {
		atmo._refSat[isys] = pAtmo->second._refSat[isys];
	}

	for (auto pSta : pAtmo->second._staAtmos) {
		cnt++;
		AtmoSite tmpSite;
		tmpSite._staInfo._name = pSta.first;
		for (int isys = 0; isys < NUMSYS; isys++) {
			tmpSite._staInfo._supSys[isys] = pSta.second._staInfo._supSys[isys];
		}
		tmpSite._staInfo._ID = pSta.second._staInfo._ID;
		for (int i = 0; i < 3; i++) {
			tmpSite._staInfo._xyz[i] = pSta.second._staInfo._xyz[i];
			tmpSite._staInfo._blh[i] = pSta.second._staInfo._blh[i];
			tmpSite._staInfo._std_xyz[i] = pSta.second._staInfo._std_xyz[i];
		}
		for (int isys = 0; isys < NUMSYS; isys++) {
			if (!(_cursys & _sysidx[isys])) {
				continue;
			}

			int nsat = 0;
			for (const auto& pSat : _satList[isys]) {
				int prn = pSat;
				int ref = pAtmo->second._refSat[isys];

				int nep = 0;
				map<int, double> SD_Stec;
				nep = findSatStec(tnow, pSta.first, prn, ref, isys, group, SD_Stec);
				if (nep <= 0) {
					//printf("%s %c%02d %c%02d\n", pSta.first.c_str(), idx2sys(isys), prn, idx2sys(isys), ref);
					continue;
				}

				auto pSat = pSta.second._satInfos[isys].find(prn);
				if (pSat == pSta.second._satInfos[isys].end()) {
					continue;
				}

				tmpSite._satInfos[isys].emplace(prn, pSat->second);
				nsat++;
			}
			tmpSite._staInfo._satNum[isys] = nsat;
		}
		atmo._staAtmos.emplace(pSta.first, tmpSite);
	}

	for (auto& pSta : atmo._staAtmos) {
		for (int isys = 0; isys < NUMSYS; isys++) {
			if (!(_cursys & _sysidx[isys])) {
				continue;
			}


		}
	}


	return true;
}

const StecModEpoch* StecModel::StecModInLastEpoch()
{
	return _stecModList.size() > 0 ? &_stecModList.rbegin()->second : NULL;
}