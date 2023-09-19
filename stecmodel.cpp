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

int StecModel::getAtmo(IN Gtime tnow, IN AtmoEpochs& group, OUT AtmoEpoch& atmo)
{
	int cnt = 0;
	AtmoEpochs::iterator pAtmo;
	if ((pAtmo = group.find(tnow)) == group.end()) {
		return 0;
	}

	atmo._stanum = 0;
	atmo._time = tnow;
	time2epoch(tnow, atmo._ep);
	for (int isys = 0; isys < NUMSYS; isys++) {
		atmo._refSat[isys] = pAtmo->second._refSat[isys];
	}

	for (auto pSta : pAtmo->second._staAtmos) {
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
		cnt++;
	}

	return cnt;
}

int StecModel::uniformRef(IO AtmoEpoch& atmo)
{
	for (auto& pSta : atmo._staAtmos) {
		bool buse = false;
		for (int isys = 0; isys < NUMSYS; isys++) {
			if (!(_cursys & _sysidx[isys])) {
				continue;
			}
			if (pSta.second._satInfos[isys].size() == 0) {
				continue;
			}

			auto pRef = pSta.second._satInfos[isys].find(atmo._refSat[isys]);
			if (pRef == pSta.second._satInfos[isys].end() || atmo._refSat[isys] == 0) {
				pSta.second._staInfo._satNum[isys] = 0;
				pSta.second._satInfos[isys].clear();
			}
			else {
				double refion = pRef->second._iono;
				for (auto& pSat : pSta.second._satInfos[isys]) {
					pSat.second._iono -= refion;
				}
				buse = true;
			}
		}
		atmo._stanum += buse ? 1 : 0;
	}

	return atmo._stanum;
}

bool StecModel::preCheckSatModel(IN Gtime tnow, IN AtmoEpochs& group, OUT AtmoEpoch& atmo)
{
	/* calculate ROTI */
	_stecRoti.procRoti(tnow, group);
	/* get current atmo data */
	if (!getAtmo(tnow, group, atmo)) { return false; }
	/* uniform the ref sat and recovery the SD STEC */
	if (!uniformRef(atmo)) { return false; }

	return true;
}

const StecModEpoch* StecModel::StecModInLastEpoch()
{
	return _stecModList.size() > 0 ? &_stecModList.rbegin()->second : NULL;
}

void StecModel::initStecMod(IN AtmoEpoch& atmo)
{
	/* skip GLONASS */
	if (_cursys & SYS_GLO && _stecModRes._stanum > 0) {
		return;
	}

	/* init stecmod */
	for (int isys = 0; isys < NUMSYS; isys++) {
		_stecModCur._satNum[isys] = 0;
		_stecModCur._refsat[isys] = atmo._refSat[isys];
		_stecModCur._stecmodgnss[isys].clear();
	}
	_stecModCur._modetype = 1;
	_stecModCur._time = atmo._time;

	/* init stecres */
	_stecModRes._stanum = atmo._stanum;
	time2epoch(atmo._time, _stecModRes._ep);
	_stecModRes._time = atmo._time;
	_stecModRes._staRes.clear();

	return;
}

void StecModel::satModEst(IN AtmoEpoch& atmo, IN GridInfo& grid)
{
	/* 1.计算基线距离 */
	_siteGridDist.calcAllDist(atmo, grid);

	/* 2.参数/残差估计 */


}