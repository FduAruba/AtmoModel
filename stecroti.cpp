#include "stecroti.h"

RotKey::RotKey(IN string site, IN char system, IN int prn)
{
	_site = site;
	_system = system;
	_prn = prn;
}

void StecRoti::reset()
{
	_siteSatStec.clear();
	for (int i = 0;i < NUMSYS;i++) {
		_cellSatsRot[i].clear();
	}
}

int StecRoti::pushTEC(IN Gtime t, IN string site, IN AtmoSite& siteatmo)
{
	for (int isys = 0;isys < NUMSYS;isys++) {
		for (auto& it : siteatmo._satInfos[isys]) {
			RotKey key(site, idx2sys(isys), it.first);
			_siteSatStec[key].emplace(t, it.second._iono);
		}
	}
	return (int)_siteSatStec.size();
}

int StecRoti::calcROT()
{
	for (auto& it : _siteSatStec) {
		if (it.second.size() < 2) { continue; }

		int isys = satid2idx(it.first._system, it.first._prn);
		int prn = it.first._prn;

		auto in = it.second.end();
		for (auto ie = it.second.begin();ie != it.second.end();++ie) {
			in++ = ie;
			if (in != it.second.end()) {
				double diffsec = in->first - ie->first;
				if (diffsec > 0.0) {
					double rot = (in->second - ie->second) * 60.0 / diffsec; // TECu/min
					_cellSatsRot[isys][prn]._rotArray.push_back(rot);
					_cellSatsRot[isys][prn]._sum += rot;
				}
			}
		}
	}

	return (int)_cellSatsRot->size();
}

bool StecRoti::procRoti(IN Gtime tnow, IN AtmoEpochs& group)
{
	set<string> siteList;

	StecRoti::reset();
	if (group.size() < 2) {
		return true;
	}

	for (auto& it : group) {
		auto& sites = it.second;

		for (auto& atmo : sites._staAtmos) {
			if (pushTEC(it.first, atmo.first, atmo.second)) {
				siteList.emplace(atmo.first);
			}
		}
	}

	return true;
}