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

int StecRoti::pushTEC(IN AtmoEpochs& group, OUT set<string>* list)
{
	for (auto& it : group) {
		Gtime ep = it.first;
		auto& sites = it.second;

		for (auto& atmo : sites._staAtmos) {
			string isite = atmo.first;
			AtmoSite& dat = atmo.second;
			int cnt = 0;

			for (int isys = 0; isys < NUMSYS; isys++) {
				for (auto& ie : dat._satInfos[isys]) {
					int prn = ie.first;

					RotKey key(isite, idx2sys(isys), prn);
					_siteSatStec[key].emplace(ep, ie.second._iono);
					cnt++;
				}
			}
			if (cnt > 0 && list) { list->emplace(isite); }
		}
	}

	return (int)_siteSatStec.size();
}

int StecRoti::calcROT()
{
	int cnt = 0;
	
	for (auto& it : _siteSatStec) {
		if (it.second.size() < 2) { continue; }

		int isys = satid2idx(it.first._system, it.first._prn);
		int prn = it.first._prn;

		auto in = it.second.end();
		for (auto ie = it.second.begin();ie != it.second.end();++ie) {
			in = ie; in++;
			if (in != it.second.end()) {
				double diffsec = in->first - ie->first;
				if (diffsec > 0.0) {
					double rot = (in->second - ie->second) * 60.0 / diffsec; // TECu/min
					_cellSatsRot[isys][prn]._rotArray.push_back(rot);
					_cellSatsRot[isys][prn]._sum += rot;
					cnt++;
				}
			}
		}
	}

	return cnt;
}

bool StecRoti::procRoti(IN Gtime tnow, IN AtmoEpochs& group)
{
	set<string> siteList;

	StecRoti::reset();
	if (group.size() < 2) {
		return true;
	}

	if (!pushTEC(group, &siteList)) {
		return false;
	}

	if (!calcROT()) {
		return false;
	}

	return true;
}