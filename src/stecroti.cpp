#include "stecroti.h"

RotKey::RotKey(IN string site, IN char system, IN int prn)
{
	_site = site;
	_system = system;
	_prn = prn;
}

double SatRot::calcROTI()
{
	_roti = 0.0;
	int MAXITER = 10, rotSize = 0;

	for (int i = 0; i < MAXITER; i++) {
		bool brk = true;
		if ((rotSize = (int)_rotArray.size()) == 0) {
			return 0.0;
		}
		double avgRot = _sum / rotSize; // ¾ùÖµ
		double argSum = 0.0;
		for (auto& it : _rotArray) {
			argSum += std::pow((it - avgRot), 2);
		}
		double STD = std::sqrt(argSum / rotSize); // ·½²î

		//for (auto it = _rotArray.begin(); it != _rotArray.end();) {
		//	double diff = fabs(*it - avgRot);
		//	if (diff > 2 * STD) {
		//		_sum -= (*it);
		//		it = _rotArray.erase(it);
		//		brk = false;
		//	}
		//	else {
		//		++it;
		//	}
		//}

		if (brk) {
			_roti = STD;
			break;
		}
	}
	return _roti;
}

void StecRoti::reset()
{
	_siteSatStec.clear();
	for (int i = 0;i < NUMSYS;i++) {
		_cellSatsRot[i].clear();
	}
	_badroti = 0;
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

int StecRoti::calcROT(IN Gtime tnow)
{
	int cnt = 0, bwrite = 0, cjmp = 0;
	
	for (auto& it : _siteSatStec) {
		if (it.second.size() < 2) { continue; }

		int isys = satid2idx(it.first._system, it.first._prn);
		int prn = it.first._prn;

		auto in = it.second.end();
		for (auto ie = it.second.begin(); ie != it.second.end(); ++ie) {
			in = ie; in++;
			if (in != it.second.end()) {
				double diffsec = in->first - ie->first;
				if (diffsec > 0.0) {
					double rot = (in->second - ie->second) * 60.0 / diffsec; // TECu/min
					/*if (fabs(rot) > 20.0 && in->first >= tnow) {
						printf("%s %c%02d %4s ROT=%6.2f\n", strtime(ie->first, 2).c_str(), it.first._system, 
							prn, it.first._site.c_str(), rot);
						cjmp++;
						bwrite = 1;
					}*/
					_cellSatsRot[isys][prn]._rotArray.push_back(rot);
					_cellSatsRot[isys][prn]._sum += rot;
					cnt++;
				}
			}
		}
	}
	/*if (bwrite && cjmp > 50) {
		printf("%s %d\n", strtime(tnow, 2).c_str(), cjmp);
	}*/
	
	return cnt;
}

int StecRoti::calcROTI(IN Gtime tnow)
{
	int cnt = 0;

	for (int isys = 0; isys < NUMSYS; isys++) {
		for (auto& it : _cellSatsRot[isys]) {
			double roti = it.second.calcROTI();
			if (roti != 0.0) {
				if (roti > 4.0) {
					_badroti++;
				}
				cnt++;
			}
		}
	}

	return cnt;
}

bool StecRoti::procRoti(IN Gtime tnow, IN AtmoEpochs& group)
{
	set<string> siteList;

	StecRoti::reset();

	if (group.size() < 2) { return true; }
	// push tec data
	if (!pushTEC(group, &siteList)) { return false; }
	// calculate ROT
	if (!calcROT(tnow)) { return false; }
	// calculate ROTI
	if (!calcROTI(tnow)) { return false; }

	return true;
}