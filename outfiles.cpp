#include "outfiles.h"

bool checkRef(IN int sys, IN int ref, IN SiteAtmo& roviono)
{
	if (ref == 0) {
		return false;
	}
	if (roviono._satIon[sys].find(ref) == roviono._satIon[sys].end()) {
		return false;
	}
	const auto& pRef = roviono._satIon[sys][ref];
	if (fabs(pRef._iono) < DBL_EPSILON) {
		return false;
	}
	if (pRef._fixflag != 0 && pRef._fixflag != 1) {
		return false;
	}

	return true;
}

bool checkSat(IN int sys, IN int ref, IN int prn, IN SiteAtmo& roviono, IN double minEl)
{
	if (prn == ref) {
		return false;
	}
	const auto& pSat = roviono._satIon[sys][prn];
	const auto& pRef = roviono._satIon[sys][ref];
	if (pSat._azel[1] < minEl) {
		return false;
	}
	if (pSat._fixflag != 0 && pSat._fixflag != 1) {
		return false;
	}
	if (sys != IDX_GLO && (pSat._fixflag != pRef._fixflag)) {
		return false;
	}
	if (fabs(pSat._iono) < DBL_EPSILON || pSat._covIono > CUT_COV_STEC[sys]) {
		return false;
	}
	//TODO: bamboo数据引入ROT/ROTI时启用
	/*if (fabs(pSat._quality) < DBL_EPSILON || fabs(pSat._quality) > 0.2 || fabs(pSat._rot) > 0.8) {
		return false;
	}*/
	
	return true;
}

bool findStecModelSat(IN int sys, IN int prn, IN ProStecMod& stecmod, OUT ProStecModSat& satmod)
{
	if (stecmod._stecmod[sys].find(prn) != stecmod._stecmod[sys].end()) {
		satmod.deepcopy(stecmod._stecmod[sys][prn]);
		return true;
	}
	return false;
}

void rovStecDiff(IN Coption& cfg, IN GridInfo& grid, IN SiteAtmo& roviono, IN ProStecMod& stecmod)
{
	double thresEl = cfg._minel * D2R;
	double latO = cfg._center[0] * D2R;
	double lonO = cfg._center[1] * D2R;
	
	for (int isys = 0; isys < NUMSYS; isys++) {
		int ref = stecmod._satRef[isys];
		if (!checkRef(isys, ref, roviono)) {
			continue;
		}

		for (const auto& pSat : roviono._satIon[isys]) {
			int prn = pSat.first;
			if (!checkSat(isys, ref, prn, roviono, thresEl)) {
				continue;
			}

			ProStecModSat stecmodsat;
			ProStecModSat stecmodref;

			if (!findStecModelSat(isys, prn, stecmod, stecmodsat) ||
				!findStecModelSat(isys, ref, stecmod, stecmodref)) {
				continue;
			}
			
			double f1 = 0.0;




			int tt = 1;

		}
	}
}

void outRovStec(IN Coption& cfg, IN GridInfo& grid, IN SiteAtmos& rovaug, IN ProStecMod& stecmod, IN FileFps& rovfps, IN int type)
{
	char buff[MAXOUTCHARS] = { '/0' }, * p = buff;

	for (const auto& pRov : rovaug) {
		int sid = pRov.first;
		string rov = pRov.second._name;
		if (rovfps.find(rov) == rovfps.end())			 { continue; }
		if (rovfps[rov].find(type) == rovfps[rov].end()) { continue; }
		if (rovfps[rov][type] == NULL)					 { continue; }

		rovStecDiff(cfg, grid, rovaug[sid], stecmod);

		int tt = 1;
	}
}

void createRovFile(IN Coption& cfg, OUT FileFps& fps)
{
	if (0 != _access(cfg._pathou.c_str(), 0)) {
		_mkdir(cfg._pathou.c_str());
	}
	
	for (const auto& pSta : cfg._rov) {
		string rov = pSta._name;
		string rovpath = cfg._pathou + "\\" + rov;

		if (cfg._modeltype & 1) {
			string stecpath = rovpath + "-stec.txt";
			FILE* fp = fopen(stecpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(1, fp);
			}
		}
		if (cfg._modeltype & 2) {
			string ztdpath = rovpath + "-ztd.txt";
			FILE* fp = fopen(ztdpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(2, fp);
			}
		}
		if (cfg._modeltype & 4) {
			string stdpath = rovpath + "-std.txt";
			FILE* fp = fopen(stdpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(3, fp);
			}
		}
		if (cfg._modeltype & 8) {
			string vtecpath = rovpath + "-vtec.txt";
			FILE* fp = fopen(vtecpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(4, fp);
			}
		}
	}
}