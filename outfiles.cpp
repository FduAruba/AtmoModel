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

double rovRes(IN double lat, IN double lon, IN GridInfo& grid, IN ProStecModSat& satmod, IN ProStecModSat& refmod)
{
	double res = 0.0;
	StaDistIonArr stalist;
	stalist.reserve(grid._gridNum);

	for (int i = 0; i < grid._gridNum; i++) {
		if (fabs(satmod._stecpergrid[i]._stec - ERROR_VALUE) < DBL_EPSILON ||
			fabs(refmod._stecpergrid[i]._stec - ERROR_VALUE) < DBL_EPSILON) {
			continue;
		}
		double ion = satmod._stecpergrid[i]._stec - refmod._stecpergrid[i]._stec;
		double latG = satmod._stecpergrid[i]._lat;
		double lonG = satmod._stecpergrid[i]._lon;
		double dist = sphereDist(latG, lonG, lat, lon);
		if (dist < 10.0) {
			dist = 10.0;
		}
		if (gridVaild(latG, lonG, lat, lon, grid._step[0] * D2R, grid._step[1] * D2R)) {
			stalist.emplace_back(dist, ion);
		}
	}

	if (stalist.size() < 1) { return 0.0; }
	stable_sort(stalist.begin(), stalist.end());

	res = modelIDW(stalist, 4, 999000.0, 2);
	res = fabs(res - ERROR_VALUE) < DBL_EPSILON ? 0.0 : res;

	return res;
}

void rovStecDiff(IN Coption& cfg, IN GridInfo& grid, IN SiteAtmo& roviono, IN ProStecMod& stecmod)
{
	double thresEl = cfg._minel * D2R;
	double latr = roviono._blh[0];
	double lonr = roviono._blh[1];
	double lat0 = cfg._center[0] * D2R;
	double lon0 = cfg._center[1] * D2R;
	
	int n0 = 0, n1 = 0;
	double dmax = 0.0, dno = 0.0;

	for (int isys = 0; isys < NUMSYS; isys++) {
		int ref = stecmod._satRef[isys];
		if (!checkRef(isys, ref, roviono)) {
			continue;
		}
		const auto& pRef = roviono._satIon[isys][ref];

		dmax = dno = 0.0;
		for (const auto& pSat : roviono._satIon[isys]) {
			int prn = pSat.first;
			if (!checkSat(isys, ref, prn, roviono, thresEl)) {
				continue;
			}

			ProStecModSat modsat;
			ProStecModSat modref;

			if (!findStecModelSat(isys, prn, stecmod, modsat) ||
				!findStecModelSat(isys, ref, stecmod, modref)) {
				continue;
			}
			
			double f1 = satfreq(isys, 0);
			double fact = 40.3E16 / pow(f1, 2);

			double ppp_stec = fact * (pSat.second._iono - pRef._iono);
			double mod_stec = fact * ((modsat._coff[0] - modref._coff[0]) +
									  (modsat._coff[1] - modref._coff[1]) * (latr - lat0) +
									  (modsat._coff[2] - modref._coff[2]) * (lonr - lon0));

			double res = 0.0;
			if (cfg._useres) {
				res = rovRes(latr, lonr, grid, modsat, modref);
			}

			// debug
			double diff0 = fabs(ppp_stec - mod_stec);
			double diff1 = fabs(ppp_stec - mod_stec - fact * res);
			if (diff1 <= diff0) {
				//printf("%s %c%02d nores:%5.3f withres:%5.3f [yes]\n", roviono._name.c_str(), idx2sys(isys), prn, diff0, diff1);
				n0++;
			}
			else {
				//printf("%s %c%02d nores:%5.3f withres:%5.3f [no]\n", roviono._name.c_str(), idx2sys(isys), prn, diff0, diff1);
				if (diff1 > dmax) {
					dmax = diff1;
					dno = diff0;
				}
			}
			n1++;
			// debug

			mod_stec += fact * res;
		}
		printf("%c dno=%5.3f dmax=%5.3f\n", idx2sys(isys), dno, dmax);
	}
	printf("%s nsmall=%3d nall=%3d\n", roviono._name.c_str(), n0, n1);

	return;
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