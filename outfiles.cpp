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
	/*if (prn == ref) {
		return false;
	}*/
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

bool findStecModelSat(IN int sys, IN int prn, IN ProStecMod& stecmod, OUT ProStecModSat* satmod)
{
	if (stecmod._stecmod[sys].find(prn) != stecmod._stecmod[sys].end()) {
		satmod->deepcopy(stecmod._stecmod[sys][prn]);
		return true;
	}
	return false;
}

double rovRes(IN double lat, IN double lon, IN GridInfo& grid, IN ProStecModSat* satmod, IN ProStecModSat* refmod, OUT int* n)
{
	double res = 0.0;
	StaDistIonArr stalist;
	stalist.reserve(grid._gridNum);

	for (int i = 0; i < grid._gridNum; i++) {
		if (fabs(satmod->_stecpergrid[i]._stec - ERROR_VALUE) < DBL_EPSILON ||
			fabs(refmod->_stecpergrid[i]._stec - ERROR_VALUE) < DBL_EPSILON) {
			continue;
		}
		double ion = satmod->_stecpergrid[i]._stec - refmod->_stecpergrid[i]._stec;
		double latG = satmod->_stecpergrid[i]._lat;
		double lonG = satmod->_stecpergrid[i]._lon;
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

	res = modelIDW(stalist, 4, 999000.0, 2, n);
	res = fabs(res - ERROR_VALUE) < DBL_EPSILON ? 0.0 : res;

	return res;
}

void rovStecDiff(IN Coption& cfg, IN GridInfo& grid, IN FILE* fp, IN SiteAtmo& roviono, IN ProStecMod& stecmod, OUT OutRovStec& rovout)
{
	char buff[MAXOUTCHARS] = { '\0' }, * p = buff;
	double thresEl = cfg._minel * D2R;
	double latr = roviono._blh[0];
	double lonr = roviono._blh[1];
	double lat0 = cfg._center[0] * D2R;
	double lon0 = cfg._center[1] * D2R;
	string ROV = roviono._name;
	bool stat = false;
	string strt = strtime(stecmod._time, 2);

	ProStecModSat* modsat = new ProStecModSat;
	ProStecModSat* modref = new ProStecModSat;
	
	/*int n0 = 0, n1 = 0;
	double dmax = 0.0, dno = 0.0;*/

	rovout._time = stecmod._time;
	rovout._name = ROV;

	for (int isys = 0; isys < NUMSYS; isys++) {
		char SYS = idx2sys(isys);
		int ref = stecmod._satRef[isys];
		if (!checkRef(isys, ref, roviono)) {
			continue;
		}
		const auto& pRef = roviono._satIon[isys][ref];
		rovout._satRef[isys] = ref;

		//dmax = dno = 0.0;
		int ngood = 0, nall = 0;
		for (const auto& pSat : roviono._satIon[isys]) {
			int prn = pSat.first;
			if (!checkSat(isys, ref, prn, roviono, thresEl)) {
				continue;
			}

			modsat->reset(); modref->reset();
			if (!findStecModelSat(isys, prn, stecmod, modsat) ||
				!findStecModelSat(isys, ref, stecmod, modref)) {
				continue;
			}

			double el = pSat.second._azel[1] * R2D;
			double f1 = satfreq(isys, 0);
			double fact = 40.3E16 / pow(f1, 2);
			double ppp_stec = fact * (pSat.second._iono - pRef._iono);
			double mod_stec = fact * ((modsat->_coff[0] - modref->_coff[0]) +
									  (modsat->_coff[1] - modref->_coff[1]) * (latr - lat0) +
									  (modsat->_coff[2] - modref->_coff[2]) * (lonr - lon0));
			int ngrid = 0;
			double res = 0.0;
			if (cfg._useres) {
				res = fact * rovRes(latr, lonr, grid, modsat, modref, &ngrid);
			}
			double diff0 = fabs(ppp_stec - mod_stec);
			double diff1 = fabs(ppp_stec - mod_stec - res);
			if (diff1 <= diff0) {
				ngood++;
			}
			nall++;

			if (diff1 <= diff0 && diff1 < 0.08) {
				cfg._rovstatic[ROV][0]++;
				cfg._ngoodres++;
			}
			else {
				if (diff1 >= 0.08) {
					double dtec = diff1 / fact;
					printf("%s %s %c%02d OUT: res0=%6.3f res1=%6.3f dstec=%6.2f QI=%6.2f el=%5.1f nsta=%2d ngrid=%2d\n",
						strt.c_str(), ROV.c_str(), SYS, prn, diff0, diff1, dtec, modsat->_QI[1], el, modsat->_nsta, ngrid);
					cfg._rovstatic[ROV][2]++;
					cfg._noutl++;
				}
				else {
					cfg._rovstatic[ROV][1]++;
					cfg._nbadres++;
				}
			}
			cfg._rovstatic[ROV][3]++;
			cfg._nvali++;


			/*if (diff1 > diff0 && diff1 > 0.1) {
				printf("%s %s %c%02d outlier: nores=%6.3f withres=%6.3f\n", strt.c_str(), ROV.c_str(), SYS, prn, diff0, diff1);
				cfg._noutl++;
			}
			if (diff1 <= diff0) {
				cfg._ngoodres++;
			}
			else {
				cfg._nbadres++;
			}
			cfg._nvali++;*/
			///* debug ---------------------------------------------------------------------*/ 
			//if (diff1 <= diff0) {
			//	printf("%s %c%02d nores:%5.3f withres:%5.3f [yes]\n", roviono._name.c_str(), idx2sys(isys), prn, diff0, diff1);
			//	n0++;
			//}
			//else {
			//	printf("%s %c%02d nores:%5.3f withres:%5.3f [no]\n", roviono._name.c_str(), idx2sys(isys), prn, diff0, diff1);
			//	if (diff1 > dmax) {
			//		dmax = diff1;
			//		dno = diff0;
			//	}
			//}
			//n1++;
			///* debug ---------------------------------------------------------------------*/

			OutSatVeri satdat;
			satdat._sys = SYS;
			satdat._prn = prn;
			satdat._fixflag = pSat.second._fixflag;
			satdat._resflag = diff1 <= diff0 ? 1 : 0;
			satdat._ngrid = ngrid;
			satdat._el = el;
			satdat._dstec[0] = ppp_stec;
			satdat._dstec[1] = mod_stec;
			satdat._dstec[2] = mod_stec + res;
			satdat._dDstec[0] = diff0;
			satdat._dDstec[1] = diff1;
			rovout._satveris[isys].emplace(prn, satdat);
			stat = true;
		}
		//printf("%c dno=%5.3f dmax=%5.3f\n", idx2sys(isys), dno, dmax);
		rovout._rate[isys] = (1.0 * ngood) / (1.0 * nall);
	}
	//printf("%s nsmall=%3d nall=%3d\n", roviono._name.c_str(), n0, n1);

	if (stat) {
		for (int i = 0; i < NUMSYS; i++) {
			rovout._satNum[i] = (int)rovout._satveris[i].size();
		}
	}

	delete modsat; delete modref;
	return;
}

void outRovStec(IN Coption& cfg, IN GridInfo& grid, IN SiteAtmos& rovaug, IN ProStecMod& stecmod, IN FileFps& rovfps, IN int type)
{
	vector<OutRovStec> rovs;

	for (const auto& pRov : rovaug) {
		int sid = pRov.first;
		string rov = pRov.second._name;
		if (rovfps.find(rov) == rovfps.end())			 { continue; }
		if (rovfps[rov].find(type) == rovfps[rov].end()) { continue; }
		if (rovfps[rov][type] == NULL)					 { continue; }

		if (cfg._rovstatic.find(rov) == cfg._rovstatic.end()) {
			cfg._rovstatic.emplace(rov, vector<double>(4, 0.0));
		}
		
		OutRovStec rovout;
		rovStecDiff(cfg, grid, rovfps[rov][type], rovaug[sid], stecmod, rovout);
		rovs.push_back(rovout);
	}
	int tt = 1;
}

void createRovFile(IN Coption& cfg, OUT FileFps& fps)
{
	if (0 != _access(cfg._pathou.c_str(), 0)) {
		if (!_mkdir(cfg._pathou.c_str())) {
			printf("***ERROR: create dir %s failed, please check it!\n", cfg._pathou.c_str());
			return;
		}
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