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

double rovRes(IN double lat, IN double lon, IN GridInfo& grid, IN ProStecModSat* satmod, IN ProStecModSat* refmod, OUT int* n, OUT int* stas)
{
	double res = 0.0;
	StaDistIonArr stalist;
	stalist.reserve(grid._gridNum);
	map<double, int> cntgrid;

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
			cntgrid.emplace(dist, satmod->_stecpergrid[i]._nsta);
		}
	}

	if (stalist.size() < 1) { return 0.0; }
	stable_sort(stalist.begin(), stalist.end());

	res = modelIDW(stalist, 4, 999000.0, 2, n);
	res = fabs(res - ERROR_VALUE) < DBL_EPSILON ? 0.0 : res;

	if (n && stas) {
		int i = 0;
		for (auto it = cntgrid.begin(); it != cntgrid.end(); ++it) {
			if (i >= *n) {
				break;
			}
			stas[i++] = it->second;
		}
	}

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
	double THRES_STEC = 0.05;
	string ROV = roviono._name;
	bool stat = false;
	string strt = strtime(stecmod._time, 2);

	ProStecModSat* modsat = new ProStecModSat;
	ProStecModSat* modref = new ProStecModSat;
	
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
			int perGrid[4] = { 0 };
			if (cfg._useres) {
				res = fact * rovRes(latr, lonr, grid, modsat, modref, &ngrid, perGrid);
			}
			double diff0 = fabs(ppp_stec - mod_stec);
			double diff1 = fabs(ppp_stec - mod_stec - res);
			double stec0 = diff0 / fact;
			double stec1 = diff1 / fact;
			
			// DEBUG
			if (diff1 <= diff0 && diff1 < THRES_STEC) {
				cfg._rovstatic[ROV][0]++;
				cfg._ngoodres++;
			}
			else {
				if (diff1 >= THRES_STEC) {
					//double dtec = diff1 / fact;
					/*printf("%s %s %c%02d OUT: res0=%6.3f res1=%6.3f dstec=%6.2f QI=%6.2f el=%5.1f nsta=%2d ngrid=%2d",
						strt.c_str(), ROV.c_str(), SYS, prn, diff0, diff1, dtec, modsat->_QI[1], el, modsat->_nsta, ngrid);
					printf(" %1d %1d %1d %1d\n", perGrid[0], perGrid[1], perGrid[2], perGrid[3]);*/
					cfg._rovstatic[ROV][2]++;
					cfg._noutres++;
				}
				else {
					cfg._rovstatic[ROV][1]++;
					cfg._nbadres++;
				}
			}
			cfg._rovstatic[ROV][3]++;
			cfg._nvali++;

			bool st1 = diff1 <= diff0 ? true : false;
			bool st2 = diff1 >= THRES_STEC && stec1 > modsat->_QI[1] && ngrid >= 1 ? true : false;
			//if (st1 || !st2) { continue; }
			// DEBUG

			OutSatVeri satdat;
			satdat._sys        = SYS;
			satdat._prn        = prn;
			satdat._el         = el;
			satdat._fixflag    = pSat.second._fixflag;
			satdat._ngrid      = ngrid;
			satdat._pergrid[0] = perGrid[0];
			satdat._pergrid[1] = perGrid[1];
			satdat._pergrid[2] = perGrid[2];
			satdat._pergrid[3] = perGrid[3];
			satdat._dstec[2]   = ppp_stec;
			satdat._dstec[0]   = mod_stec;
			satdat._dstec[1]   = mod_stec + res;
			satdat._dDstec[0]  = (ppp_stec - mod_stec) / fact;
			satdat._dDstec[1]  = (ppp_stec - mod_stec - res) / fact;
			satdat._QI[0]      = modsat->_QI[0];
			satdat._QI[1]      = modsat->_QI[1];
			satdat._resflag    = st1 ? 1 : 0;
			satdat._outflag    = st2 ? 0 : 1;

			rovout._satveris[isys].emplace(prn, satdat);
			stat = true;
		}
	}

	if (stat) {
		for (int i = 0; i < NUMSYS; i++) {
			rovout._satNum[i] = (int)rovout._satveris[i].size();
		}
	}

	delete modsat; delete modref;
	return;
}

int printRovEpoch(IN OutRovStec& rovdat, IN FILE* fp)
{
	int nsat = 0, satnum[NUMSYS] = { 0 };
	char buff[MAXOUTCHARS] = { '\0' }, * p = buff;
	string tstr = strtime(rovdat._time, 1);
	
	for (int isys = 0; isys < NUMSYS; isys++) {
		satnum[isys] = rovdat._satNum[isys];
		nsat += rovdat._satNum[isys];
	}
	if (nsat == 0) { return 0; }

	p += sprintf(p, ">%s %4s", tstr.c_str(), rovdat._name.c_str());
	p += sprintf(p, " %2d %2d %2d %2d %2d %2d", nsat, satnum[0], satnum[1], satnum[2], satnum[3], satnum[4]);
	p += sprintf(p, "\n");

	fwrite(buff, (int)(p - buff), sizeof(char), fp);
	return nsat;
}

void printSatStec(IN OutSatVeri& dat, IN FILE* fp)
{
	char buff[MAXOUTCHARS] = { '\0' }, * p = buff;

	p += sprintf(p, "%c%02d %1d", dat._sys, dat._prn, dat._fixflag);
	p += sprintf(p, " %4.1f %7.3f %7.3f %7.3f", dat._el, dat._dstec[0], dat._dstec[1], dat._dstec[2]);
	p += sprintf(p, " %1d %1d %5.2f %5.2f", dat._resflag, dat._outflag, dat._dDstec[0], dat._dDstec[1]);
	p += sprintf(p, " %5.2f %5.2f", dat._QI[0], dat._QI[1]);
	p += sprintf(p, " %1d %1d %1d %1d %1d\n", dat._ngrid, dat._pergrid[0], dat._pergrid[1], dat._pergrid[2], dat._pergrid[3]);

	fwrite(buff, (int)(p - buff), sizeof(char), fp);
}

void printRovStec(IN vector<OutRovStec> rovs, IN FileFps& rovfps, IN int type)
{
	char buff[MAXOUTCHARS] = { '\0' }, * p = buff;

	for (auto& it : rovs) {
		string rov = it._name;
		
		int ns = printRovEpoch(it, rovfps[rov][type]);

		for (int isys = 0; isys < NUMSYS; isys++) {
			for (auto& isat : it._satveris[isys]) {
				printSatStec(isat.second, rovfps[rov][type]);
			}
		}
	}
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

	printRovStec(rovs, rovfps, type);

	return;
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