#include "readfiles.h"

static int getStations(IN set<string> rovset, OUT Coption& config)
{
	DIR* dir = NULL;
	dirent* file = NULL;
	char tmp[MAXNAME] = { '\0' }, * ext = NULL;
	int i = 0;
	struct stat statbuf;

	if (!(dir = opendir(config._pathin.c_str()))) {
		printf("***ERROR: open dir failed, please check it!\n");
		return 0;
	}

	while ((file = readdir(dir)) != NULL) {
		if (strncmp(file->d_name, ".", 1) == 0)  { continue; } //skip ".",".." for linux
		if (!(ext = strrchr(file->d_name, '.'))) { continue; }
		if (!strstr(ext + 1, "aug"))             { continue; }
		if (strstr(file->d_name, "NRTK"))        { continue; }
		
		Cstation sta;
		sta._path = config._pathin + "\\" + file->d_name;

		stat(sta._path.c_str(), &statbuf);
		size_t filezize = statbuf.st_size;
		if (filezize <= 0) {
			continue; }

		strncpy(tmp, file->d_name, 4); 
		tmp[4] = '\0';
		sta._name.append(tmp);
		sta._ID = ++i;

		if (rovset.find(sta._name) != rovset.end()) {
			config._rov.emplace_back(sta);
		}
		else {
			config._sta.emplace_back(sta);
		}
	}

	return static_cast<int>(config._sta.size());
}

static bool configFile(IN FILE* fp, OUT Coption& config)
{
	char buff[MAXCHARS] = { '\0' };
	char path[MAXCHARS] = { '\0' };
	char* p = buff;
	int ret = 0;
	set<string> rovs;

	while (fgets(buff, sizeof(buff), fp)) {
		if (strstr(buff, "[end]")) {
			break;
		}
		else {
			p = strstr(buff, "=");
			if (p == NULL || (++p) == NULL) { continue; }

			if (strstr(buff, "pathin")) {
				ret = sscanf(p, "%s", path);
				if (ret != 1 || strlen(path) < 4) {
					printf("***ERROR: read config input path fail!\n");
					return false;
				}
				modifyPath(path);
				config._pathin.append(path);
			}
			else if (strstr(buff, "pathou")) {
				ret = sscanf(p, "%s", path);
				if (ret != 1 || strlen(path) < 4) {
					printf("***ERROR: read config output path fail!\n");
					return false;
				}
				modifyPath(path);
				config._pathou.append(path);
			}
			else if (strstr(buff, "rover")) {
				int n = (int)str2num(p + 1, 0, 2);
				for (int i = 0, j = 3; i < n; i++, j += 5) {
					char tmp[MAXNAME] = { '\0' };
					strncpy(tmp, p + j, 4);
					tmp[4] = '\0';

					std::string name;
					name.append(tmp);
					rovs.emplace(name);
				}
			}
		}
	}

	config._pathou = config._pathin + "\\" + config._pathou;

	return getStations(rovs, config) > 0 ? true : false;
}

static bool configTime(IN FILE* fp, OUT Coption& config)
{
	char buff[MAXCHARS] = { '\0' };
	char* p = buff;
	int ret = 0;

	while (fgets(buff, sizeof(buff), fp)) {
		if (strstr(buff, "[end]")) {
			break;
		}
		else {
			p = strstr(buff, "=");
			if (p == NULL || (++p) == NULL) { continue; }

			if (strstr(buff, "ts")) {
				ret = sscanf(p, "%d %d %d %d %d %d", &config._ts[0], &config._ts[1], &config._ts[2], 
													 &config._ts[3], &config._ts[4], &config._ts[5]);
				if (ret != 6) {
					printf("***ERROR: read config beg time fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "te")) {
				ret = sscanf(p, "%d %d %d %d %d %d", &config._te[0], &config._te[1], &config._te[2],
													 &config._te[3], &config._te[4], &config._te[5]);
				if (ret != 6) {
					printf("***ERROR: read config end time fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "ti")) {
				ret = sscanf(p, "%d", &config._ti);
				if (ret != 1) {
					printf("***ERROR: read config time interval fail!\n");
					return false;
				}
			}
		}
	}

	Gtime ts = epoch2time(config._ts);
	Gtime te = epoch2time(config._te);

	return (ts <= te) && (config._ti > 0);
}

static bool configSystem(IN FILE* fp, OUT Coption& config)
{
	char buff[MAXCHARS] = { '\0' };
	char* p = buff;
	int ret = 0;

	while (fgets(buff, sizeof(buff), fp)) {
		if (strstr(buff, "[end]")) {
			break;
		}
		else {
			p = strstr(buff, "=");
			if (p == NULL || (++p) == NULL) { continue; }

			if (strstr(buff, "useres")) {
				ret = sscanf(p, "%d", &config._useres);
				if (ret != 1) {
					printf("***ERROR: read config useres fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "usesys")) {
				ret = sscanf(p, "%d", &config._usesys);
				if (ret != 1) {
					printf("***ERROR: read config usesys fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "fixsys")) {
				ret = sscanf(p, "%d", &config._fixsys);
				if (ret != 1) {
					printf("***ERROR: read config fixsys fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "ressys")) {
				ret = sscanf(p, "%d", &config._ressys);
				if (ret != 1) {
					printf("***ERROR: read config ressys fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "maxsatres")) {
				ret = sscanf(p, "%d", &config._maxsatres);
				if (ret != 1) {
					printf("***ERROR: read config maxsatres fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "refsatsmooth")) {
				ret = sscanf(p, "%d", &config._refsatsmooth);
				if (ret != 1) {
					printf("***ERROR: read config refsatsmooth fail!\n");
					return false;
				}
			}
		}
	}

	return true;
}

static bool configModel(IN FILE* fp, OUT Coption& config)
{
	char buff[MAXCHARS] = { '\0' };
	char* p = buff;
	int ret = 0;

	while (fgets(buff, sizeof(buff), fp)) {
		if (strstr(buff, "[end]")) {
			break;
		}
		else {
			p = strstr(buff, "=");
			if (p == NULL || (++p) == NULL) { continue; }

			if (strstr(buff, "algotype")) {
				ret = sscanf(p, "%d", &config._algotype);
				if (ret != 1) {
					printf("***ERROR: read config algotype fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "modeltype")) {
				ret = sscanf(p, "%d", &config._modeltype);
				if (ret != 1) {
					printf("***ERROR: read config modeltype fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "fittype")) {
				ret = sscanf(p, "%d", &config._fittype);
				if (ret != 1) {
					printf("***ERROR: read config fittype fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "minel")) {
				ret = sscanf(p, "%lf", &config._minel);
				if (ret != 1) {
					printf("***ERROR: read config minel fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "qimulti")) {
				ret = sscanf(p, "%lf", &config._qimulti);
				if (ret != 1) {
					printf("***ERROR: read config qimulti fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "qibase")) {
				ret = sscanf(p, "%lf", &config._qibase);
				if (ret != 1) {
					printf("***ERROR: read config qibase fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "qicoeff")) {
				ret = sscanf(p, "%lf", &config._qicoeff);
				if (ret != 1) {
					printf("***ERROR: read config qicoeff fail!\n");
					return false;
				}
			}
		}
	}

	return true;
}

static bool configRegion(IN FILE* fp, OUT Coption& config)
{
	char buff[MAXCHARS] = { '\0' };
	char* p = buff;
	int ret = 0;

	while (fgets(buff, sizeof(buff), fp)) {
		if (strstr(buff, "[end]")) {
			break;
		}
		else {
			p = strstr(buff, "=");
			if (p == NULL || (++p) == NULL) { continue; }

			if (strstr(buff, "latgrid")) {
				ret = sscanf(p, "%lf %lf", &config._latgrid[0], &config._latgrid[1]);
				if (ret != 2 || config._latgrid[0] >= config._latgrid[1]) {
					printf("***ERROR: read config latgrid fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "longrid")) {
				ret = sscanf(p, "%lf %lf", &config._longrid[0], &config._longrid[1]);
				if (ret != 2 || config._longrid[0] >= config._longrid[1]) {
					printf("***ERROR: read config longrid fail!\n");
					return false;
				}
			}
			else if (strstr(buff, "step")) {
				ret = sscanf(p, "%lf %lf", &config._step[0], &config._step[1]);
				if (ret != 2 || config._step[0] <= 0.0 || config._step[1] <= 0.0) {
					printf("***ERROR: read config lat/lon step fail!\n");
					return false;
				}
			}
		}
	}

	config._latcell[0] = config._latgrid[0] - config._step[0] / 2.0;
	config._latcell[1] = config._latgrid[1] + config._step[0] / 2.0;
	config._loncell[0] = config._longrid[0] - config._step[1] / 2.0;
	config._loncell[1] = config._longrid[1] + config._step[1] / 2.0;
	config._length[0]  = config._latcell[1] - config._latcell[0];
	config._length[1]  = config._loncell[1] - config._loncell[0];
	config._center[0]  = (config._latgrid[0] + config._latgrid[1]) / 2.0;
	config._center[1]  = (config._longrid[0] + config._longrid[1]) / 2.0;

	return true;
}

extern bool readConfigFile(IN const char* fname, OUT Coption& config)
{
	FILE* fp = NULL;
	char buff[MAXCHARS] = { '\0' };
	char comm[MAXCOMMENT] = "ERROR_EMPTY";
	char* p = buff;
	int type = 0;

	if (fname == NULL || strlen(fname) == 0) {
		return false;
	}
	if ((fp = fopen(fname, "r")) == NULL) {
		return false;
	}

	while (fgets(buff, sizeof(buff), fp)) {
		if (buff[0] == '#' || strlen(buff) <= 1) {
			continue;
		}
		else if (strstr(buff, "[beg]")) {
			p = strstr(buff, "]");
			if (p == NULL || (++p) == NULL) { continue; }

			if (2 == sscanf(p, "%d %s", &type, comm)) {
				printf("read config: [%d]%s...\n", type, comm);
				switch (type)
				{
				case 1: { // file
					if (!configFile(fp, config))   { return false; }
					break;
				}
				case 2: { // time
					if (!configTime(fp, config))   { return false; }
					break;
				}
				case 3: { // system
					if (!configSystem(fp, config)) { return false; }
					break;
				}
				case 4: { // model
					if (!configModel(fp, config))  { return false; }
					break;
				}
				case 5: { // region
					if (!configRegion(fp, config)) { return false; }
					break;
				}
				default: {
					printf("***WARNNING: invalid config type [%d]\n", type);
					break;
				}
				}
			}
			else {
				printf("***ERROR: decode config head [%d]%s fail!\n", type, comm);
				return false;
			}
		}
	}

	return true;
}

static int decodeEpoch(IN Gtime tnow, IN char* line, OUT SiteAtmo& sta, OUT Gtime& tnext)
{
	int ret = 0, ns = 0, id = 0;
	int ep[6] = { 0 };
	char st[5] = { '\0' };
	double blh[3] = { 0.0 };
	double zhd = 0.0, zwd = 0.0, std = 0.0;

	ret = sscanf(line, "* %d %d %d %d %d %d.00 %d %s %d %lf %lf %lf %lf %lf %lf",
				 ep, ep + 1, ep + 2, ep + 3, ep + 4, ep + 5, &ns, st, &id,
				 blh, blh + 1, blh + 2, &zhd, &zwd, &std);
	if (ret != 15) { return false; }

	blh[0] *= D2R; blh[1] *= D2R;

	Gtime t = epoch2time(ep);
	if (t > tnow) {
		tnext = t;
		return 0;
	}
	else if (t == tnow) {
		tnext = t;
		for (int i = 0; i < 3; i++) {
			sta._blh[i] = blh[i];
		}
		pos2ecef(blh, sta._xyz);
		sta._zhd = zhd;
		sta._zwd = zwd;
		sta._std_zwd = std;
	}
	else {
		return 0;
	}
	
	return ns;
}

static int decodeSatData(IN Gtime tnow, IN FILE* fp, IN int ns, OUT SatAtmos* satinfo, OUT int& nbad)
{
	char buff[MAXCHARS] = { '\0' };
	int ret = 0, fix = -1, prn = 0, cnt = 0;
	char sid[4] = { '\0' };
	double tec = 0.0, ele = 0.0, amb[2] = { 0.0 };

	for (int i = 0; i < ns; i++) {
		if (fgets(buff, sizeof(buff), fp) && checksys(buff[0])) {
			ret = sscanf(buff, "%s %d %lf %lf %lf %lf", sid, &fix, &ele, &tec, amb, amb + 1);
			if (ret == 6) {
				int sat = satid2no(sid);
				int sys = satsys(sat, &prn);

				if (!sat || !sys) {
					printf("***ERROR: sat=%d sys=%d\n", sat, sys);
					continue;
				}

				// skip BDS3 IGSO
				if (sys == SYS_BDS3 && (prn == 38 || prn == 39 || prn == 40)) {
					cnt++;
					continue;
				}

				SatAtmo info;
				info._sys     = sysstr(sys);
				info._prn     = prn;
				info._fixflag = fix;
				info._azel[1] = ele * D2R;
				info._iono	  = tec;
				info._nflo[0] = amb[0];
				info._nflo[1] = amb[1];

				switch (sys)
				{
				case SYS_GPS:  { satinfo[0].emplace(prn, info); break; }
				case SYS_GLO:  { satinfo[1].emplace(prn, info); break; }
				case SYS_GAL:  { satinfo[2].emplace(prn, info); break; }
				case SYS_BDS2: { satinfo[3].emplace(prn, info); break; }
				case SYS_BDS3: { satinfo[4].emplace(prn, info); break; }
				default: { break; }
				}

				cnt++;
			}
		}
		else {
			nbad++;
		}
	}

	return cnt;
}

static int decodeData(IN Gtime tnow, OUT vector<Cstation>& sta, OUT SiteAtmos& stas)
{
	char buff[MAXCHARS] = { '\0' };
	
	for (auto& pSta : sta) { // loop stations
		if (pSta._fp == NULL) {
			pSta._fp = fopen(pSta._path.c_str(), "r");
		}
		if (pSta._fp == NULL) { continue; }

		Gtime t = { 0 };
		SiteAtmo station;
		station._name = pSta._name;
		station._ID = pSta._ID;

		/* decode epoch */
		int ns = 0;
		if (!pSta._line.empty()) {
			ns = decodeEpoch(tnow, (char*)pSta._line.c_str(), station, t);
		}
		else {
			while (fgets(buff, sizeof(buff), pSta._fp)) {
				if (buff[0] == '*') {
					ns = decodeEpoch(tnow, buff, station, t);

					if (ns == 0 && t > tnow) {
						//int dep = (int)(timediff(t, tnow) / config._ti);
						//printf("%s: %s [jump] dt=%d \n", strtime(tnow, 2).c_str(), pSta._name.c_str(), dep);
						pSta._line = buff;
						break;
					}
					else if (ns > 0) {
						break;
					}
				}
			}
		}

		/* decode sat data */
		if (ns > 0 && tnow == t) {
			pSta._line.clear();
			pSta._nepo++;
			int nb = 0;
			int nc = decodeSatData(tnow, pSta._fp, ns, station._satIon, nb);
			if (nb && nc == 0) {
				pSta._nbad++;
				//printf("%s: %s [lack] ns=%2d nb=%2d\n", strtime(tnow, 2).c_str(), pSta._name.c_str(), ns, nb);
			}
		}

		/* insert site data */
		if (station._satIon->size() > 0) {
			stas.emplace(station._ID, station);
		}
	}

	return (int)stas.size();
}

extern bool readAugmentData(IN Gtime tnow, IN Coption& config, OUT SiteAtmos& stas, OUT SiteAtmos& rovs)
{
	int nall = (int)config._sta.size();

	int nsta = decodeData(tnow, config._sta, stas);
	int nrov = decodeData(tnow, config._rov, rovs);

	double thres = 50.0;
	double rate = (double)nsta / (double)nall * 100.0;

	if (rate < thres) {
		//printf("%s nsta=%3d nrov=%3d %5.1f%%\n", strtime(tnow, 2).c_str(), nsta, nrov, rate);
		config._nlack++;
	}

	return nsta > 0 ? true : false;
}

extern void movOption(IN Coption& config, OUT ProOption& opt)
{
	opt._maxroti[0] = 0.5;

	opt._usesys[0] = config._usesys & SYS_GPS  ? 1 : 0;
	opt._usesys[1] = config._usesys & SYS_GLO  ? 1 : 0;
	opt._usesys[2] = config._usesys & SYS_GAL  ? 1 : 0;
	opt._usesys[3] = config._usesys & SYS_BDS2 ? 1 : 0;
	opt._usesys[4] = config._usesys & SYS_BDS3 ? 1 : 0;

	opt._fixsys[0] = config._fixsys & SYS_GPS  ? 1 : 0;
	opt._fixsys[1] = config._fixsys & SYS_GLO  ? 1 : 0;
	opt._fixsys[2] = config._fixsys & SYS_GAL  ? 1 : 0;
	opt._fixsys[3] = config._fixsys & SYS_BDS2 ? 1 : 0;
	opt._fixsys[4] = config._fixsys & SYS_BDS3 ? 1 : 0;

	opt._ressys[0] = config._ressys & SYS_GPS  ? 1 : 0;
	opt._ressys[1] = config._ressys & SYS_GLO  ? 1 : 0;
	opt._ressys[2] = config._ressys & SYS_GAL  ? 1 : 0;
	opt._ressys[3] = config._ressys & SYS_BDS2 ? 1 : 0;
	opt._ressys[4] = config._ressys & SYS_BDS3 ? 1 : 0;

	for (int i = 0; i < NUMSYS; i++) {
		if (opt._usesys[i] == 0) {
			opt._fixsys[i] = opt._ressys[i] = 0;
		}
	}

	opt._ti        = config._ti;
	opt._minel     = config._minel;
	opt._maxsatres = config._maxsatres;
	opt._qimulti   = config._qimulti;
	opt._qibase    = config._qibase;
	opt._qicoeff   = config._qicoeff;
	opt._refsatsmooth = config._refsatsmooth;
	opt._algotype  = config._algotype;
	opt._fittype   = config._fittype;
}

extern bool movGrids(IN Coption& config, OUT GridInfo& grid)
{
	for (int i = 0; i < 2; i++) {
		grid._latgrid[i] = config._latgrid[i];
		grid._latcell[i] = config._latcell[i];
		grid._longrid[i] = config._longrid[i];
		grid._loncell[i] = config._loncell[i];
		grid._step[i]    = config._step[i];
		grid._center[i]  = config._center[i];
		grid._length[i]  = config._length[i];
	}

	int idx = 0;
	double dlat = grid._step[0];
	double dlon = grid._step[1];
	for (double ilat = grid._latgrid[0]; ilat <= grid._latgrid[1]; ilat += dlat) {
		for (double ilon = grid._longrid[0]; ilon <= grid._longrid[1]; ilon += dlon) {
			grid._grids[idx]._id = idx + 1;
			grid._grids[idx]._lat = ilat * D2R;
			grid._grids[idx]._lon = ilon * D2R;

			if (++idx > MAX_GRID) {
				printf("***WARNNING: too many grids %d>%d\n", idx, MAX_GRID);
				printf("   latmin=%6.1f latmax=%6.1f dlat=%6.1f\n", grid._latgrid[0], grid._latgrid[1], dlat);
				printf("   lonmin=%6.1f lonmax=%6.1f dlon=%6.1f\n", grid._longrid[0], grid._longrid[1], dlon);
				return false;
			}
		}
	}
	grid._gridNum = idx;

	return true;
}

extern bool movAtmos(IN Gtime tnow, SiteAtmos& stas, OUT AtmoInfo& stecinf)
{
	for (auto& pSta : stas) {
		int ID = pSta.first;
		auto oneatmoinf = &stecinf._sitesols[ID];
		
		oneatmoinf->_time = tnow;
		oneatmoinf->_name = pSta.second._name;
		oneatmoinf->_ID   = pSta.second._ID;
		for (int i = 0; i < 3; i++) {
			oneatmoinf->_xyz[i] = pSta.second._xyz[i];
			oneatmoinf->_blh[i] = pSta.second._blh[i];
		}
		oneatmoinf->_zhd = pSta.second._zhd;
		oneatmoinf->_zwd = pSta.second._zwd;
		oneatmoinf->_std_zwd = pSta.second._std_zwd;

		for (int isys = 0; isys < NUMSYS; isys++) {
			if (stas[ID]._satIon[isys].size() <= 0) {
				continue;
			}

			for (auto& pSat : stas[ID]._satIon[isys]) {
				int prn = pSat.first;
				auto onesatinf = &oneatmoinf->_satsols[isys][prn];

				onesatinf->_sys     = pSat.second._sys;
				onesatinf->_prn     = pSat.second._prn;
				onesatinf->_fixflag = pSat.second._fixflag;
				onesatinf->_iono    = pSat.second._iono;
				onesatinf->_covIono = pSat.second._covIono;
				onesatinf->_trop    = pSat.second._trop;
				onesatinf->_covTrop = pSat.second._covTrop;
				onesatinf->_quality = pSat.second._quality;
				onesatinf->_rot     = pSat.second._rot;
				for (int i = 0; i < 2; i++) {
					onesatinf->_azel[i] = pSat.second._azel[i];
				}
				for (int i = 0; i < 3; i++) {
					onesatinf->_xyz[i] = pSat.second._xyz[i];
				}
			}
			oneatmoinf->_satNum[isys] = (int)oneatmoinf->_satsols[isys].size();
		}
	}
	stecinf._time = tnow;
	stecinf._sitenum = (int)stecinf._sitesols.size();

	return stecinf._sitenum > 0 ? true: false;
}

extern void freeFps(IO Coption& cfg, IO FileFps& fps)
{
	for (auto& pSta : cfg._sta) {
		if (pSta._fp != NULL) {
			fclose(pSta._fp); pSta._fp = NULL;
		}
	}

	for (auto& pSta : cfg._rov) {
		if (pSta._fp != NULL) {
			fclose(pSta._fp); pSta._fp = NULL;
		}
	}

	for (auto& pSta : fps) {
		for (auto& iFp : pSta.second) {
			if (iFp.second != NULL) {
				fclose(iFp.second); iFp.second = NULL;
			}
		}
	}
}