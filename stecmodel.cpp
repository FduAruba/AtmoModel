#include "stecmodel.h"

#define		THRESHOLD_MAX_EXPAND_TIME	(120.0)

bool StecGrid::isVaild(IN double lat, IN double lon, IN double dlat, IN double dlon) const
{
	return (_lat > lat - dlat) && (_lat < lat + dlat) && (_lon > lon - dlon) && (_lon < lon + dlon);
}

void StecModSat::reset()
{
	_system = '\0';
	_sat = 0;
	for (int i = 0; i < 4; i++) {
		_coeff[i] = _coeff_rms[i] = 0.0;
	}
	for (int i = 0; i < 2; i++) {
		_QI[i] = 0.0;
	}
	_gridNum = 0;
	_satreslevel = 0;
	_ele = 0.0;
	_stecpergrid.clear();
}

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

bool StecModel::setCurSys(IN int* usesys, IN int symbol)
{
	_cursys = SYS_NONE;

	switch (symbol)
	{
	case 0: {
		_cursys |= usesys[IDX_GPS]  ? SYS_GPS  : SYS_NONE;
		_cursys |= usesys[IDX_GAL]  ? SYS_GAL  : SYS_NONE;
		_cursys |= usesys[IDX_BDS2] ? SYS_BDS2 : SYS_NONE;
		_cursys |= usesys[IDX_BDS3] ? SYS_BDS3 : SYS_NONE;
		break;
	}
	case 1: {
		_cursys |= usesys[IDX_GLO]  ? SYS_GLO  : SYS_NONE;
		break;
	}
	default: { break;}
	}

	return _cursys ? true : false;
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
	/* 计算 ROTI */
	if (!_stecRoti.procRoti(tnow, group)) {
		return false;
	}
	/* 获取当前数据 */
	if (!getAtmo(tnow, group, atmo)) {
		return false;
	}
	/* 统一参考星&计算星间单差STEC */
	if (!uniformRef(atmo)) {
		return false;
	}

	return true;
}

const StecModEpoch* StecModel::StecModInLastEpoch()
{
	return _stecModList.size() > 0 ? &_stecModList.rbegin()->second : NULL;
}

void StecModel::initStecMod(IN AtmoEpoch& atmo)
{
	/* GLONASS不历元初始化 */
	if ((_cursys & SYS_GLO) && _stecModRes._stanum > 0) {
		return;
	}

	/* 初始化建模 */
	_stecModCur._time = atmo._time;
	_stecModCur._modetype = 1;
	for (int isys = 0; isys < NUMSYS; isys++) {
		_stecModCur._satNum[isys] = 0;
		_stecModCur._refsat[isys] = atmo._refSat[isys];
		_stecModCur._stecmodgnss[isys].clear();
	}
	
	/* 初始化残差 */
	_stecModRes._time = atmo._time;
	time2epoch(atmo._time, _stecModRes._ep);
	_stecModRes._stanum = atmo._stanum;
	_stecModRes._staRes.clear();

	return;
}

bool StecModel::checkSatStatus(IN AtmoEpoch& atmo, IN int sys, IN int prn)
{
	int nsta = 0, badquality = 0, unfixsta = 0;
	double rate = 0.0;
	string tstr = strtime(atmo._time, 2);
	char SYS = idx2sys(sys);
	
	/* 1.统计当前卫星roti大且浮点解的站点数 */
	for (const auto& pSta : atmo._staAtmos) {
		const auto pSat = pSta.second._satInfos[sys].find(prn);
		if (pSat != pSta.second._satInfos[sys].end()) {
			nsta++;
			if (_roti[0] > 0.0 && pSat->second._quality > _roti[0]) {
				badquality++;
			}
			if (pSat->second._fixflag != 1) {
				unfixsta++;
			}
		}
	}

	/* TODO:2.roti大的站点数>10% */
	rate = (1.0 * badquality) / (1.0 * atmo._stanum);
	if (rate > 0.1) {
		printf("%s bad roti rate = %.2f%%\n", tstr.c_str(), rate * 100.0);
		return false;
	}

	/* 3.浮点解站点数>80% */
	rate = (1.0 * unfixsta) / (1.0 * atmo._stanum);
	if (_fixsys[sys] && rate > 0.8) {
		printf("%s unfix rate = %.2f%%\n", tstr.c_str(), rate * 100.0);
		return false;
	}

	/* 4.观测站点数>50% */
	rate = (1.0 * nsta) / (1.0 * atmo._stanum);
	if (rate < THRES_USESTA_PCT) {
		//printf("%s %c%02d usable site rate = %6.2f%%\n", tstr.c_str(), SYS, prn, rate * 100.0);
		return false;
	}

	//printf("%s %c%02d usable site rate = %6.2f%%\n", tstr.c_str(), SYS, prn, rate * 100.0);
	return true;
}

bool StecModel::getStecObs(IN AtmoEpoch& atmo, IN int sys, IN int prn, OUT vector<stecOBS>& obss)
{
	int nobs = 0;
	int ref = _stecModCur._refsat[sys];
	
	for (const auto& pSta : atmo._staAtmos) {
		if (!pSta.second._staInfo._supSys[sys]) {
			continue; // 系统支持
		}

		if (nobs > MAX_STEC_OBS) {
			continue; // obs数量上限
		}

		auto pSat = pSta.second._satInfos[sys].find(prn);
		if (pSat == pSta.second._satInfos[sys].end()) {
			continue; // 卫星存在
		}

		double stec = pSat->second._iono;
		double el = pSat->second._azel[1];
		if ((fabs(stec) < DBL_EPSILON && prn != ref) || 
			(pSat->second._quality > _roti[0])) {
			continue; // 卫星stec
		}

		double weight = 1.0;
		if (_fixsys[sys] && pSat->second._fixflag != 1) {
			continue; // 固定解
		}
		if (!_fixsys[sys] && pSat->second._fixflag != 1) {
			weight /= 3.0;
		}

		stecOBS obs;
		obs._id   = pSta.second._staInfo._ID;
		obs._lat  = pSta.second._staInfo._blh[0];
		obs._lon  = pSta.second._staInfo._blh[1];
		obs._stec = stec;
		obs._wgt  = weight;
		obss.push_back(obs);
		nobs++;
	}
	
	return (nobs < (STECNX + 1)) ? false : true;
}

int StecModel::markUnhalthSites(IO vector<stecOBS>& obss)
{
	int nbad[2] = { 0 };
	double vmean = 0.0, vrms = 0.0;
	vector<double> dat;

	do {
		nbad[1] = nbad[0];
		
		for (auto& pOBS : obss) { dat.push_back(pOBS._stec); }
		calcMeanStd(dat, vmean, vrms);

		for (auto it = obss.begin(); it != obss.end();) {
			Pos pos;
			pos._lat = it->_lat;
			pos._lon = it->_lon;
			if ((it->_stec - vmean) > 4 * vrms) {
				_badsta.push_back(pos);
				nbad[0]++;
				it = obss.erase(it);
			}
			else {
				++it;
			}

		}
	} while (nbad[0] > nbad[1]);

	return nbad[0];
}

int StecModel::markUnhalthSitesRes(IN int sys, IN int prn, OUT int& nsta)
{
	int nbad[2] = { 0 };
	double vmean = 0.0, vrms = 0.0;
	vector<double> dat;
	vector<double> cleandat;

	do {
		nsta = 0;
		nbad[1] = nbad[0];

		// 计算所有站的均值&方差
		for (const auto& pSta : _stecModRes._staRes) {
			auto pSat = pSta.second._satInfos[sys].find(prn);
			if (pSat == pSta.second._satInfos[sys].end()) {
				continue;
			}
			dat.push_back(pSat->second._iono);
		}
		calcMeanStd(dat, vmean, vrms);

		// 根据一次方差剔除粗差，重新计算均值&方差
		for (const auto& pSta : _stecModRes._staRes) {
			auto pSat = pSta.second._satInfos[sys].find(prn);
			if (pSat == pSta.second._satInfos[sys].end()) {
				continue;
			}
			if (fabs(pSat->second._iono - vmean) >= 2.0 * vrms) {
				continue;
			}
			cleandat.push_back(pSat->second._iono);
		};
		calcMeanStd(cleandat, vmean, vrms);

		// 根据二次方差剔除粗差站点
		for (auto& pSta : _stecModRes._staRes) {
			Pos pos(pSta.second._staInfo._blh[0], pSta.second._staInfo._blh[1]);
			auto pSat = pSta.second._satInfos[sys].find(prn);
			if (pSat == pSta.second._satInfos[sys].end()) {
				continue;
			}

			double diff = fabs(pSat->second._iono - vmean);
			if (diff > 1.0 && diff > 3.0 * vrms) {
				_badsta.push_back(pos);
				nbad[0]++;
				pSta.second._satInfos[sys].erase(prn);
			}
			else {
				nsta++;
			}
		}
	} while (nbad[0] > nbad[1]);

	return nbad[0];
}

bool StecModel::estCoeff(IN vector<stecOBS>& obss, IN GridInfo& grid, IN int sys, IN int prn, OUT StecModSat& dat)
{
	int nobs = 0;
	double rms = 0.0, rmsmax = 0.0;
	double coeff[STECNX] = { 0.0 }, coeff_rms[STECNX] = { 0.0 };
	VecXd* Y = NULL, * P = NULL, * V = NULL;
	VecXd* X = NULL, * dx = NULL;
	MatXd* H = NULL, * Q = NULL;
	
	dat.reset();

	Y = new VecXd(MAX_STEC_OBS);
	P = new VecXd(MAX_STEC_OBS);
	H = new MatXd(MAX_STEC_OBS, STECNX);
	Y->setZero(); P->setZero(); H->setZero();

	/* 1.矩阵构建 */
	for (auto& pOBS : obss) {
		if (nobs > MAX_STEC_OBS) {
			continue;
		}

		Pos pos(pOBS._lat, pOBS._lon);
		if (find(_badsta.begin(), _badsta.end(), pos) != _badsta.end()) {
			continue;
		}

		double dlat = pos._lat - grid._center[0] * D2R;
		double dlon = pos._lon - grid._center[1] * D2R;

		(*Y)(nobs) = pOBS._stec;
		(*P)(nobs) = pOBS._wgt;
		(*H)(nobs, 0) = 1.0;
		(*H)(nobs, 1) = dlat;
		(*H)(nobs, 2) = dlon;
		nobs++;
	}

	Y->conservativeResize(nobs);
	P->conservativeResize(nobs);
	H->conservativeResize(nobs, STECNX);

	/* 2.最小二乘求解系数 */
	V  = new VecXd(nobs);
	X  = new VecXd(STECNX);
	dx = new VecXd(STECNX);
	Q  = new MatXd(STECNX, STECNX);
	V->setZero(); X->setZero(); dx->setZero(); Q->setZero();

	(*V) = (*Y) - (*H) * (*X);
	for (int i = 0; i < MAX_ITER; i++) {
		// 最小二乘
		(*dx) = wlsq_LU(*H, *V, *P, Q);
		(*X) += (*dx);
		(*V) = (*Y) - (*H) * (*X);

		// 更新权重
		rmsmax = (*V).dot((*V));
		rms = sqrt(rmsmax / (nobs - STECNX));
		for (int j = 0; j < nobs; j++) {
			(*P)(j) *= robust((*V)(j) * sqrt((*P)(j)), rms);
		}

		// 残差判定
		if ((*dx).dot((*dx)) < 1.0E-3) {
			break;
		}
	}

	if (Q->cols() != STECNX || fabs((*X)(0)) > 200.0) {
		Y->resize(0); P->resize(0); V->resize(0);
		X->resize(0); dx->resize(0);
		H->resize(0, 0); Q->resize(0, 0);
		delete Y; delete P; delete V;
		delete X; delete dx; delete H; delete Q;
		return false;
	}

	/* 3.更新rms */
	VecXd V2 = sliceVecByRate(*V);
	rmsmax = V2.dot(V2);
	nobs = (int)V2.rows();
	rms = sqrt(rmsmax / (nobs - STECNX));

	vector<double> arr;
	for (int i = 0; i < V->rows(); i++) {
		arr.push_back((*V)(i));
	}
	stable_sort(arr.begin(), arr.end());

	if (_qi_multi > 0.9) {
		int idx = (int)(arr.size() * _qi_multi);
		if (idx > 0 && idx < arr.size()) {
			rms = arr[idx - 1];
		}
	}
	rms *= _qi_coeff;
	rms += _qi_base;
	rms = rms > 31.51 ? 31.51 : rms;

	for (int i = 0; i < STECNX; i++) {
		coeff[i] = (*X)(i);
		coeff_rms[i] = (*Q)(i, i);
	}

	/* 4.保存建模系数 */
	dat._system = idx2sys(sys);
	dat._sat = prn;
	dat._coeff[0] = coeff[0];
	dat._coeff[1] = coeff[1];
	dat._coeff[2] = coeff[2];
	dat._coeff[3] = 0.0;
	dat._coeff_rms[0] = coeff_rms[0];
	dat._coeff_rms[1] = coeff_rms[1];
	dat._coeff_rms[2] = coeff_rms[2];
	dat._coeff_rms[3] = 0.0;
	dat._QI[0] = rms;
	dat._QI[1] = rms;

	Y->resize(0); P->resize(0); V->resize(0);
	X->resize(0); dx->resize(0); 
	H->resize(0, 0); Q->resize(0, 0);
	delete Y; delete P; delete V;
	delete X; delete dx; delete H; delete Q;

	return true;
}

bool StecModel::oneSatModelEst(IN AtmoEpoch& atmo, IN GridInfo& grid, IN int sys, IN int prn, OUT StecModSat& dat)
{
	string tstr = strtime(atmo._time, 2);
	int SYS = idx2sys(sys);
	bool stat = false;
	vector<stecOBS> obss;

	_badsta.clear();
	
	/* 1.检测当前卫星可否建模 */
	if (!this->checkSatStatus(atmo, sys, prn)) {
		return false;
	}

	/* 2.获取观测值 */
	if (!this->getStecObs(atmo, sys, prn, obss)) {
		return false;
	}

	for (int i = 0; i < 1; i++) {
		/* 3.剔除粗差观测值，当占比>30%，返回false */
		int nbad = markUnhalthSites(obss);
		if (1.0 * nbad / obss.size() > 0.3) {
			printf("%s %c%02d nbad = %2d nall = %2d\n", tstr.c_str(), SYS, prn, nbad, (int)obss.size());
			break;
		}

		/* 4.参数估计 */
		stat = false;
		if (!(stat = estCoeff(obss, grid, sys, prn, dat))) {
			break;
		}
	}
	
	return stat;
}

double StecModel::calcGridTecRes(IN int sys, IN int prn, IN GridEach& grid)
{
	int ref = _stecModCur._refsat[sys];
	double res = ERROR_VALUE;

	StaDistIonArr stalist;
	stalist.reserve(_stecModRes._staRes.size());

	for (const auto& pSta : _stecModRes._staRes) {
		string site = pSta.first;
		auto pSat = pSta.second._satInfos[sys].find(prn);
		if (pSat == pSta.second._satInfos[sys].end()) {
			continue;
		}
		// distance
		double dist = _siteGridDist.getDist(grid._id, site);
		if (pSat->second._fixflag != 1 && sys == IDX_GPS) {
			dist *= 2.0;
		}
		if (dist > 0.5 * CUT_DIST) {
			continue;
		}
		if (dist < 10.0) {
			dist = 10.0;
		}
		// iono
		double ion = pSat->second._iono;
		if (fabs(ion) < DBL_EPSILON && prn != ref) {
			continue;
		}
		stalist.emplace_back(dist, ion);
	}
	if (stalist.size() < 1) { return ERROR_VALUE; }

	sort(stalist.begin(), stalist.end());
	res = modelIDW(stalist, 4, 70000.0, 2);

	return res;
}

double StecModel::calcRovTecRes(IN string site, IN const double* blh, IN GridInfo& grid, IN StecModSat& dat)
{
	double lat = blh[0], lon = blh[1], res = 0.0;
	StaDistIonArr stalist;
	stalist.reserve(dat._gridNum);

	for (int i = 1; i <= dat._gridNum; i++) {
		const auto& onegrid = dat._stecpergrid[i];

		double ion = onegrid._stec;
		if (fabs(ion - ERROR_VALUE) < DBL_EPSILON) {
			continue;
		}

		double dist = _siteGridDist.getDist(i, site);
		if (dist < 10.0) {
			dist = 10.0;
		}

		if (onegrid.isVaild(lat, lon, grid._step[0], grid._step[1])) {
			stalist.emplace_back(dist, ion);
		}
	}

	if (stalist.size() < 1) { return 0.0; }
	stable_sort(stalist.begin(), stalist.end());

	res = modelIDW(stalist, 4, 999000.0, 2);
	res = fabs(res - ERROR_VALUE) < DBL_EPSILON ? 0.0 : res;

	return res;
}

bool StecModel::oneSatStecRes(IN AtmoEpoch& atmo, IN GridInfo& grid, IN int sys, IN int prn, IO StecModSat& dat)
{
	int nres = 0, nbig = 0;
	int nsta = 0, nbad = 0;
	double mean_res = 0.0, mean_ele = 0.0;
	double mean_rot = 0.0, mean_roti = 0.0;
	double rate = 0.0;
	char SYS = idx2sys(sys);

	/* 更新站点残差 */
	for (auto& pSta : atmo._staAtmos) {
		const string site = pSta.first;
		const AtmoSite& atmosta = pSta.second;

		// 跳过问题站
		Pos pos(atmosta._staInfo._blh[0], atmosta._staInfo._blh[1]);
		if (find(_badsta.begin(), _badsta.end(), pos) != _badsta.end()) {
			continue;
		}
		// 跳过问题星
		auto pSat = atmosta._satInfos[sys].find(prn);
		if (pSat == atmosta._satInfos[sys].end()) {
			continue;
		}
		// 跳过浮点解
		if (_fixsys[sys] && pSat->second._fixflag != 1) {
			continue;
		}

		double dlat = pos._lat - grid._center[0] * D2R;
		double dlon = pos._lon - grid._center[1] * D2R;
		double stec = dat._coeff[0] + dat._coeff[1] * dlat + dat._coeff[2] * dlon;
		double diff = pSat->second._iono - stec;

		auto& sta_res = _stecModRes._staRes[site];
		sta_res._staInfo = pSta.second._staInfo;
		sta_res._satInfos[sys][prn]._iono = diff;
		sta_res._satInfos[sys][prn]._fixflag = pSat->second._fixflag;

		double thres = sys == IDX_GLO ? CUT_STEC_RES : 10.0 * CUT_STEC_RES;
		if (fabs(diff) > thres) {
			nbig++;
		}
		mean_res += fabs(diff);
		nres++;
		//printf("%s %c%02d %6.2f\n", site.c_str(), SYS, prn, diff);
		mean_ele  += pSat->second._azel[1];
		mean_rot  += pSat->second._rot;
		mean_roti += pSat->second._quality;
	}
	if (nres <= 0) { return false; }
	//printf("%d\n", nres);
	mean_ele  /= nres;
	mean_rot  /= nres;
	mean_roti /= nres;
	dat._ele = mean_ele;

	/* 粗差筛除 */
	nbad = markUnhalthSitesRes(sys, prn, nsta);
	rate = 1.0 * nbig / nsta;
	if (rate > 0.3) { // 粗差率
		dat._QI[0] = dat._QI[1] = 31.51;
	}
	rate = 1.0 * nres / atmo._staAtmos.size();
	if (rate <= 0.2) { // 站点率
		dat._QI[0] = dat._QI[1] = 31.51; 
	}
	dat._ele = mean_ele;
	if (mean_ele < _minel * D2R) { // 高度角
		//printf("LOW EL: %c%02d %5.2f(deg)\n", SYS, prn, mean_ele * R2D);
		return false;
	}

	/* 保存数据 */
	dat._gridNum = grid._gridNum;
	//int nerr = 0;
	for (int i = 0; i < grid._gridNum; i++) {
		double res = calcGridTecRes(sys, prn, grid._grids[i]);
		if (fabs(res - ERROR_VALUE) < DBL_EPSILON) {
			res = ERROR_VALUE;
			//nerr++;
		}

		auto& stecgrid = dat._stecpergrid[i + 1];
		stecgrid._gridid = i + 1;
		stecgrid._stec = res;
		stecgrid._lat = grid._grids[i]._lat;
		stecgrid._lon = grid._grids[i]._lon;
		stecgrid._rms = 0.0;
		//printf("%c%02d #%02d Res:%8.3f\n", SYS, prn, i + 1, res);
	}
	//printf("%02d\n", nerr);

	return true;
}

bool StecModel::recalculateQI(IN AtmoEpoch& atmo, IN int sys, IN int prn, IN GridInfo& grid, OUT StecModSat& dat)
{
	int nbig = 0, nres = 0;
	double mean_res = 0.0, max_res = 0.0;
	double res = 0.0, rate = 0.0;
	VecXd V = VecXd::Zero(0);
	char SYS = idx2sys(sys);

	vector<double> staRes;
	staRes.reserve(atmo._staAtmos.size());
	
	/* 计算所有格网到所有建模站的残差 */
	for (auto& pSta : atmo._staAtmos) {
		const string site = pSta.first;
		const AtmoSite& atmosta = pSta.second;

		// 跳过问题站
		Pos pos(atmosta._staInfo._blh[0], atmosta._staInfo._blh[1]);
		if (find(_badsta.begin(), _badsta.end(), pos) != _badsta.end()) {
			continue;
		}
		// 跳过问题星
		auto pSat = atmosta._satInfos[sys].find(prn);
		if (pSat == atmosta._satInfos[sys].end()) {
			continue;
		}
		// 跳过浮点解
		if (_fixsys[sys] && pSat->second._fixflag != 1) {
			continue;
		}

		double dlat = pos._lat - grid._center[0] * D2R;
		double dlon = pos._lon - grid._center[1] * D2R;
		double stec = dat._coeff[0] + dat._coeff[1] * dlat + dat._coeff[2] * dlon;
		double res  = calcRovTecRes(site, atmosta._staInfo._blh, grid, dat);
		stec += res;
		double absdiff = fabs(pSat->second._iono - stec);

		double thres = sys == IDX_GLO ? CUT_STEC_RES : 10.0 * CUT_STEC_RES;
		if (absdiff > thres) {
			nbig++;
		}
		nres++;
		mean_res += absdiff;

		staRes.emplace_back(absdiff);
	}
	if (nres <= 0) { return false; }

	/* 更新QI */
	stable_sort(staRes.begin(), staRes.end());
	V = VecXd::Zero(nres);
	for (int i = 0; i < nres; i++) {
		V(i) = staRes[i];
	}
	max_res = V.dot(V);

	rate = 1.0 * nbig / nres;
	if (rate > 0.3) {
		dat._QI[1] = 31.51;
		return true;
	}

	rate = 1.0 * nres / atmo._staAtmos.size();
	if (rate > 0.2) {
		res = sqrt(max_res / (nres - STECNX));
		if (_qi_multi > 0.9) {
			int idx = (int)(nres * _qi_multi);
			if (idx > 0 && idx < nres) {
				res = staRes[idx];
			}
		}
		res *= _qi_coeff;
		res += _qi_base;
		res = res > 31.51 ? 31.51 : res;
	}
	else {
		res = 31.51;
	}
	dat._QI[1] = res;
	return true;
}

bool StecModel::checkSatContinuous(IN Gtime tnow, IN int sys, IN int prn, IN int ref, IN GridInfo& grid, IN StecModSat& dat)
{
	bool isfind = false;
	Gtime tpre = { 0 };
	const StecModSat* pre_sat = NULL;
	const StecModSat* pre_ref = NULL;

	if (dat._stecpergrid.empty()) {
		return false;
	}

	/* 1.遍历卫星建模结果(时间倒序)，查找在规定时间阈值之内是否有当前卫星&参考星，若均存在则返回true */



	/* 2.若无法找到(时间非连续)，遍历格网点 */

	return true;
}

void StecModel::satModEst(IN AtmoEpoch& atmo, IN GridInfo& grid)
{
	/* 1.计算基线距离 */
	_siteGridDist.calcAllDist(atmo, grid);

	/* 2.参数/残差估计 */
	for (int isys = 0; isys < NUMSYS; isys++) {
		if (!(_cursys & _sysidx[isys])) {
			continue;
		}

		for (const auto& pSat : _satList[isys]) {
			StecModSat satdata;

			// 单星建模估计系数
			if (this->oneSatModelEst(atmo, grid, isys, pSat, satdata)) {

				// 单星建模估计格网残差
				if (!this->oneSatStecRes(atmo, grid, isys, pSat, satdata)) {
					continue;
				}
				// 重新计算QI
				if (_useres && satdata._QI[0] < 13.75) {
					this->recalculateQI(atmo, isys, pSat, grid, satdata);
				}

				// 检查建模结果是否连续

				// 插入当前历元系统建模结果

			}
		}
	}

	return;
}