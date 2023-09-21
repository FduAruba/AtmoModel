#include "stecmodel.h"

#define		THRESHOLD_MAX_EXPAND_TIME	(120.0)

void StecModSat::reset()
{
	_system = '\0';
	_sat = 0;
	for (int i = 0; i < 4; i++) {
		_coff[i] = _coff_res[i] = 0.0;
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

		double weight = 0.0;
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
		
		for (auto& pOBS : obss) {
			dat.push_back(pOBS._stec);
		}
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

bool StecModel::estCoeff(IN vector<stecOBS>& obss, IN GridInfo& grid, IN int sys, IN int prn, OUT StecModSat& dat)
{
	dat.reset();

	int nobs = 0;
	MatXd* Y = new MatXd(MAX_STEC_OBS, 1);
	MatXd* P = new MatXd(MAX_STEC_OBS, MAX_STEC_OBS);
	MatXd* H = new MatXd(MAX_STEC_OBS, STECNX);
	Y->setZero();
	P->setZero();
	H->setZero();

	/* 1.矩阵构建 */
	for (auto& pOBS : obss) {
		if (nobs > MAX_STEC_OBS) {
			continue;
		}

		Pos pos; pos._lat = pOBS._lat; pos._lon = pOBS._lon;
		if (find(_badsta.begin(), _badsta.end(), pos) != _badsta.end()) {
			continue;
		}

		double dlat = pos._lat - grid._center[0] * D2R;
		double dlon = pos._lon - grid._center[1] * D2R;

		(*Y)(nobs, 0) = pOBS._stec;
		(*H)(nobs, 0) = 1.0;
		(*H)(nobs, 1) = dlat;
		(*H)(nobs, 2) = dlon;
		(*P)(nobs, nobs) = pOBS._wgt;
		nobs++;
	}


	/* 2.least square */

	/* 3.更新rms */

	/* 4.保存建模系数 */

	Y->resize(0, 0);
	P->resize(0, 0);
	H->resize(0, 0);
	delete Y;
	delete P;
	delete H;
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
		if (!estCoeff(obss, grid, sys, prn, dat)) {
			break;
		}
	}
	
	return stat;
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
			oneSatModelEst(atmo, grid, isys, pSat, satdata);

			// 单星建模估计格网残差

			// 检查建模结果是否连续

			// 插入最新历元建模结果
		}
	}

	return;
}