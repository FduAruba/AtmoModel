#include "localmodel.h"

void LocalAtmoModel::setBasic(IN int res)
{
	this->_stecPro.setBasicOption(this->_proOption, res);
	this->_ztdPro.setBasicOption(this->_proOption, res);
}

void LocalAtmoModel::setReigon(IN GridInfo& grid)
{
	_gridinfo.deepcopy(grid);
}

void LocalAtmoModel::setOption(IN ProOption& opt)
{
	_proOption.deepcopy(opt);
}

bool LocalAtmoModel::setRefSites(IN SiteAtmos& stas)
{
	_stanumGEC = _stanumR = 0;
	for (auto pSta : stas) {
		int id = pSta.first;

		if (pSta.second._blh[0] * R2D < _gridinfo._latcell[0] ||
			pSta.second._blh[0] * R2D > _gridinfo._latcell[1] ||
			pSta.second._blh[1] * R2D < _gridinfo._loncell[0] ||
			pSta.second._blh[1] * R2D > _gridinfo._loncell[1]) {
			//printf("%s %8.2f%8.2f\n", pSta.second._name.c_str(), pSta.second._blh[0] * R2D, pSta.second._blh[1] * R2D);
			continue;
		}

		if (_allsites.find(id) != _allsites.end()) {
			for (int i = 0; i < NUMSYS; i++) {
				if (pSta.second._satIon[i].size() > 0) {
					_allsites[id]._satNum[i] = (int)pSta.second._satIon[i].size();
					_allsites[id]._supSys[i] = 1;
				}
				else {
					_allsites[id]._satNum[i] = 0;
					_allsites[id]._supSys[i] = 0;
				}
			}
			if (pSta.second._satIon[0].size() > 0 ||
				pSta.second._satIon[2].size() > 0 ||
				pSta.second._satIon[3].size() > 0 ||
				pSta.second._satIon[4].size() > 0) {
				_stanumGEC++;
			}
			if (pSta.second._satIon[1].size() > 0) {
				_stanumR++;
			}
		}
		else {
			SiteInfo onesite;
			onesite._name = pSta.second._name;
			onesite._ID = id;
			for (int i = 0; i < 3; i++) {
				onesite._xyz[i] = pSta.second._xyz[i];
				onesite._blh[i] = pSta.second._blh[i];
			}
			for (int i = 0; i < NUMSYS; i++) {
				if (pSta.second._satIon[i].size() > 0) {
					onesite._satNum[i] = (int)pSta.second._satIon[i].size();
					onesite._supSys[i] = 1;
				}
				else {
					onesite._satNum[i] = 0;
					onesite._supSys[i] = 0;
				}
			}
			if (pSta.second._satIon[0].size() > 0 ||
				pSta.second._satIon[2].size() > 0 ||
				pSta.second._satIon[3].size() > 0 ||
				pSta.second._satIon[4].size() > 0) {
				_stanumGEC++;
			}
			if (pSta.second._satIon[1].size() > 0) {
				_stanumR++;
			}
			_allsites.emplace(id, onesite);
		}
	}
	
	return _allsites.size() > 0 ? true : false;
}

bool LocalAtmoModel::getAtmoEpoch(IN Gtime tnow, IN int symbol, OUT AtmoEpochs::iterator& it)
{
	AtmoEpochs* group = NULL;
	switch (symbol)
	{
	case 0: { group = &_stecGroupGEC; break;}
	case 1: { group = &_stecGroupR;   break;}
	default: { return false;}
	}

	if ((it = group->find(tnow)) == group->end()) {
		return false;
	}

	return true;
}

int LocalAtmoModel::stationSys(IN SiteSol& sol)
{
	int ret = 0;
	if (sol._satsols[0].size() > 0 ||
		sol._satsols[2].size() > 0 ||
		sol._satsols[3].size() > 0 ||
		sol._satsols[4].size() > 0) {
		ret += 1;
	}
	if (sol._satsols[1].size() > 0) {
		ret += 2;
	}
	return ret;
}

void LocalAtmoModel::cntFixSat(IN SatSols& sol, OUT map<int, int>& track)
{
	for (auto pSat : sol) {
		int prn = pSat.first;
		if (pSat.second._fixflag == 1) {
			track[prn]++;
		}
	}
}

int LocalAtmoModel::fixDecision(IN AtmoInfo& stecinf, OUT map<int, int>* track)
{
	int nfixGEC = 0, nfixR = 0, nstaGEC = 0, nstaR = 0, fixflag = 0;

	for (auto pSta : stecinf._sitesols) {

		int ret = stationSys(pSta.second);

		if (ret == 1 || ret == 3) {
			nstaGEC++;
			if (pSta.second._satsols[IDX_GPS].size() > 0) {
				cntFixSat(pSta.second._satsols[0], track[0]);
			}
			if (pSta.second._satsols[IDX_GAL].size() > 0) {
				cntFixSat(pSta.second._satsols[2], track[2]);
			}
			if (pSta.second._satsols[IDX_BDS2].size() > 0) {
				cntFixSat(pSta.second._satsols[3], track[3]);
			}
			if (pSta.second._satsols[IDX_BDS3].size() > 0) {
				cntFixSat(pSta.second._satsols[4], track[4]);
			}

			if (track[0].size() > 0 || track[2].size() > 0 || 
				track[3].size() > 0 || track[4].size() > 0) {
				nfixGEC++;
			}
		}

		if (ret == 2 || ret == 3) {
			nstaR++;
			if (pSta.second._satsols[IDX_GLO].size() > 0) {
				cntFixSat(pSta.second._satsols[1], track[1]);
			}

			if (track[1].size() > 0) {
				nfixR++;
			}
		}
	}

	if ((1.0 * nfixGEC / nstaGEC) > THRES_FIXSOL_PCT) {
		fixflag = 1;
	}

	return fixflag;
}

void LocalAtmoModel::sortSiteNumber(IN int isys, IN map<int, EleStation>* numstas, OUT map<int, set<EleStation>>& vprn)
{
	for (auto& pSat : numstas[isys]) {
		if (isys == IDX_BDS2 && pSat.first < 6) {
			continue;
		}

		int prn = pSat.first;
		pSat.second._el  /= (double)pSat.second._num;
		pSat.second._rot /= (double)pSat.second._num;

		if (pSat.second._el < _proOption._minel * D2R) {
			//printf("skip %d %2d el=%.2f(deg)\n", isys, prn, pSat.second._el * R2D);
			continue;
		}

		EleStation tmp = pSat.second;
		if (vprn.find(pSat.second._num) == vprn.end()) {
			set<EleStation> vSat;
			vSat.emplace(tmp);
			vprn.emplace(pSat.second._num, vSat);
		}
		else {
			vprn[pSat.second._num].emplace(tmp);
		}
	}
}

void LocalAtmoModel::findRefSat(IN int isys, IN map<int, EleStation>* numstas, IN map<int, set<EleStation>>& vprn, OUT int* refsat)
{
	int stanum = isys == IDX_GLO ? _stanumR : _stanumGEC;
	bool isfind = false;
	double pctold = 0.0, pctnew = 0.0;
	refsat[isys] = _stecPro._stecModCur._refsat[isys];
	map<int, EleStation>::iterator pold, pnew;
	const StecModEpoch* modLastEpoch = _stecPro.StecModInLastEpoch();

	pold = numstas[isys].find(refsat[isys]);
	pnew = numstas[isys].end();
	if (pold == numstas[isys].end() || pold->second._el < 30.0 * D2R) {
		refsat[isys] = 0;
		pold = numstas[isys].end();
	}

	for (auto pNum = vprn.rbegin(); pNum != vprn.rend(); ++pNum) {
		// 1>旧参考星的共视测站数>当前卫星且高度角>40
		if (pold != numstas[isys].end()) {
			pctold = (1.0 * pold->second._num) / stanum;
			pctnew = (1.0 * pNum->first) / stanum;
			if (pctold - pctnew > 0.1 && pold->second._el > 40.0 * D2R) {
				break;
			}
		}
		// 2>找出高度角>40，且当前共视站最多的卫星
		for (auto pSat = pNum->second.rbegin(); pSat != pNum->second.rend(); ++pSat) {
			int prn = pSat->_prn;

			if (pSat->_el > 40.0 * D2R) {
				pold = numstas[isys].find(refsat[isys]);
				pnew = numstas[isys].find(prn);

				if (pold == numstas[isys].end()) {
					refsat[isys] = prn;
					pold = numstas[isys].find(prn);
					isfind = true;
				}
				else if (pold != numstas[isys].end() && pnew != numstas[isys].end()) {
					pctold = (1.0 * pold->second._num) / stanum;
					pctnew = (1.0 * pnew->second._num) / stanum;
					//TODO: rot判定
					if (pctnew - pctold > 0.1 || pold->second._el < 40.0 * D2R) {
						if (modLastEpoch != NULL && 
							modLastEpoch->_stecmodgnss[isys].find(prn) != modLastEpoch->_stecmodgnss[isys].end()) {
							refsat[isys] = prn;
							isfind = true;
						}
					}
				}
			}
		}
		// 3>若所有卫星均高度角<40，直接取最多共视站的卫星作为参考星
		if (refsat[isys] == 0 || pctnew - pctold > 0.1) {
			refsat[isys] = vprn.rbegin()->second.rbegin()->_prn;
			pold = numstas[isys].find(refsat[isys]);
		}
	}
	// 4>共视站比差距没有>10%，直接取最多共视站的卫星作为参考星
	if (!isfind && refsat[isys] == 0) {
		refsat[isys] = vprn.rbegin()->second.rbegin()->_prn;
	}
}

void LocalAtmoModel::uniformDatum(IN Gtime tnow, IN int* refsat, IN int symbol)
{
	AtmoEpochs::iterator pAtmo;
	if (!getAtmoEpoch(tnow, symbol, pAtmo)) {
		return;
	}
	
	for (auto& pSta : pAtmo->second._staAtmos) {
		string ista = pSta.first;

		for (int isys = 0; isys < NUMSYS; isys++) {
			auto pRef = pSta.second._satInfos[isys].find(refsat[isys]);

			if (refsat[isys] == 0 || pSta.second._satInfos[isys].size() == 0) {
				pSta.second._staInfo._satNum[isys] = 0;
				pSta.second._satInfos[isys].clear();
			}
			else if (pRef == pSta.second._satInfos[isys].end()) {
				//printf("%s: %1d %2d erase %s\n", strtime(tnow, 2).c_str(), isys, refsat[isys], ista.c_str());
				pSta.second._staInfo._satNum[isys] = 0;
				pSta.second._satInfos[isys].clear();
			}
			else {
				_stecPro._stecModCur._refsat[isys] = refsat[isys];
			}
		}
	}

	for (int isys = 0; isys < NUMSYS; isys++) {
		pAtmo->second._refSat[isys] = refsat[isys];
	}
}

bool LocalAtmoModel::selectOneRefSat(IN map<int, EleStation>* numstas, OUT int* refsat)
{
	bool stat = false;
	map<int, set<EleStation>> vPrn;
	
	for (int isys = 0; isys < NUMSYS; isys++) {
		if (_proOption._usesys[isys] == 0 || numstas[isys].empty()) {
			continue;
		}
		// 整理卫星对应的站点数
		vPrn.clear();
		sortSiteNumber(isys, numstas, vPrn);

		if (vPrn.size() <= 0) {
			//printf("no valid sat %d %d\n", isys, (int)numstas[isys].size());
			continue;
		}

		// 筛选参考星
		findRefSat(isys, numstas, vPrn, refsat);

		stat |= refsat[isys] > 0 ? true : false;
	}

	return stat;
}

void LocalAtmoModel::inputAtmoSat(IN int isys, IN SatSols* src, OUT SatInfos* dst, OUT map<int, EleStation>* numstas)
{
	for (auto& pSat : src[isys]) {
		if (_stecPro._fixsys[isys] == 1 && pSat.second._fixflag != 1) {
			continue;
		}
		
		int prn = pSat.first;
		SatInfo info(pSat.second);

		if (isnan(info._iono) || isinf(info._iono)) {
			continue;
		}
		//TODO: 判定roti是否大于阈值0.5，中海达暂时无法输出

		// 插入观测数据
		dst[isys].emplace(prn, info);
		// 插入卫星列表
		_stecPro._satList[isys].emplace(prn);
		// 统计参考站信息
		if (numstas[isys].find(prn) == numstas[isys].end()) {
			EleStation tmp;
			tmp._prn = prn;
			tmp._num = 1;
			tmp._el  = pSat.second._azel[1];
			tmp._rot = pSat.second._rot;
			numstas[isys].emplace(prn, tmp);
		}
		else {
			//TODO: 判定rot是否小于3，才插入
			numstas[isys][prn]._el  += pSat.second._azel[1];
			numstas[isys][prn]._rot += pSat.second._rot;
			numstas[isys][prn]._num += 1;
		}
	}
}

bool LocalAtmoModel::inputAtmoSite(IN Gtime tnow, IN AtmoInfo& stecinf, OUT map<int, int>* track, OUT map<int, EleStation>* numstas)
{
	bool stat = false;
	// 定义新数据
	AtmoEpoch atmoGEC(tnow, 0);
	AtmoEpoch atmoR(tnow, 0);
	// 统计各系统卫星的测站数
	fixDecision(stecinf, track);
	// 遍历测站导入数据
	for (auto& pSta : stecinf._sitesols) {
		int staid = pSta.first;
		if (_allsites.find(staid) == _allsites.end()) {
			continue;
		}

		AtmoSite station(pSta.second);
		int bGEC = 0, bR = 0;
		// GEC
		if (_proOption._usesys[IDX_GPS] && station._staInfo._supSys[IDX_GPS]) {
			inputAtmoSat(IDX_GPS, pSta.second._satsols, station._satInfos, numstas);
			bGEC |= 1;
		}
		if (_proOption._usesys[IDX_GAL] && station._staInfo._supSys[IDX_GAL]) {
			inputAtmoSat(IDX_GAL, pSta.second._satsols, station._satInfos, numstas);
			bGEC |= 1;
		}
		if (_proOption._usesys[IDX_BDS2] && station._staInfo._supSys[IDX_BDS2]) {
			inputAtmoSat(IDX_BDS2, pSta.second._satsols, station._satInfos, numstas);
			bGEC |= 1;
		}
		if (_proOption._usesys[IDX_BDS3] && station._staInfo._supSys[IDX_BDS3]) {
			inputAtmoSat(IDX_BDS3, pSta.second._satsols, station._satInfos, numstas);
			bGEC |= 1;
		}
		if (bGEC) {
			stat |= true;
			atmoGEC._staAtmos.emplace(station._staInfo._name, station);
		}
		// R
		if (_proOption._usesys[IDX_GLO] && station._staInfo._supSys[IDX_GLO]) {
			inputAtmoSat(IDX_GLO, pSta.second._satsols, station._satInfos, numstas);
			bR |= 1;
		}
		if (bR) {
			for (int i = 0; i < NUMSYS; i++) {
				if (i == 1) { continue; }
				station._satInfos[i].clear();
			}
			stat |= true;
			atmoR._staAtmos.emplace(station._staInfo._name, station);
		}
	}
	_stecGroupGEC.emplace(tnow, atmoGEC);
	_stecGroupR.  emplace(tnow, atmoR);

	return stat;
}

bool LocalAtmoModel::inputAtmoEpoch(IN Gtime tnow, IN AtmoInfo& stecinf, IN bool selectRef)
{
	int refsat[NUMSYS] = { 0,0,0,0,0 };
	map<int, int> trackStas[NUMSYS];		// sys->prn->num
	map<int, EleStation> numStas[NUMSYS];	// sys->prn->el
	string tstr = strtime(tnow, 2);
	int ep[6] = { 0 };
	time2epoch(tnow, ep);

	/* 1.预处理 */
	// 前历元数据是否重复
	if (_stecGroupGEC.find(tnow) != _stecGroupGEC.end() ||
		_stecGroupR.  find(tnow) != _stecGroupR.  end()) {
		return false;
	}
	// 剔除最早数据
	while (_stecGroupGEC.size() > MAX_EPOCH_STORE) { _stecGroupGEC.erase(_stecGroupGEC.begin()); }
	while (_stecGroupR.  size() > MAX_EPOCH_STORE) { _stecGroupR.  erase(_stecGroupR.  begin()); }
	// 卫星列表清零
	for (int i = 0; i < NUMSYS; i++) { _stecPro._satList[i].clear(); }

	/* 2.数据输入 */
	if (!inputAtmoSite(tnow, stecinf, trackStas, numStas)) {
		return false;
	}

	/* 3.根据共视测站数剔除卫星(50%) */
	for (int isys = 0; isys < NUMSYS; isys++) {
		for (auto& pSat : trackStas[isys]) {
			int prn = pSat.first;
			if ((isys != 1 && pSat.second < _stanumGEC * THRES_USESTA_PCT) ||
				(isys == 1 && pSat.second < _stanumR   * THRES_USESTA_PCT)) {
				//char SYS = idx2sys(isys);
				//printf("%s erase %c%02d\n", tstr.c_str(), SYS, prn);
				_stecPro._satList[isys].erase(prn);
				numStas[isys].erase(prn);
			}
		}
	}

	/* 4.选共视卫星 */
	if (!selectRef) { return true; }
	if (selectOneRefSat(numStas, refsat)) {
		uniformDatum(tnow, refsat, 0);
		uniformDatum(tnow, refsat, 1);
	}

	return true;
}

bool LocalAtmoModel::doStecModSys(IN int symbol)
{
	AtmoEpoch proAtmo;
	AtmoEpochs* groupAtmo = symbol == 0 ? &_stecGroupGEC : &_stecGroupR;

	/* 设置当前系统 */
	if (!_stecPro.setCurSys(_proOption._usesys, symbol)) {
		//printf("1\n");
		return false;
	}

	/* 数据预处理 */
	if (!_stecPro.preCheckSatModel(_stecPro._tnow, *groupAtmo, proAtmo)) {
		//printf("2\n");
		return false;
	}

	/* 初始化STEC模型&残差 */
	_stecPro.initStecMod(proAtmo);

	/* stec单星建模 */
	if (!_stecPro.satModEst(proAtmo, _gridinfo)) {
		//printf("3  %2d\n", _stecPro._satList[0].size());
		return false;
	}

	if (_stecPro._stecRoti._badroti > 4) {
		_nbadroti++;
	}

	return true;
}

void LocalAtmoModel::resetRefSatBD2(IO StecModEpoch& mod)
{
	map<double, int> QIlist;
	
	if (mod._refsat[IDX_BDS2] > 0) {
		mod._stecmodgnss[IDX_BDS2].erase(mod._refsat[IDX_BDS2]);
		mod._satNum[IDX_BDS2] -= 1;
	}

	for (const auto& pSat : mod._stecmodgnss[IDX_BDS2]) {
		QIlist.emplace(pSat.second._QI[0], pSat.first);
	}

	for (auto& pSat : QIlist) {
		if (pSat.first > 0.0) {
			mod._refsat[IDX_BDS2] = pSat.second;
			break;
		}
	}
}

void LocalAtmoModel::setSatResLevel(IO StecModEpoch& mod)
{
	for (int isys = 0, nres = 0; isys < NUMSYS; isys++) {
		if (_proOption._ressys[isys]) {
			map<double, int> QIlist;
			for (const auto& pSat : mod._stecmodgnss[isys]) {
				QIlist.emplace(pSat.second._QI[1], pSat.first);
			}

			for (auto& pSat : QIlist) {
				if (nres < _proOption._maxsatres) {
					mod._stecmodgnss[isys][pSat.second]._satreslevel = 1;
					nres++;
				}
			}
		}
	}
}

void LocalAtmoModel::copyStecSat(IN int sys, IN int ref, IN map<int, StecModSat>& src, OUT map<int, ProStecModSat>& dst)
{
	for (const auto& pSat : src) {
		int prn = pSat.first;
		const StecModSat& dat = pSat.second;

		ProStecModSat tmp;
		tmp._sys = dat._system;
		tmp._sat = dat._sat;
		tmp._nsta = dat._nsta;
		for (int i = 0; i < 4; i++) {
			tmp._coff[i] = dat._coeff[i];
			tmp._coff_res[i] = dat._coeff_rms[i];
		}
		for (int i = 0; i < 2; i++) {
			tmp._QI[i] = dat._QI[i];
		}
		if (_proOption._algotype) {
			if (prn == ref) {
				if (sys == IDX_BDS2) {
					tmp._QI[0] = tmp._QI[1] = 0.08;
				}
				else {
					tmp._QI[0] = tmp._QI[1] = 0.05;
				}
			}
			else {
				if (sys == IDX_BDS2) {
					tmp._QI[0] = tmp._QI[0] > 0.125 ? tmp._QI[0] : 0.15;
					tmp._QI[1] = tmp._QI[1] > 0.125 ? tmp._QI[1] : 0.15;
				}
				else {
					tmp._QI[0] = tmp._QI[0] > 0.0625 ? tmp._QI[0] : 0.1;
					tmp._QI[1] = tmp._QI[1] > 0.0625 ? tmp._QI[1] : 0.1;
				}
			}
		}
		else {
			if (prn == ref) {
				if (sys == IDX_BDS2) {
					tmp._QI[0] = tmp._QI[1] = 0.35;
				}
				else {
					tmp._QI[0] = tmp._QI[1] = 0.35;
				}
			}
			else {
				if (sys == IDX_BDS2) {
					tmp._QI[0] = tmp._QI[0] > 0.5 ? tmp._QI[0] : 0.55;
					tmp._QI[1] = tmp._QI[1] > 0.5 ? tmp._QI[1] : 0.55;
				}
				else {
					tmp._QI[0] = tmp._QI[0] > 0.3125 ? tmp._QI[0] : 0.4;
					tmp._QI[1] = tmp._QI[1] > 0.3125 ? tmp._QI[1] : 0.4;
				}
			}
		}
		tmp._satreslevel = dat._satreslevel;
		tmp._gridNum = dat._gridNum;
		tmp._ele = dat._ele * R2D;
		for (int i = 0; i < tmp._gridNum; i++) {
			tmp._stecpergrid[i]._gridid = dat._stecpergrid.at(i + 1)._gridid;
			tmp._stecpergrid[i]._lat    = dat._stecpergrid.at(i + 1)._lat;
			tmp._stecpergrid[i]._nsta   = dat._stecpergrid.at(i + 1)._nsta;
			tmp._stecpergrid[i]._lon    = dat._stecpergrid.at(i + 1)._lon;
			tmp._stecpergrid[i]._stec   = dat._stecpergrid.at(i + 1)._stec;
			tmp._stecpergrid[i]._gridreslevel = dat._satreslevel > 0 ? 1 : 0;
		}

		dst.emplace(prn, tmp);
	}
}

void LocalAtmoModel::copyStecMod(IN StecModEpoch& src, OUT ProStecMod& dst)
{
	dst._time = src._time;
	dst._modetype = src._modetype;
	for (int isys = 0; isys < NUMSYS; isys++) {
		dst._satNum[isys] = src._satNum[isys];
		dst._satRef[isys] = dst._satNum[isys] > 0 ? src._refsat[isys] : 0;

		copyStecSat(isys, dst._satRef[isys], src._stecmodgnss[isys], dst._stecmod[isys]);
	}
}

bool LocalAtmoModel::doStecMod(IN Gtime tnow, IN AtmoInfo& stecinf, OUT ProStecMod& stecmod)
{
	bool stat = false;
	string tstr = strtime(tnow, 2);
	int ep[6] = { 0 };
	time2epoch(tnow, ep);
	_stecPro.settime(tnow);

	/* 1.输入当前历元的STEC数据 */
	if (!inputAtmoEpoch(tnow, stecinf, true)) {
		return false;
	}

	/* 2.STEC建模 */
	if (_stanumGEC > 0) { stat |= doStecModSys(0); }
	if (_stanumR   > 0) { stat |= doStecModSys(1); }
	if (!stat) {
		//printf("%s do stecmod fail\n", tstr.c_str());
		return false;
	}

	/* 3.参考星平滑 */
	if (_stecPro._refsatsmooth) {
		_stecPro.refSatSmoothing(tnow);
	}

	/* 4.更新历史建模数据 */
	_stecPro.addStecMod(tnow);

	/* 5.若采用混合解建模，BD2重新选择QI最小的卫星作为参考星 */
	StecModEpoch stecModNow = _stecPro._stecModCur;
	if (_proOption._algotype == 0) {
		resetRefSatBD2(stecModNow);
	}

	/* 6.遍历系统/卫星，设置残差等级 */
	this->setSatResLevel(stecModNow);

	/* 7.保存当前建模信息 */
	this->copyStecMod(stecModNow, stecmod);

	return true;
}

int LocalAtmoModel::inputZtdEpoch(IN Gtime tnow, IN AtmoInfo& ztdinf)
{
	string tstr = strtime(tnow, 2);
	double maxlat = -999, maxlon = -999, minlat = 999, minlon = 999;
	/* 1.导入数据 */
	for (auto& iSta : ztdinf._sitesols) {
		ZtdInfo inf;
		inf._time    = iSta.second._time;
		inf._name    = iSta.second._name;
		inf._ID      = iSta.second._ID;
		inf._zhd     = iSta.second._zhd;
		inf._zwd     = iSta.second._zwd;
		inf._std_zwd = iSta.second._std_zwd;
		for (int i = 0; i < 3; i++) {
			inf._xyz[i] = iSta.second._xyz[i];
			inf._blh[i] = iSta.second._blh[i];
		}

		if (inf._blh[0] > _gridinfo._latcell[0] * D2R &&
			inf._blh[0] < _gridinfo._latcell[1] * D2R &&
			inf._blh[1] > _gridinfo._loncell[0] * D2R &&
			inf._blh[1] < _gridinfo._loncell[1] * D2R)
		{
			_ztds.emplace(inf._name, inf);
		}
	}
	/* 2.站点分布稀疏化 */
	if (_proOption._bsparse) { sparseSites(); }

	/* 3.统计站点分布 */
	//for (auto& iSta : _ztds) {
	//	if (iSta.second._blh[0] < minlat) { minlat = iSta.second._blh[0]; }
	//	if (iSta.second._blh[0] > maxlat) { maxlat = iSta.second._blh[0]; }
	//	if (iSta.second._blh[1] < minlon) { minlon = iSta.second._blh[1]; }
	//	if (iSta.second._blh[1] > maxlon) { maxlon = iSta.second._blh[1]; }
	//}
	//printf("%s %7.2f %7.2f %7.2f %7.2f\n", tstr.c_str(), minlat * R2D, maxlat * R2D, minlon * R2D, maxlon * R2D);
	
	return (int)_ztds.size();
}

int LocalAtmoModel::sparseSites()
{
	int k = 0;
	map<double, string> dists;
	ZtdInfos ztd_sp;

	for (auto& iCel : _gridinfo._cells) {
		k = 0;
		dists.clear();
		//printf("%2d: ", iCel.first);
		
		for (auto& iSta : _ztds) {
			if (iSta.second._blh[0] > iCel.second[0] && iSta.second._blh[0] < iCel.second[1] &&
				iSta.second._blh[1] > iCel.second[3] && iSta.second._blh[1] < iCel.second[4]) 
			{
				double dist = sphereDist(iSta.second._blh[0], iSta.second._blh[1], iCel.second[2], iCel.second[5]);
				dists.emplace(dist, iSta.first);
			}
		}

		for (auto iDist : dists) {
			if (k++ < 2) {
				ztd_sp.emplace(iDist.second, _ztds[iDist.second]);
				//printf("%4s ", iDist.second.c_str());
			}
			_ztds.erase(iDist.second);
		}
		//printf("\n");
	}
	_ztds = ztd_sp;

	return (int)_ztds.size();
}

void LocalAtmoModel::copyZtdMod(IN int npara, IN ZtdModEpoch& src, OUT ProZtdMod& dst)
{
	dst._time = src._time;
	dst._zhd  = src._zhd;
	dst._qi   = src._qi;
	dst._ncoeff = npara;
	dst._nsta = src._nsta;
	for (int i = 0; i < npara; i++) {
		dst._coeff[i] = src._coeff[i];
		dst._coeff_rms[i] = src._coeff_rms[i];
	}
}

bool LocalAtmoModel::doZtdMod(IN Gtime tnow, IN AtmoInfo& ztdinf, OUT ProZtdMod& ztdmod)
{
	int nsit = 0, bmean = _proOption._meanzhd;
	bool stat = false;
	string tstr = strtime(tnow, 2);

	_ztdPro.settime(tnow);
	_ztds.clear();
	
	/* 1.输入当前历元的ZTD数据 */
	if (!(nsit = inputZtdEpoch(tnow, ztdinf))) {
		return false;
	}

	/* 2.ztd对齐至海平面msl */
	_ztdPro.ztd2msl(bmean, _ztds);

	/* 3.zwd建模 */
	stat = _ztdPro.ofcModel(tnow, _ztds, _gridinfo, bmean);
	//if (!stat) { cout << tstr << " ofc fail " << _ztds.size() << endl; }

	/* 4.保存当前建模信息 */
	ZtdModEpoch ztdModNow = _ztdPro._ztdModCur;
	this->copyZtdMod(_ztdPro._ncoeff, ztdModNow, ztdmod);

	return stat;
}