#include "localmodel.h"

void LocalAtmoModel::setUseres(IN int res)
{
	_useres = res;
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
	for (auto pSta : stas) {
		int id = pSta.first;

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
		// 1>�ɲο��ǵĹ��Ӳ�վ��>��ǰ�����Ҹ߶Ƚ�>40
		if (pold != numstas[isys].end()) {
			pctold = (1.0 * pold->second._num) / stanum;
			pctnew = (1.0 * pNum->first) / stanum;
			if (pctold - pctnew > 0.1 && pold->second._el > 40.0 * D2R) {
				break;
			}
		}
		// 2>�ҳ��߶Ƚ�>40���ҵ�ǰ����վ��������
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
					//TODO: rot�ж�
					if (pctnew - pctold > 0.1 || pold->second._el < 40.0 * D2R) {
						if (modLastEpoch != NULL && modLastEpoch->_stecmodgnss[isys].find(prn) != modLastEpoch->_stecmodgnss[isys].end()) {
							refsat[isys] = prn;
							isfind = true;
						}
					}
				}
			}
		}
		// 3>���������Ǿ��߶Ƚ�<40��ֱ��ȡ��๲��վ��������Ϊ�ο���
		if (refsat[isys] == 0 || pctnew - pctold > 0.1) {
			refsat[isys] = vprn.rbegin()->second.rbegin()->_prn;
			pold = numstas[isys].find(refsat[isys]);
		}
	}
	// 4>����վ�Ȳ��û��>10%��ֱ��ȡ��๲��վ��������Ϊ�ο���
	if (!isfind && refsat[isys] == 0) {
		refsat[isys] = vprn.rbegin()->second.rbegin()->_prn;
	}
}

void  LocalAtmoModel::uniformDatum(IN Gtime tnow, IN int* refsat, IN int symbol)
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
		// �������Ƕ�Ӧ��վ����
		vPrn.clear();
		sortSiteNumber(isys, numstas, vPrn);

		if (vPrn.size() <= 0) {
			printf("no valid sat %d %d\n", isys, (int)numstas[isys].size());
			continue;
		}

		// ɸѡ�ο���
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
		//TODO: �ж�roti�Ƿ������ֵ0.5���к�����ʱ�޷����

		// ����۲�����
		dst[isys].emplace(prn, info);
		// ���������б�
		_stecPro._satList[isys].emplace(prn);
		// ͳ�Ʋο�վ��Ϣ
		if (numstas[isys].find(prn) == numstas[isys].end()) {
			EleStation tmp;
			tmp._prn = prn;
			tmp._num = 1;
			tmp._el  = pSat.second._azel[1];
			tmp._rot = pSat.second._rot;
			numstas[isys].emplace(prn, tmp);
		}
		else {
			//TODO: �ж�rot�Ƿ�С��3���Ų���
			numstas[isys][prn]._el  += pSat.second._azel[1];
			numstas[isys][prn]._rot += pSat.second._rot;
			numstas[isys][prn]._num += 1;
		}
	}
}

bool LocalAtmoModel::inputAtmoSite(IN Gtime tnow, IN AtmoInfo& stecinf, OUT map<int, int>* track, OUT map<int, EleStation>* numstas)
{
	bool stat = false;
	// ����������
	AtmoEpoch atmoGEC(tnow, 0);
	AtmoEpoch atmoR(tnow, 0);
	// ͳ�Ƹ�ϵͳ���ǵĲ�վ��
	fixDecision(stecinf, track);
	// ������վ��������
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
	
	/* 1.Ԥ���� */
	// ǰ��Ԫ�����Ƿ��ظ�
	if (_stecGroupGEC.find(tnow) != _stecGroupGEC.end() ||
		_stecGroupR.  find(tnow) != _stecGroupR.  end()) {
		return false;
	}
	// �޳���������
	while (_stecGroupGEC.size() > MAX_EPOCH_STORE) { _stecGroupGEC.erase(_stecGroupGEC.begin()); }
	while (_stecGroupR.  size() > MAX_EPOCH_STORE) { _stecGroupR.  erase(_stecGroupR.  begin()); }
	// �����б�����
	for (int i = 0; i < NUMSYS; i++) { _stecPro._satList[i].clear(); }

	/* 2.�������� */
	if (!inputAtmoSite(tnow, stecinf, trackStas, numStas)) {
		return false;
	}

	/* 3.���ݹ��Ӳ�վ���޳�����(50%) */
	for (int isys = 0; isys < NUMSYS; isys++) {
		for (auto& pSat : trackStas[isys]) {
			int prn = pSat.first;
			if ((isys != 1 && pSat.second < _stanumGEC * THRES_USESTA_PCT) ||
				(isys == 1 && pSat.second < _stanumR   * THRES_USESTA_PCT)) {
				//char SYS = idx2sys(isys);
				//printf("%s erase %c%02d\n", tstr.c_str(), SYS, prn);
				_stecPro._satList[isys].erase(prn);
			}
		}
	}

	/* 4.ѡ�������� */
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

	/* ���õ�ǰϵͳ */
	if (!_stecPro.setCurSys(_proOption._usesys, symbol)) {
		return false;
	}
	/* ����Ԥ���� */
	if (!_stecPro.preCheckSatModel(_stecPro._tnow, *groupAtmo, proAtmo)) {
		return false;
	}
	/* ��ʼ��STECģ��&�в� */
	_stecPro.initStecMod(proAtmo);
	/* STEC���ǽ�ģ */
	_stecPro.satModEst(proAtmo, _gridinfo);



	if (_stecPro._stecRoti._badroti > 4) {
		_nbadroti++;
	}

	return true;
}

bool LocalAtmoModel::doStecMod(IN Gtime tnow, IN AtmoInfo& stecinf, OUT ProStecMod& stecmod)
{
	_stecPro._tnow = tnow;

	/* 1.���뵱ǰ��Ԫ��STEC����  */
	if (!inputAtmoEpoch(tnow, stecinf, true)) {
		return false;
	}

	/* 2.STEC��ģ */
	if (_stanumGEC > 0) { doStecModSys(0); }
	if (_stanumR   > 0) { doStecModSys(1); }

	/* 3.�����ø���/�̶����Ͻ�ģ��BD2����ѡ��QI��С��������Ϊ�ο��� */

	/* 4.����ϵͳ/���ǣ��������ǲв���Ϣ */

	return true;
}
