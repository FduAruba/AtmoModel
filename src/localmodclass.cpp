#include "localmodclass.h"

void AtmoInfo::reset()
{
	_time = { 0 };
	_sitenum = 0;
	_sitesols.clear();
}

void ProStecGrid::deepcopy(IN ProStecGrid& src)
{
	this->_gridid = src._gridid;
	this->_nsta = src._nsta;
	this->_lat = src._lat;
	this->_lon = src._lon;
	this->_stec = src._stec;
	this->_gridreslevel = src._gridreslevel;
}

void ProStecGrid::reset()
{
	this->_gridid = 0;
	this->_nsta = 0;
	this->_lat = 0.0;
	this->_lon = 0.0;
	this->_stec = 0.0;
	this->_gridreslevel = 0;
}

void ProStecModSat::deepcopy(IN ProStecModSat& src)
{
	this->_sys = src._sys;
	this->_sat = src._sat;
	this->_nsta = src._nsta;
	this->_ele = src._ele;
	for (int i = 0; i < 4; i++) {
		this->_coff[i] = src._coff[i];
		this->_coff_res[i] = src._coff_res[i];
	}
	for (int i = 0; i < 2; i++) {
		this->_QI[i] = src._QI[i];
	}
	this->_gridNum = src._gridNum;
	for (int i = 0; i < this->_gridNum; i++) {
		this->_stecpergrid[i].deepcopy(src._stecpergrid[i]);
	}
	this->_satreslevel = src._satreslevel;
}

void ProStecModSat::reset()
{
	this->_sys = '\0';
	this->_sat = 0;
	this->_nsta = 0;
	this->_ele = 0.0;
	for (int i = 0; i < 4; i++) {
		this->_coff[i] = 0.0;
		this->_coff_res[i] = 0.0;
	}
	for (int i = 0; i < 2; i++) {
		this->_QI[i] = 0;
	}
	for (int i = 0; i < this->_gridNum; i++) {
		this->_stecpergrid[i].reset();
	}
	this->_gridNum = 0;
	this->_satreslevel = 0;
}

void ProStecMod::reset()
{
	_time = { 0 };
	_modetype = 1;
	for (int i = 0; i < NUMSYS; i++) {
		_satNum[i] = _satRef[i] = 0;
		_stecmod[i].clear();
	}
}

void ProZtdMod::reset()
{
	_time = { 0 };
	_zhd = 0;
	_nsta = _ncoeff = _qi = 0;
	for (int i = 0; i < 10; i++) {
		_coeff[i] = _coeff_rms[i] = 0;
	}
}

SatInfo::SatInfo(IN SatSol src)
{
	_sys     = src._sys;
	_prn     = src._prn;
	_fixflag = src._fixflag;
	_iono    = src._iono;
	_covIono = src._covIono;
	_trop    = src._trop;
	_covTrop = src._covTrop;
	_quality = src._quality;
	_rot     = src._rot;
	_wgt     = 0.0;
	for (int i = 0; i < 2; i++) {
		_azel[i] = src._azel[i];
	}
	for (int i = 0; i < 3; i++) {
		_xyz[i] = src._xyz[i];
	}
}

void SiteInfo::coypSiteSol(IN SiteSol src)
{
	_name = src._name;
	_ID = src._ID;
	for (int i = 0; i < 3; i++) {
		_xyz[i] = src._xyz[i];
		_blh[i] = src._blh[i];
	}
	for (int i = 0; i < NUMSYS; i++) {
		_satNum[i] = src._satNum[i];
		_supSys[i] = src._satNum[i] > 0 ? 1 : 0;
	}
}

AtmoSite::AtmoSite(IN SiteSol sitesol)
{
	_pppflag = 6;
	for (int i = 0; i < NUMSYS; i++) {
		_refSat[i] = 0;
	}
	_staInfo.coypSiteSol(sitesol);
}

AtmoEpoch::AtmoEpoch(IN Gtime t, IN int nsta)
{
	_time = t;
	time2epoch(t, _ep);
	_stanum = nsta;
	memset(_refSat, 0, sizeof(_refSat));
}

void GridInfo::deepcopy(const GridInfo& src)
{
	for (int i = 0; i < 2; i++) {
		_latgrid[i] = src._latgrid[i];
		_longrid[i] = src._longrid[i];
		_latcell[i] = src._latcell[i];
		_loncell[i] = src._loncell[i];
		_step[i]    = src._step[i];
		_center[i]  = src._center[i];
		_length[i]  = src._length[i];
	}

	_gridNum = src._gridNum;
	for (int i = 0; i < _gridNum; i++) {
		_grids[i] = src._grids[i];
	}
}

void ProOption::deepcopy(IN ProOption& src)
{
	_ti = src._ti;
	for (int i = 0; i < NUMSYS; i++) {
		_usesys[i] = src._usesys[i];
		_fixsys[i] = src._fixsys[i];
		_ressys[i] = src._ressys[i];
	}
	_maxsatres    = src._maxsatres;
	_algotype     = src._algotype;
	_ionotype     = src._ionotype;
	_diffmode     = src._diffmode;
	_refsatsmooth = src._refsatsmooth;
	_minel        = src._minel;
	for (int i = 0; i < 5; i++) {
		_maxroti[i] = src._maxroti[i];
	}
	_qimulti = src._qimulti;
	_qibase  = src._qibase;
	_qicoeff = src._qicoeff;
	_troptype = src._troptype;
	_bsparse = src._bsparse;
	_meanzhd = src._meanzhd;
}

void SiteGridDist::calcAllDist(IN AtmoEpoch& atmo, IN GridInfo& grid)
{
	for (const auto& pSta : atmo._staAtmos) {
		string site = pSta.first;
		double latB = pSta.second._staInfo._blh[0];
		double lonB = pSta.second._staInfo._blh[1];

		if (_dist.find(1) != _dist.end() && 
			_dist[1].find(site) != _dist[1].end()) {
			continue;
		}

		for (int i = 0; i < grid._gridNum; i++) {
			int id      = grid._grids[i]._id;
			double latG = grid._grids[i]._lat;
			double lonG = grid._grids[i]._lon;
			_dist[id][site] = sphereDist(latG, lonG, latB, lonB);
		}
	}
}

double SiteGridDist::getDist(IN int id, IN string site)
{
	auto it = _dist.find(id);
	if (it != _dist.end()) {
		auto it_site = it->second.find(site);
		if (it_site != it->second.end()) {
			return it_site->second;
		}
	}
	return 0.0;
}

void SatRoti::reinit(IN Gtime tnow, IN double stec)
{
	_stec.clear();
	_rot.clear();
	_roti = 0.0;
	_tlast = tnow;
	_stec.emplace(tnow, stec);
}

double SatRoti::proRoti(IN Gtime tnow, IN SatAtmo& satinf)
{
	_rot.clear();
	_tlast = tnow;
	
	if (_stec.size() < 1) {
		_stec.emplace(tnow, satinf._iono);
		return 0.0;
	}
	
	while (_stec.size() > 10) {
		_stec.erase(_stec.begin());
	}
	
	Gtime tlast = _stec.rbegin()->first;
	string t0 = strtime(tlast, 2);
	string t1 = strtime(tnow, 2);
	if (tnow <= tlast || fabs(tnow - tlast) > 120.0) {
		this->reinit(tnow, satinf._iono);
		//printf("tnow=%s tlast=%s \n", t1.c_str(), t0.c_str());
		return 0.0;
	}

	for (auto it = _stec.begin(); it != _stec.end();) {
		if (fabs(tnow - it->first) >= 240.0) {
			it = _stec.erase(it);
		}
		break;
	}

	_stec.emplace(tnow, satinf._iono);

	for (auto it = _stec.begin(); it != _stec.end(); ++it) {
		auto in = it; ++in;
		if (in != _stec.end()) {
			double rot = (in->second - it->second) * 60.0 / (in->first - it->first);
			_rot.emplace(in->first, rot);
		}
	}

	int sz = _rot.size();
	double avgsum = 0.0, argsum = 0.0, roti = 0.0;
	for (auto it = _rot.begin(); it != _rot.end(); ++it) {
		avgsum += it->second;
	}
	avgsum /= sz;
	for (auto it = _rot.begin(); it != _rot.end(); ++it) {
		argsum += pow(it->second - avgsum, 2);
	}
	roti = sqrt(argsum / sz);

	_roti = roti;

	return roti;
}

void StaRoti::proRoti(IN Gtime tnow, IN SiteAtmo& siteinf)
{
	string tstr = strtime(tnow, 2);
	
	for (int isys = 0; isys < NUMSYS; isys++) {
		char SYS = idx2sys(isys);
	
		if (siteinf._satIon[isys].size() <= 0) {
			continue;
		}

		for (auto& pSat : siteinf._satIon[isys]) {
			int prn = pSat.first;
			
			if (_rotis[isys].find(prn) == _rotis[isys].end()) {
				SatRoti sat(tnow, prn);
				double roti = sat.proRoti(tnow, pSat.second);
				_rotis[isys].emplace(prn, sat);
				//printf("%s %c%02d %5.2f\n", tstr.c_str(), SYS, prn, roti);
			}
			else {
				SatRoti& sat = _rotis[isys].find(prn)->second;
				double roti = sat.proRoti(tnow, pSat.second);
				//printf("%s %c%02d %5.2f\n", tstr.c_str(), SYS, prn, roti);
			}
		}
	}
}