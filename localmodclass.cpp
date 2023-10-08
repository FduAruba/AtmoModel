#include "localmodclass.h"

void AtmoInfo::reset()
{
	_time = { 0 };
	_sitenum = 0;
	_sitesols.clear();
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
	_diffmode     = src._diffmode;
	_refsatsmooth = src._refsatsmooth;
	_minel        = src._minel;
	for (int i = 0; i < 5; i++) {
		_maxroti[i] = src._maxroti[i];
	}
	_qimulti = src._qimulti;
	_qibase  = src._qibase;
	_qicoeff = src._qicoeff;
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