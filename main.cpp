#include <iostream>
#include <bitset>
#include "outfiles.h"

using namespace std;

void rovRoti(IN Gtime tnow, IN SiteAtmos& rovinf, OUT StaRotiMap& rovs)
{
	for (auto& pSta : rovinf) {
		string site = pSta.second._name;

		if (rovs.find(site) == rovs.end()) {
			StaRoti sta(site);
			sta.proRoti(tnow, pSta.second);
			rovs.emplace(site, sta);
		}
		else {
			StaRoti& sta = rovs[site];
			sta.proRoti(tnow, pSta.second);
		}
	}
}

void outRovRoti(IN Gtime tnow, IN StaRotiMap& rovs, OUT FileFps& fps)
{
	char buff0[MAXCHARS] = { '\0' }, * p = buff0;
	char buff1[MAXOUTCHARS] = { '\0' }, * q = buff1;
	string tstr = strtime(tnow, 1);

	for (auto& pSta : rovs) {
		string site = pSta.first;

		int n = 0;
		for (int isys = 0; isys < NUMSYS; isys++) {
			char SYS = idx2sys(isys);

			for (auto& pSat : pSta.second._rotis[isys]) {
				int prn = pSat.first;
				Gtime t = pSat.second._tlast;
				double rot = pSat.second._rot.size() > 0 ? pSat.second._rot.rbegin()->second : 0.0;
				double roti = pSat.second._roti;

				if (t == tnow && roti != 0.0) {
					q += sprintf(q, "%c%02d %8.2f %8.2f\n", SYS, prn, rot, roti);
					n++;
				}
			}
		}

		if (n > 0) {
			p += sprintf(p, ">%s %02d\n", tstr.c_str(), n);
			fwrite(buff0, (int)(p - buff0), sizeof(char), fps[site][5]);
			fwrite(buff1, (int)(q - buff1), sizeof(char), fps[site][5]);
		}
		p = buff0; q = buff1;
	}
}

void outdebug(IN Coption& config, IN LocalAtmoModel* localMod)
{
	for (auto pSta : config._sta) {
		double r1 = pSta._nepo > 0 ? (double)pSta._nbad / (double)pSta._nepo * 100.0 : 0.0;
		double r2 = pSta._nepo > 0 ? (double)pSta._nepo / 17280.0 * 100.0 : 0.0;
		printf("%s #%2d nep=%6d nbad=%6d r1=%6.2f%% r2=%6.2f%%\n", pSta._name.c_str(),
			pSta._ID, pSta._nepo, pSta._nbad, r1, r2);
	}
	printf("nlack=%5d\n", config._nlack);
	printf("nbadroti=%5d\n", localMod->_nbadroti);

	for (auto it = config._rovstatic.begin(); it != config._rovstatic.end(); ++it) {
		string rov_t = it->first;
		double c1 = it->second[0] / it->second[3] * 100.0;
		double c2 = it->second[1] / it->second[3] * 100.0;
		double c3 = it->second[2] / it->second[3] * 100.0;
		printf("%s nall=%7d ngood=%7d nbad=%7d nout=%7d %.2f%% %5.2f%% %5.2f%%\n",
			rov_t.c_str(), (int)it->second[3], (int)it->second[0], 
			(int)it->second[1], (int)it->second[2], c1, c2, c3);
	}

	double r1 = (1.0 * config._ngoodres) / (1.0 * config._nvali) * 100.0;
	double r2 = (1.0 * config._nbadres)  / (1.0 * config._nvali) * 100.0;
	double r3 = (1.0 * config._noutres)  / (1.0 * config._nvali) * 100.0;
	printf("---- nall=%7d ngood=%7d nbad=%7d nout=%7d %.2f%% %5.2f%% %5.2f%%\n",
		config._nvali, config._ngoodres, config._nbadres, config._noutres, r1, r2, r3);
}

int main()
{
	int dbg = 1;
	double t1, t2;
	char fname[MAXCHARS] = "C:\\Users\\shuqh\\Desktop\\rsim.cfg";
	Coption config;									// 配置文件
	ProOption* popt          = new ProOption;		// 处理设置
	GridInfo* grid           = new GridInfo;		// 格网点信息
	AtmoInfo* stecinf        = new AtmoInfo;		// stec数据
	ProStecMod* stecmod      = new ProStecMod;		// stec建模
	LocalAtmoModel* localMod = new LocalAtmoModel;	// 大气建模类
	FileFps outfps;
	StaRotiMap stamap;
	
	t1 = clock();
	/* 读取配置文件 */
	if (readConfigFile(fname, config)) {
		movOption(config, *popt);
		movGrids(config, *grid);
	}
	else {
		return -1;
	}
	
	Gtime ts = epoch2time(config._ts);
	Gtime te = epoch2time(config._te);
	for (Gtime t = ts; t <= te; t.time += (time_t)config._ti) {
		string tstr = strtime(t, 2);
		printf("\r%s: processing...%c", tstr.c_str(), t == te ? '\n' : '\0');

		/* 读取单站数据 */
		SiteAtmos stas, rovs;
		if (!readAugmentData(t, config, stas, rovs)) {
			continue;
		}

		/* 分配历元数据 */
		stecinf->reset(); stecmod->reset();
		if (!movAtmos(t, stas, *stecinf)) {
			continue;
		}

		/* 历元基础设置 */
		if (t == ts) {
			localMod->setUseres(config._useres);
			localMod->setReigon(*grid);
			localMod->setOption(*popt);
			localMod->_stecPro.setBasicOption(*popt, config._useres);
			createRovFile(config, outfps);
			writeModelHead(config, outfps["MODEL"][1]);
		}
		localMod->setRefSites(stas);

		/* 验证站ROTI计算输出 */
		rovRoti(t, rovs, stamap);
		outRovRoti(t, stamap, outfps);
		
		/* 大气建模 */
		if (config._modeltype & 1) {
			if (localMod->doStecMod(t, *stecinf, *stecmod)) {
				outRovStec(config, *grid, rovs, *stecmod, outfps, 1);
				outStecModel(t, config, *grid, *stecmod, outfps["MODEL"][1]);
			}
		}
		if (config._modeltype & 2) { //TODO: ZTD建模
		}
		if (config._modeltype & 4) { //TODO: STD建模
		}
		if (config._modeltype & 8) { //TODO: VTEC建模
		}
	}

	if (dbg) { outdebug(config, localMod); }

	t2 = clock();
	printf("\n* Total program time: %4.1f minute\n", (t2 - t1) / CLOCKS_PER_SEC / 60.0);

	/* 释放内存空间 */
	delete popt; delete grid; delete stecinf; delete stecmod; delete localMod;

	return 0;
}