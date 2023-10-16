#include <iostream>
#include <bitset>
#include "outfiles.h"

using namespace std;

int main()
{
	double t1, t2;
	char fname[MAXCHARS] = "C:\\Users\\shuqh\\Desktop\\rsim.cfg";
	Coption config;									// �����ļ�
	ProOption* popt          = new ProOption;		// ��������
	GridInfo* grid           = new GridInfo;		// ��������Ϣ
	AtmoInfo* stecinf        = new AtmoInfo;		// stec����
	ProStecMod* stecmod      = new ProStecMod;		// stec��ģ
	LocalAtmoModel* localMod = new LocalAtmoModel;	// ������ģ��
	FileFps outfps;

	t1 = clock();
	/* ��ȡ�����ļ� */
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
		/* ��ȡ��վ���� */
		SiteAtmos stas, rovs;
		if (!readAugmentData(t, config, stas, rovs)) {
			continue;
		}

		/* ������Ԫ���� */
		stecinf->reset(); stecmod->reset();
		if (!movAtmos(t, stas, *stecinf)) {
			continue;
		}

		/* ��Ԫ�������� */
		if (t == ts) {
			localMod->setUseres(config._useres);
			localMod->setReigon(*grid);
			localMod->setOption(*popt);
			localMod->_stecPro.setBasicOption(*popt, config._useres);
			createRovFile(config, outfps);
		}
		localMod->setRefSites(stas);
		
		/* ������ģ */
		if (config._modeltype & 1) {
			if (localMod->doStecMod(t, *stecinf, *stecmod)) {
				outRovStec(config, *grid, rovs, *stecmod, outfps, 1);
			}
		}
		if (config._modeltype & 2) { //TODO: ZTD��ģ
		}
		if (config._modeltype & 4) { //TODO: STD��ģ
		}
		if (config._modeltype & 8) { //TODO: VTEC��ģ
		}

		/* �����ģ��� */

	}

#if 1
	for (auto pSta : config._sta) {
		double r1 = pSta._nepo > 0 ? (double)pSta._nbad / (double)pSta._nepo * 100.0 : 0.0;
		double r2 = pSta._nepo > 0 ? (double)pSta._nepo / 17280.0 * 100.0 : 0.0;
		printf("%s #%2d nep=%6d nbad=%6d r1=%6.2f%% r2=%6.2f%%\n", pSta._name.c_str(), pSta._ID, pSta._nepo, pSta._nbad, r1, r2);
	}
	printf("nlack=%5d\n", config._nlack);
	printf("nbadroti=%5d\n", localMod->_nbadroti);

	for (auto it = config._rovstatic.begin(); it != config._rovstatic.end(); ++it) {
		string rov_t = it->first;
		double c1 = it->second[0] / it->second[3]*100.0;
		double c2 = it->second[1] / it->second[3]*100.0;
		double c3 = it->second[2] / it->second[3]*100.0;
		printf("%s nall=%7d ngood_res=%7d nbad_res=%7d nout=%7d %.2f%% %5.2f%% %5.2f%%\n",
			rov_t.c_str(),(int)it->second[3], (int)it->second[0], (int)it->second[1], (int)it->second[2],
			c1, c2, c3);
	}

	double r1 = (1.0 * config._ngoodres) / (1.0 * config._nvali) * 100.0;
	double r2 = (1.0 * config._nbadres)  / (1.0 * config._nvali) * 100.0;
	double r3 = (1.0 * config._noutl)    / (1.0 * config._nvali) * 100.0;
	printf("---- nall=%7d ngood_res=%7d nbad_res=%7d nout=%7d %.2f%% %5.2f%% %5.2f%%\n", 
		config._nvali, config._ngoodres, config._nbadres, config._noutl, r1, r2, r3);
#endif

	t2 = clock();
	printf("\n* Total program time: %4.1f minute\n", (t2 - t1) / CLOCKS_PER_SEC / 60.0);
	/* �ͷ��ڴ�ռ� */
	delete popt; delete grid; delete stecinf; delete stecmod; delete localMod;

	return 0;
}