#include <iostream>
#include <bitset>
#include "readfiles.h"

using namespace std;

int main()
{
	char fname[MAXCHARS] = "C:\\Users\\QHY\\Desktop\\rsim.cfg";
	Coption config;									// 配置文件
	ProOption* popt          = new ProOption;		// 处理设置
	GridInfo* grid           = new GridInfo;		// 格网点信息
	AtmoInfo* stecinf        = new AtmoInfo;		// stec数据
	ProStecMod* stecmod      = new ProStecMod;		// stec建模
	LocalAtmoModel* localMod = new LocalAtmoModel;	// 大气建模类

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
		/* 读取单站数据 */
		SiteAtmos stas, rovs;
		if (!readAugmentData(t, config, stas, rovs)) {
			continue;
		}

		/* 分配历元数据 */
		stecinf->reset();
		stecmod->reset();
		if (!movAtmos(t, stas, *stecinf)) {
			continue;
		}

		/* 历元基础设置 */
		if (t == ts) {
			localMod->setUseres(config._useres);
			localMod->setReigon(*grid);
			localMod->setOption(*popt);
			localMod->_stecPro.setBasicOption(*popt, config._useres);
		}
		localMod->setRefSites(stas);
		
		/* 大气建模 */
		if (config._modeltype & 1) {
			if (!localMod->doStecMod(t, *stecinf, *stecmod)) {
				continue;
			}
		}


		/* 输出建模结果 */

	}

#if 1
	for (auto pSta : config._sta) {
		double r1 = pSta._nepo > 0 ? (double)pSta._nbad / (double)pSta._nepo * 100.0 : 0.0;
		double r2 = pSta._nepo > 0 ? (double)pSta._nepo / 17280.0 * 100.0 : 0.0;
		printf("%s #%2d nep=%6d nbad=%6d r1=%6.2f%% r2=%6.2f%%\n", pSta._name.c_str(), pSta._ID, pSta._nepo, pSta._nbad, r1, r2);
	}
	printf("nlack=%d\n", config._nlack);
	printf("nbadroti=%d\n", localMod->_nbadroti);
#endif

	/* 释放内存空间 */

	return 0;
}