#include "outfiles.h"

void outRovStec(IN Gtime tnow, IN Coption& cfg, IN GridInfo& grid, IN SiteAtmos& rovaug, IN ProStecMod& stecmod, IN FileFps& rovfps)
{
	char buff[MAXOUTCHARS] = { '/0' }, * p = buff;

	for (const auto& pRov : rovaug) {
		string rov = pRov.second._name;
		if (rovfps.find(rov) == rovfps.end()) {
			//
		}
		
	}
}

void createRovFile(IN Coption& cfg, OUT FileFps& fps)
{
	if (0 != _access(cfg._pathou.c_str(), 0)) {
		_mkdir(cfg._pathou.c_str());
	}
	
	for (const auto& pSta : cfg._rov) {
		string rov = pSta._name;
		string rovpath = cfg._pathou + "\\" + rov;

		if (cfg._modeltype & 1) {
			string stecpath = rovpath + "-stec.txt";
			FILE* fp = fopen(stecpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(1, fp);
			}
		}
		if (cfg._modeltype & 2) {
			string ztdpath = rovpath + "-ztd.txt";
			FILE* fp = fopen(ztdpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(2, fp);
			}
		}
		if (cfg._modeltype & 4) {
			string stdpath = rovpath + "-std.txt";
			FILE* fp = fopen(stdpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(3, fp);
			}
		}
		if (cfg._modeltype & 8) {
			string vtecpath = rovpath + "-vtec.txt";
			FILE* fp = fopen(vtecpath.c_str(), "w");
			if (fp != NULL) {
				fps[rov].emplace(4, fp);
			}
		}
	}
}