#ifndef OUTFILES_H
#define OUTFILES_H

#include "readfiles.h"

using namespace std;

void outRovStec(IN Coption& cfg, IN GridInfo& grid, IN SiteAtmos& rovaug, IN ProStecMod& stecmod, IN FileFps& rovfps, IN int type);
void createRovFile(IN Coption& cfg, OUT FileFps& fps);


#endif
