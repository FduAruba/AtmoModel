#ifndef COMM_FUNCTION_H
#define COMM_FUNCTION_H

#include "const.h"
#include "commtime.h"

extern void modifyPath(IO char* src);
extern void trimSpace(IO char* src);
extern void cutSep(IO char* str);

extern void pos2ecef(IN double* pos, OUT double* r);

extern bool checksys(IN char s);
extern string sysstr(IN int sys);
extern char syschar(IN string system);
extern char idx2sys(IN int idx);
extern int satsys(IN int sat, OUT int* prn);
extern int satno(IN int sys, IN int prn);
extern int satid2no(IN char* id);
extern int satid2idx(IN char code, IN int num);

extern string strtime(IN Gtime t, IN int opt);
extern void time2str(IN Gtime t, IN char* s, IN int opt);
extern double str2num(IN char* s, IN int i, IN int n);

#endif
