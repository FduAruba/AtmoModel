#include "commfun.h"

/* string */
extern void modifyPath(IO char* src)
{
	trimSpace(src);
	cutSep(src);
}

extern void trimSpace(IO char* src)
{
	int i = 0, ps = 0, pe = 0;
	int len = (int)strlen(src);
	char str[MAXCHARS + 1] = { '\0' };

	if (len <= 0) { return; }
	strcpy(str, src);

	for (i = 0, ps = 0; i < len; i++) {
		if (str[i] != ' ' && str[i] != '\t') {
			ps = i;
			break;
		}
	}

	for (i = len - 1, pe = ps; i >= 0; i--) {
		if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n') {
			pe = i;
			break;
		}
	}

	i = pe == ps ? 0 : 1;
	str[pe + i] = '\0';

	strcpy(src, str + ps);
}

extern void cutSep(IO char* str)
{
	int i, len;

	for (i = 0; i < 4; i++) {
		len = (int)strlen(str);

		if (len <= 0) { break; }

		if (str[len - 1] == '\\') {
			str[len - 1] = '\0';
		}
		else { break; }
	}
}

extern double str2num(IN char* s, IN int i, IN int n)
{
	double value;
	char str[256];
	char* p = str;

	if ((i < 0) || ((int)strlen(s) < i) || ((int)sizeof(str) - 1 < n)) {
		return 0.0;
	}
	for (s += i; *s && --n >= 0; s++) {
		*p++ = (*s == 'd' || *s == 'D') ? 'E' : *s;
	}
	*p = '\0';
	if (sscanf_s(str, "%lf", &value) == 1) { return value; }
	else { return 0.0; }
}

/* position */
extern void pos2ecef(IN double* pos, OUT double* r)
{
	double sinp = sin(pos[0]), cosp = cos(pos[0]), sinl = sin(pos[1]), cosl = cos(pos[1]);
	double e2 = FE_WGS84 * (2.0 - FE_WGS84), v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);

	r[0] = (v + pos[2]) * cosp * cosl;
	r[1] = (v + pos[2]) * cosp * sinl;
	r[2] = (v * (1.0 - e2) + pos[2]) * sinp;
}

extern double sphereDist(IN double latG, IN double lonG, IN double latB, IN double lonB)
{
	double m;
	m = sin(latB) * sin(latG) + cos(latB) * cos(latG) * cos((lonB - lonG));
	if (fabs(m) >= 1.0) {
		return 0.0;
	}
	else {
		return RE_WGS84 * acos(m);
	}
}

extern bool gridVaild(IN double latG, IN double lonG, IN double lat, IN double lon, IN double dlat, IN double dlon)
{
	return (latG > lat - dlat) && (latG < lat + dlat) && (lonG > lon - dlon) && (lonG < lon + dlon);
}

/* sat system */
extern bool checksys(IN char s)
{
	char sys[] = "GREC";
	int len = (int)strlen(sys);

	for (int i = 0; i < len; i++) {
		if (s == sys[i]) {
			return true;
		}
	}

	return false;
}

extern string sysstr(IN int sys)
{
	string system = "NONE";

	switch (sys)
	{
	case SYS_GPS:  {system = "GPS";  break; }
	case SYS_GLO:  {system = "GLO";  break; }
	case SYS_GAL:  {system = "GAL";  break; }
	case SYS_BDS2: {system = "BDS2"; break; }
	case SYS_BDS3: {system = "BDS3"; break; }
	default:       {system = "NONE"; break; }
	}

	return system;
}

extern char syschar(IN string system)
{
	if		(system == "GPS") { return 'G'; }
	else if (system == "GLO") { return 'R'; }
	else if (system == "GAL") { return 'E'; }
	else if (system == "BDS2" || system == "BDS3") 
							  { return 'C'; }

	return '\0';
}

extern char idx2sys(IN int idx)
{
	switch (idx)
	{
	case 0:  { return 'G'; }
	case 1:  { return 'R'; }
	case 2:  { return 'E'; }
	case 3:  { return 'C'; }
	case 4:  { return 'C'; }
	default: { return 'X'; }
	}
	return 'X';
}

extern int satsys(IN int sat, OUT int* prn)
{
	int sys = SYS_NONE;

	if (sat <= 0 || MAXSAT < sat) {
		sys = SYS_NONE; sat = 0;
	}
	else if (sat <= NSATGPS) {
		sys = SYS_GPS; sat += MINPRNGPS - 1;
	}
	else if ((sat -= NSATGPS) <= NSATGLO) {
		sys = SYS_GLO; sat += MINPRNGLO - 1;
	}
	else if ((sat -= NSATGLO) <= NSATGAL) {
		sys = SYS_GAL; sat += MINPRNGAL - 1;
	}
	else if ((sat -= NSATGAL) <= NSATBDS) {
		sys = sat <= NSATBD2 ? SYS_BDS2 : SYS_BDS3;		
		sat += MINPRNBDS - 1;
	}
	else { 
		sys = SYS_NONE; sat = 0; 
	}

	if (prn) { *prn = sat; }

	return sys;
}

extern int satno(IN int sys, IN int prn)
{
	if (prn <= 0) { return 0; }

	switch (sys)
	{
	case SYS_GPS: {
		if (prn < MINPRNGPS || MAXPRNGPS < prn) { return 0; }
		return prn - MINPRNGPS + 1;
	}
	case SYS_GLO: {
		if (prn < MINPRNGLO || MAXPRNGLO < prn) { return 0; }
		return NSATGPS + prn - MINPRNGLO + 1;
	}	
	case SYS_GAL: {
		if (prn < MINPRNGAL || MAXPRNGAL < prn) { return 0; }
		return NSATGPS + NSATGLO + prn - MINPRNGAL + 1;
	}
	case SYS_BDS: {
		if (prn < MINPRNBDS || MAXPRNBDS < prn) { return 0; }
		return NSATGPS + NSATGLO + NSATGAL + prn - MINPRNBDS + 1;
	}
	}

	return 0;
}

extern int satid2no(IN char* id)
{
	int sys, prn;
	char code;

	if (sscanf(id, "%c%d", &code, &prn) < 2) { return 0; }

	switch (code)
	{
	case 'G': {sys = SYS_GPS; prn += MINPRNGPS - 1; break; }
	case 'R': {sys = SYS_GLO; prn += MINPRNGLO - 1; break; }
	case 'E': {sys = SYS_GAL; prn += MINPRNGAL - 1; break; }
	case 'C': {sys = SYS_BDS; prn += MINPRNBDS - 1; break; }
	default: { return 0; }
	}

	return satno(sys, prn);
}

extern int satid2idx(IN char code, IN int num) {
	int idx = -1;

	switch (code)
	{
	case 'G': { 
		if (MINPRNGPS <= num && num <= MAXPRNGPS) {
			idx = IDX_GPS;
		}
		break; 
	}
	case 'R': {
		if (MINPRNGLO <= num && num <= MAXPRNGLO) {
			idx = IDX_GLO;
			break;
		}
	}
	case 'E': {
		if (MINPRNGAL <= num && num <= MAXPRNGAL) {
			idx = IDX_GAL;
			break;
		}
	}
	case 'C': {
		if (MINPRNBDS <= num && num <= MAXPRNBDS) {
			idx = num <= NSATBD2 ? IDX_BDS2 : IDX_BDS3;
			break;
		}
	}
	default: {break;}
	}

	return idx;
}

extern double satfreq(IN int sys, IN int f)
{
	return FREQ_ALL[sys][f];
}

/* gps time */
extern string strtime(IN Gtime t, IN int opt)
{
	char buff[64] = { '\0' };
	string tstr;
	time2str(t, buff, opt);
	tstr = buff;
	return tstr;
}

extern void time2str(IN Gtime t, IN char* s, IN int opt)
{
	int ep[6];
	time2epoch(t, ep);
	switch (opt) 
	{
	case 1: {
		sprintf(s, "%04d/%02d/%02d %02d:%02d:%02d\0", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
		break;
	}
	case 2: {
		sprintf(s, "%02d:%02d:%02d\0", ep[3], ep[4], ep[5]);
		break;
	}
	default: {
		break;
	}
	}
}

/* math */
extern void calcMeanStd(IN vector<double> data, OUT double& vmean, OUT double& vstd)
{
	vmean = vstd = 0.0;
	
	int n = data.size();
	if (n > 0) {
		double sum = accumulate(data.begin(), data.end(), 0.0);
		vmean = sum / n;

		double var = 0.0;
		for (int i = 0; i < data.size(); i++) {
			var += pow(data[i] - vmean, 2);
		}
		vstd = sqrt(var / n);
	}
}

extern double robust(IN double V, OUT double rms)
{
	double k0 = 1.5;
	double k1 = 3.0;
	double scale = 1.0;

	double v = fabs(V / rms);
	if (v > k0 && v <= k1) {
		scale = (k0 / v) * pow((k1 - v) / (k1 - k0), 2);
	}
	else if (v > k1) {
		scale = 0;
	}

	return scale;
}

extern double modelIDW(IN StaDistIonArr& list, IN int nused, IN double maxdist, IN int k, OUT int* n)
{
	int cnt = 0;
	double sumwgt = 0.0, stawgt[MAXSTA] = { 0.0 };
	double res = 0.0;

	for (auto& p : list) {
		if (cnt >= nused || p.dist > maxdist) {
			break;
		}

		stawgt[cnt] = pow(100.0 / p.dist, k);
		sumwgt += stawgt[cnt];
		cnt++;
	}

	if (n) { *n = cnt; }

	if (cnt < 1) { return ERROR_VALUE; }

	for (int i = 0; i < cnt; i++) {
		res += stawgt[i] / sumwgt * list[i].ion;
	}

	return res;
}