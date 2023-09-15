#include "commtime.h"


extern double timediff(IN Gtime t1, IN Gtime t2)
{
	return (difftime(t1.time, t2.time) + t1.sec - t2.sec);
}

extern void time2epoch(IN Gtime t, OUT int* ep)
{
	const int mday[] = { /* # of days in a month */
		31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
		31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
	};
	int days, sec, mon, day;

	/* leap year if year%4==0 in 1901-2099 */
	days = (int)(t.time / 86400);
	sec = (int)(t.time - (time_t)days * 86400);
	for (day = days % 1461, mon = 0; mon < 48; mon++) {
		if (day >= mday[mon]) { day -= mday[mon]; }
		else { break; }
	}

	ep[0] = 1970 + days / 1461 * 4 + mon / 12;
	ep[1] = mon % 12 + 1;
	ep[2] = day + 1;
	ep[3] = sec / 3600;
	ep[4] = sec % 3600 / 60;
	ep[5] = sec % 60 + (int)t.sec;
}

extern Gtime epoch2time(IN int* ep)
{
	const int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
	Gtime gtime = { 0 };
	int days = 0;
	int year = (int)ep[0], month  = (int)ep[1], day = (int)ep[2];
	int hour = (int)ep[3], minute = (int)ep[4], sec = (int)ep[5];

	if (year < 1970 || 2099 < year || month < 1 || 12 < month) { return gtime; }

	/* leap year if year%4==0 in 1901-2099 */
	days = (year - 1970) * 365 + (year - 1969) / 4 + doy[month - 1] + day - 2 + (year % 4 == 0 && month >= 3 ? 1 : 0);
	gtime.time = days * 86400 + hour * 3600 + minute * 60 + sec;
	gtime.sec  = ep[5] - (int)floor(ep[5]);

	return gtime;
}