#ifndef COMM_TIME_H
#define COMM_TIME_H

#include "const.h"

typedef struct GpsTime_t
{
    time_t time;		// time (s) expressed by standard time_t
    double sec;		    // fraction of second under 1 s

    double operator - (const GpsTime_t& t) const {
        double t1 = (double)this->time + this->sec;
        double t2 = (double)t.time + t.sec;
        return (double)(t1 - t2);
    }
    bool operator < (const GpsTime_t& t) const {
        double t1 = (double)this->time + this->sec;
        double t2 = (double)t.time + t.sec;
        return (t1 < t2);
    }
    bool operator > (const GpsTime_t& t) const {
        double t1 = (double)this->time + this->sec;
        double t2 = (double)t.time + t.sec;
        return (t1 > t2);
    }
    bool operator <= (const GpsTime_t& t) const {
        double t1 = (double)this->time + this->sec;
        double t2 = (double)t.time + t.sec;
        return (t1 <= t2);
    }

    bool operator >= (const GpsTime_t& t) const {
        double t1 = (double)this->time + this->sec;
        double t2 = (double)t.time + t.sec;
        return (t1 >= t2);
    }

    bool operator == (const GpsTime_t& t) const {
        double t1 = (double)this->time + this->sec;
        double t2 = (double)t.time + t.sec;
        return (fabs(t1 - t2) < 1e-5);
    }
}Gtime;

extern double timediff(IN Gtime t1, IN Gtime t2);
extern void time2epoch(IN Gtime t, OUT int* ep);
extern Gtime epoch2time(IN int* ep);





#endif 