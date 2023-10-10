#ifndef CONST_PARA_H
#define CONST_PARA_H

#include <string>
#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <numeric>
#include <math.h>
#include <string.h>
//#include "com_interface/interface_atmos.h"

using namespace std;

#ifndef IN
#define IN
#endif
#ifndef OUT
#define OUT
#endif
#ifndef IO
#define IO
#endif

/* math ---------------------------------------------------------------------------*/
#define PI			(4.0*atan(1.0))		    // phi
#define D2R         (PI/180.0)              // deg to rad
#define R2D         (180.0/PI)              // rad to deg
#define MIN(a,b)    (((a)>(b))?(b):(a))     // select smaller
#define MAX(a,b)    (((a)<(b))?(b):(a))     // select bigger

/* sat system ---------------------------------------------------------------------*/
#define MINPRNGPS   1                       // min satellite PRN number of GPS
#define MAXPRNGPS   32                      // max satellite PRN number of GPS
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) // number of GPS satellites
#define NSYSGPS     1

#define MINPRNGLO   1                       // min satellite slot number of GLONASS
#define MAXPRNGLO   27                      // max satellite slot number of GLONASS
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) // number of GLONASS satellites
#define NSYSGLO     1

#define MINPRNGAL   1                       // min satellite PRN number of Galileo
#define MAXPRNGAL   36                      // max satellite PRN number of Galileo
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1)  // number of Galileo satellites
#define NSYSGAL     1

#define MINPRNBDS   1                       // min satellite sat number of BeiDou
#define MAXPRNBDS   64                      // max satellite sat number of BeiDou
#define NSATBDS     (MAXPRNBDS-MINPRNBDS+1) // number of BeiDou satellites
#define NSATBD2     18						// number of BDS2 satellites
#define NSATBD3     (MAXPRNBDS-NSATBD2)		// number of BDS3 satellites
#define NSYSBD2     1
#define NSYSBD3     1

#define NUMSYS      (NSYSGPS+NSYSGLO+NSYSGAL+NSYSBD2+NSYSBD3) // number of all system
#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATBDS)         // number of all satellites

#define SYS_NONE    0x00                    // navigation system: none
#define SYS_GPS     0x01                    // navigation system: GPS
#define SYS_GLO     0x02                    // navigation system: GLONASS
#define SYS_GAL     0x04                    // navigation system: Galileo
#define SYS_BDS2    0x08                    // navigation system: BeiDou-2
#define SYS_BDS3    0x10                    // navigation system: BeiDou-3
#define SYS_BDS     0x18                    // navigation system: BeiDou
#define SYS_ALL     0xFF                    // navigation system: all

#define IDX_GPS     NSYSGPS-1                                     // index: GPS
#define IDX_GLO     NSYSGPS+NSYSGLO-1                             // index: GLO
#define IDX_GAL     NSYSGPS+NSYSGLO+NSYSGAL-1                     // index: GAL
#define IDX_BDS2    NSYSGPS+NSYSGLO+NSYSGAL+NSYSBD2-1             // index: BDS2
#define IDX_BDS3    NSYSGPS+NSYSGLO+NSYSGAL+NSYSBD2+NSYSBD3-1     // index: BDS3

#define INTERFACE_MAX_GPS       15          // max obs of one site: GPS
#define INTERFACE_MAX_GLO       15          // max obs of one site: GLO
#define INTERFACE_MAX_GAL       15          // max obs of one site: GAL
#define INTERFACE_MAX_BDS2      15          // max obs of one site: BDS2
#define INTERFACE_MAX_BDS3      30          // max obs of one site: BDS3

/* freq system --------------------------------------------------------------------*/
#define FREQ1_GPS   1.57542E9      /* GPS L1 frequency (Hz) */
#define FREQ2_GPS   1.22760E9      /* GPS L2 frequency (Hz) */
#define FREQ5_GPS   1.17645E9      /* GPS L5 frequency (Hz) */

#define FREQ1_GLO   1.60200E9      /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6      /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9      /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6      /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9     /* GLONASS G3 frequency (Hz) */

#define FREQ1_GAL   1.57542E9      /* Galileo E1  frequency (Hz) */
#define FREQ2_GAL   1.17645E9      /* Galileo E5a frequency (Hz) */
#define FREQ3_GAL   1.20714E9      /* Galileo E5b frequency (Hz) */

#define FREQ1_BDS   1.561098E9     /* BDS B1-2 frequency (Hz) */
#define FREQ2_BDS   1.20714E9      /* BDS B2b  frequency (Hz) */
#define FREQ3_BDS   1.26852E9      /* BDS B3   frequency (Hz) */

#define NFREQ       3                   // number of carrier frequencies

const double FREQ_ALL[NUMSYS][NFREQ] = {{FREQ1_GPS,FREQ2_GPS,FREQ5_GPS},
										{FREQ1_GLO,FREQ2_GLO,FREQ3_GLO},
										{FREQ1_GAL,FREQ2_GAL,FREQ3_GAL},
										{FREQ1_BDS,FREQ2_BDS,FREQ3_BDS},
										{FREQ1_BDS,FREQ2_BDS,FREQ3_BDS},};

/* earth para ---------------------------------------------------------------------*/
#define OMGE        7.2921151467E-5         // earth angular velocity (IS-GPS) (rad/s)
#define RE_WGS84    6378137.0               // earth semimajor axis (WGS84) (m)
#define FE_WGS84    (1.0/298.257223563)     // earth flattening (WGS84)
#define FINV_WGS84  298.257223563           // inverse of ellipsoid flatten factor of WGS84
#define RE	        6371395.0               // the average radius of earth,for compute the IPP (m)
#define H_ION       450000.0                // the height of ionospheric thin layer (m)

/* file read/write ----------------------------------------------------------------*/
#define MAXCHARS    1024                    // max characters in one line
#define MAXCOMMENT  128                     // max characters of comment
#define MAXNAME     8                       // max characters of name
#define MAXOUTCHARS 4096                    // max characters of output

/* stec model ---------------------------------------------------------------------*/
#define STECNX               3                 // stec coffeicient number
#define MAX_GRID             200               // max number of grids
#define MAX_SITE             250               // max number of sites
#define MAX_STEC_OBS         MAX_GRID+MAX_SITE // max number of stec obs
#define MAX_EPOCH_STORE      10                // max atmo info epoch
#define MAX_ITER             10                // max iteration of LSQ
#define THRES_FIXSOL_PCT     0.5			   // threshold of fix station percentage
#define THRES_USESTA_PCT     0.5			   // threshold of use station percentage
#define CUT_STEC_RES         0.5/0.16          // cut stec residual
#define CUT_DIST             3.0E5             // cut distance from site to grid
#define ERROR_VALUE          9999.999          // stec error value


//const double CUT_DIST     = 300.0 * 1E+3;
const double GRID_STEPN   = 20;
//const double ERROR_VALUE  = 9999.999;
//const double CUT_STEC_RES = (0.5 / 0.16);
const double CUT_TECU     = (300.0 / 0.16);
const double CUT_STD_RES  = 0.3;
const double CUT_COV_STEC[] = { 0.25, 0.25, 0.25, 0.25 ,0.25 };

const int MAXSTA = 250;

const int DAY_PER_MON[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
const double TIME_DIFF_LIMIT = 0.25;

#endif // !CONST_PARA_H



