#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "edgetree.h"
#include "changestat.h"

#define ESPUTP 0
#define ESPOTP 1
#define ESPITP 2
#define ESPRTP 3
#define ESPOSP 4
#define ESPISP 5

/*ESP calculation functions*/
void espUTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espOTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espITP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espRTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espOSP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espISP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);

/*Changescore functions*/
D_CHANGESTAT_FN(d_desp);
D_CHANGESTAT_FN(d_dgwesp);
D_CHANGESTAT_FN(d_qwesp);
#endif
