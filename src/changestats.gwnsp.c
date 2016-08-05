#include "changestats.gwnsp.h"


/*****************
 changestat: d_nsp
*****************/
/*
Note that d_esp is a meta-function, dispatching actual changescore 
calculation to one of the esp*_calc routines, based on the selected shared 
partner type code.

Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Recursive two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
D_CHANGESTAT_FN(d_dnsp) { 
  int i,type;
  double *dvec,*cs_esp, *cs_dsp;
  
  /*Set things up*/
  ZERO_ALL_CHANGESTATS(i);
  type=(int)INPUT_PARAM[0];     /*Get the ESP type code to be used*/
  dvec=INPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/
  cs_esp=Calloc(N_CHANGE_STATS,double);     /*Allocate memory for the DSP changescores*/
  cs_dsp=Calloc(N_CHANGE_STATS,double);     /*Allocate memory for the DSP changescores*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: 
    espUTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspUTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPOTP: 
    espOTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspOTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPITP: 
    espITP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspITP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPRTP: 
    espRTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspRTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPOSP: 
    espOSP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspOSP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPISP: 
    espISP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspISP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
  
  for(i=0;i<N_CHANGE_STATS;i++)
    CHANGE_STAT[i]=(cs_dsp[i]-cs_esp[i]);
    
  /*Free the changescore memory*/
  Free(cs_esp);
  Free(cs_dsp);
}


/*****************
 changestat: d_gwnsp
*****************/

/*
Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Recursive two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
D_CHANGESTAT_FN(d_dgwnsp) { 
  int type;
  Vertex i,maxesp;
  double alpha, oneexpa,*dvec,*cs_esp, *cs_dsp;
  
  /*Set things up*/
  CHANGE_STAT[0] = 0.0;         /*Zero the changestat*/
  alpha = INPUT_PARAM[0];       /*Get alpha*/
  oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  type=(int)INPUT_PARAM[1];     /*Get the ESP type code to be used*/
  maxesp=(int)INPUT_PARAM[2];   /*Get the max ESP cutoff to use*/
  cs_esp=Calloc(maxesp,double);     /*Allocate memory for the ESP changescores*/
  dvec=Calloc(maxesp,double);   /*Allocate memory for the ESP vals*/
  for(i=0;i<maxesp;i++)         /*Initialize the ESP vals*/
    dvec[i]=i+1.0;
  
  cs_dsp=Calloc(maxesp,double);     /*Allocate memory for the ESP changescores*/

  /*Obtain the changescores (by type)*/
  switch(type){
    case ESPUTP: 
      espUTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspUTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPOTP: 
      espOTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspOTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPITP: 
      espITP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspITP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPRTP: 
      espRTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspRTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPOSP: 
      espOSP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspOSP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPISP: 
      espISP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspISP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
  }
  
  /*Compute the gwnsp changescore*/
  for(i=0;i<maxesp;i++)
    if((cs_dsp[i]-cs_esp[i])!=0.0)                  /*Filtering to save a few pow() calls*/
      CHANGE_STAT[0]+=(1.0-pow(oneexpa,dvec[i]))*(cs_dsp[i]-cs_esp[i]);
  CHANGE_STAT[0]*=exp(alpha);

  /*Free the changescore memory*/
  Free(cs_esp);
  Free(cs_dsp);
  Free(dvec);
}




