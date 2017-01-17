
#include "changestats.users.h"


CHANGESTAT_FN(d_gwdspattr){
  Edge e, f;
  int i, echange, ochange;
  int L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange, p1attr, p2attr, tattr, hattr, uattr;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[N_NODES+2];
  oneexpa = 1.0-exp(-alpha);
  
  p1attr = INPUT_PARAM[N_NODES];
  p2attr = INPUT_PARAM[N_NODES + 1];
  
  FOR_EACH_TOGGLE(i){
	tail = TAIL(i); 
	head = HEAD(i);
	tattr = INPUT_PARAM[tail-1];
	hattr = INPUT_PARAM[head-1];
		
    cumchange=0.0;
    ochange = IS_OUTEDGE(tail, head) ? -1 : 0;
    echange = 2*ochange + 1;
	
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (u != tail){
		uattr = INPUT_PARAM[u-1];
        L2tu=ochange;
        /* step through outedges of u */
		if(MIN(tattr,uattr) == MIN(p1attr,p2attr) && MAX(tattr,uattr) == MAX(p1attr,p2attr)){
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
	  }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head, e, u){
      if (u != tail){
		uattr = INPUT_PARAM[u-1];
        L2tu=ochange;
        /* step through outedges of u */
		if(MIN(tattr,uattr) == MIN(p1attr,p2attr) && MAX(tattr,uattr) == MAX(p1attr,p2attr)){
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
	  }
    }
    
    /* step through outedges of tail  */
    STEP_THROUGH_OUTEDGES(tail, e, u){
      if (u != head){
		uattr = INPUT_PARAM[u-1];
        L2uh=ochange;
        /* step through outedges of u */
		if(MIN(hattr,uattr) == MIN(p1attr,p2attr) && MAX(hattr,uattr) == MAX(p1attr,p2attr)){
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
	  }
    }
    /* step through inedges of tail */
    STEP_THROUGH_INEDGES(tail, e, u){
      if (u != head){
		uattr = INPUT_PARAM[u-1];
        L2uh=ochange;
        /* step through outedges of u */
		if(MIN(hattr,uattr) == MIN(p1attr,p2attr) && MAX(hattr,uattr) == MAX(p1attr,p2attr)){
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
	}
    }
	
    
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}
