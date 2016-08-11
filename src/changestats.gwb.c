/*  
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 */
#include "changestats.gwb.h"

D_CHANGESTAT_FN(d_gwb1nsp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    STEP_THROUGH_INEDGES(head, e, u){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u 
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        } */
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    
    STEP_THROUGH_INEDGES(tail, e, u){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u 
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }*/
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){
    cumchange=0.0;
    L2th=0;
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;

    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_UNDIRECTED_EDGE(u, tail)){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }
    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2th)) ;
    }else{
      cumchange += (double)L2th;
    }
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) -= cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

}

D_CHANGESTAT_FN(d_gwb2nsp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (u != tail){
        L2tu=ochange;
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    STEP_THROUGH_OUTEDGES(tail, e, u){
      if (u != head){
        L2uh=ochange;
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }

    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){
    cumchange=0.0;
    L2th=0;
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_UNDIRECTED_EDGE(u, tail)){
  L2th++;
	L2tu=ochange;
	L2uh=ochange;

	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }

    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2th)) ;
    }else{
      cumchange += (double)L2th;
    }
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) -= cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

}
