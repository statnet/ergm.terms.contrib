
CHANGESTAT_FN(d_gwnspattrhetero){
  Edge e, f;
  int i, echange, ochange;
  int L2tu, L2uh, L2th;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange, exattr, tattr, hattr, uattr, vattr;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[N_NODES+1];
  oneexpa = 1.0-exp(-alpha);
  
  exattr = INPUT_PARAM[N_NODES];
  
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
      uattr = INPUT_PARAM[u-1];
      if (u != tail && tattr != exattr && uattr != exattr && tattr != uattr){
        L2tu=ochange;
        /* step through outedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(v != head && IS_OUTEDGE(tail, v)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    
    /* step through inedges of tail */
    STEP_THROUGH_INEDGES(tail, e, v){
      vattr = INPUT_PARAM[v-1];
      if (v != head && vattr != exattr && hattr != exattr && vattr != hattr){
        L2uh=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(v, f, u){
          if(u != tail && IS_INEDGE(head, u)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
  
  alpha = INPUT_PARAM[N_NODES+1];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){
    tail = TAIL(i); 
    head = HEAD(i);
    tattr = INPUT_PARAM[tail-1];
    hattr = INPUT_PARAM[head-1];
    cumchange=0.0;
    L2th=0;
    ochange = IS_OUTEDGE(tail, head) ? -1 : 0;
    echange = 2*ochange + 1;
    
    STEP_THROUGH_OUTEDGES(head, e, u){
      uattr = INPUT_PARAM[u-1];
      if (u != tail && tattr != exattr && uattr != exattr && tattr != uattr && IS_OUTEDGE(tail, u)){
        L2tu=ochange;
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_OUTEDGE(tail, v)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head, e, u){
      uattr = INPUT_PARAM[u-1];
      if (tattr != exattr && hattr != exattr && tattr != hattr && IS_OUTEDGE(tail, u)){
        L2th++;
      }
      if (uattr != exattr && hattr != exattr && uattr != hattr && IS_OUTEDGE(u, tail)){
        L2uh=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_OUTEDGE(v, head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh) ;
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

