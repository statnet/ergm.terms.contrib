#include "changestats.maxdeg.h"

CHANGESTAT_FN(d_maxdegree) {
  Vertex t, h, node3;
  int i, maxdeg, hdeg, tdeg;
  Edge e;
  int attrflag;
  double t_nodecov, h_nodecov;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    t = TAIL(i); h = HEAD(i);
    attrflag = INPUT_PARAM[0];
    maxdeg = INPUT_PARAM[1];
    if(attrflag==0){
      tdeg = IN_DEG[t]+OUT_DEG[t];
      hdeg = IN_DEG[h]+OUT_DEG[h];
      CHANGE_STAT[0] += IS_UNDIRECTED_EDGE(t,h) ?
      /* edge exists and will be removed by toggle,
         so increment if we will drop down to maxdeg*/
         (tdeg==maxdeg+1) +  (hdeg==maxdeg+1) :   
      /* edge does not exist and will be added by toggle,
         so decrement if we will move above maxdeg*/  
         - (tdeg==maxdeg) - (hdeg==maxdeg); 
    }else{
      t_nodecov = INPUT_PARAM[t+1];
      h_nodecov = INPUT_PARAM[h+1];
      if (t_nodecov == h_nodecov) {
        tdeg = 0;
        STEP_THROUGH_OUTEDGES(t, e, node3) { /* step through outedges of tail */
          if(INPUT_PARAM[node3+1]==t_nodecov){++tdeg;}
        }
        STEP_THROUGH_INEDGES(t, e, node3) { /* step through inedges of tail */
          if(INPUT_PARAM[node3+1]==t_nodecov){++tdeg;}
        }
        hdeg = 0;
        STEP_THROUGH_OUTEDGES(h, e, node3) { /* step through outedges of head */
          if(INPUT_PARAM[node3+1]==h_nodecov){++hdeg;}
        }
        STEP_THROUGH_INEDGES(h, e, node3) { /* step through inedges of head */
          if(INPUT_PARAM[node3+1]==h_nodecov){++hdeg;}
        }
        CHANGE_STAT[0] += IS_UNDIRECTED_EDGE(t,h) ?
           (tdeg==maxdeg+1) + (hdeg==maxdeg+1) :
           -(tdeg==maxdeg) - (hdeg==maxdeg);
      }else{
        CHANGE_STAT[0] = 0;
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

