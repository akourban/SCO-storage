#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "code_source.h"

int main(){
  int i,j,k;

  sizeC = nhsc_max;
  if (sizeC < NC + 2 * NL) sizeC = NC + 2 * NL;

  smax = S_MAX;
  if (smax % 4!=0) smax += 4 -(smax % 4);

  size2D = size_2D(NL,NC, sizeC);

 	init_config();
  init_densite();
  ligcol();

  montecarlo();
  //imprimer();
  printf ( "FIN \n");
  //getchar();

  return 0;
}