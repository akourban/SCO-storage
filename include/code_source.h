#ifndef CODE_SOURCE_H
#define CODE_SOURCE_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NL 3
#define NC 3
#define NSITES NL*NC

#define nhsc_max (NC*2+(NL-2)*2)
#define NBRE_D 6
#define S_MAX   (NL*(NC-1)+NC*(NL-1))

#define NPAS_MC 100000

extern int nhs,s,nhsc;
extern int sizeC, smax, size2D;

extern double d[NSITES+1][S_MAX+1][nhsc_max+1];
extern int tab[NL+2][NC+2], lig[NSITES+1], col[NSITES+1];




void init_config();
void init_densite();
void ligcol();
void montecarlo();
void imprimer();
int size_2D(int n1, int n2, int sizeC);
int index_2D(int m, int s, int c, int n, int sizeC, int smax, int size2D);


#endif