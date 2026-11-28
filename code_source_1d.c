/* 28 fev 2017 entropic sampling version modif1:
"s" garder comme (s+S_MAX)/2 for 2D rectangular NLxNC avec NC>=NL
(nhs,s,nhsc) en travaillant avec un réseau (NL+2)*(NC+2) */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NL 3
#define NC 3
#define NSITES NL*NC

#define nhsc_max ( NC*2+(NL-2)*2 )
#define NBRE_D 6
#define S_MAX   ( NL*(NC-1)+NC*(NL-1) )

#define NPAS_MC 100000

int nhs,s,nhsc;
int sizeC, smax, size2D;

float d[NSITES+1][S_MAX+1][nhsc_max+1];
int tab[NL+2][NC+2], lig[NSITES+1], col[NSITES+1];
void init_config (void);
void init_densite(void);
void ligcol (void);
void montecarlo (void) ;
void imprimer (void) ;

int size_2D(int n1, int n2, int sizeC){
  int n, smax, ll, rs;
  n = n1*n2;
  if (sizeC<n2+2*n1) sizeC = n2+2*n1; // default value of sizeC for algorithm
  smax = n1*(n2-1)+n2*(n1-1);
  rs=smax%4;
  if (rs!=0) smax=smax+(4-rs);//necessary for even number n
  
  if (n%2 == 0){
    ll = n/2-1;
    return 4*ll*(ll+1)*sizeC + smax*sizeC + (n+1)*sizeC ;
  }
  else {
    ll = (n-1)/2-1; 
    return 4*ll*(ll+1)*sizeC + 2*smax*sizeC + (n+1)*sizeC;
  }
}

int index_2D(int m, int s, int c, int n, int sizeC, int smax, int size2D){
  int k;
  int sizeP1, ll;
  int ii,kk;
  
  if (s==0){
    k = size2D - (n+1)*sizeC + m*sizeC + c ;
  }
  else{
    
    if (n%2 == 0){
      ll = n/2-1;
      sizeP1 = 2*ll*(ll+1)*sizeC;
      if (m < (n/2)+1){
	if ((s<=4*m)&&(s<=smax)){ 
	  k = 2*(m-1)*((m-1)+1)*sizeC + (s-1)*sizeC + c;
	}
	else{
	  k = -1;
	}
      }
      else {
	ii = (m-1-n/2);
	kk = 2*(ll*(ll+1)-(ll-ii)*(ll-ii+1));
	if (s<=4*(n-m)){
	  k = sizeP1 + (smax)*sizeC + kk*sizeC + (s-1)*sizeC + c;
        }
	else{
	  k=-1;
	}
      }
    }
    else {
      ll = (n-1)/2-1;
      sizeP1 = 2*ll*(ll+1)*sizeC;
      sizeP1 += smax*sizeC;
      if (m <= (n-1)/2){
	if ((s<=4*m)&&(s<=smax)){
	  k = 2*(m-1)*((m-1)+1)*sizeC + (s-1)*sizeC + c;
	}
	else{
	  k=-1;
	}
      }
      else if (m == (n-1)/2+1) {
	if (s<=smax){
	  k = sizeP1 + (s-1)*sizeC + c;
	}
	else{
	  k=-1;
	}
      }
      else {
	ii = (m-2-(n-1)/2);
	kk = 2*(ll*(ll+1)-(ll-ii)*(ll-ii+1));
	if (s<=4*(n-m)){
	  k = sizeP1 + smax*sizeC + kk*sizeC + (s-1)*sizeC + c;
	}
	else{
	  k=-1;
	}
      }
    }
    //if (k>=(size2D-(n+1)*sizeC)) fprintf(stderr, "ERROR out of triangle %d > %d \n",k, size2D-(n+1)*sizeC);
  }
  //if (k==-1) fprintf(stderr, "ERROR out of triangle %d = (%d,%d,%d)  \n",k,m,s,c);
  return k;
}

void init_config (){
  int i,j,s_v;
  nhs=0;
  s=0;
  nhsc=0;

  for (i=0;i<=NL+1;i++)
    for (j=0;j<=NC+1;j++)
      tab[i][j]=0;

  for (i=1;i<=NL;i++)
    for (j=1;j<=NC;j++){
      tab[i][j]=-1;

      if (((float)rand()/(RAND_MAX+1.)) >= 0.5){
        tab[i][j]=1;
        nhs++;
      }
    }
    for (i=1;i<=NL;i++)
      for (j=1;j<=NC;j++){
        s_v=tab[i+1][j]+tab[i-1][j]+tab[i][j+1]+tab[i][j-1];
        s=s+tab[i][j]*s_v ;
      }
      
      s=s/2;

    for (j=1;j<=NC;j++){
      if (tab[1][j] == 1) nhsc ++;
      if (tab[NL][j] == 1) nhsc++;  
    }

    for (i=2;i<=NL-1;i++){
      if (tab[i][1]==1) nhsc++;
      if (tab[i][NC] ==1) nhsc++;
    }

    for (i=0;i<=NL+1;i++){
      printf (" \n ");
      for (j=0;j<=NC+1;j++)
        printf( "%3d ", tab[i][j]);
    }
 }

 void init_densite(){
   int i,j,k;
   for (i=0;i<=NSITES;i++)
      for (j=0;j<=S_MAX;j++)
          for (k=0; k<=nhsc_max;k++)
             d[i][j][k]=1.;
}

 void ligcol(){
  int i,j,k,ii;
  ii=1;
  for (i=1;i<=NL;i++)
    for (j=1;j<=NC;j++){
      lig[ii]=i;
      col[ii]=j;
      ii++;
    }
}

void imprimer(){
  int i,j,k,n1=0,aux;
  FILE *psortie;

  psortie=fopen("r6_6_mod.dat","w+");

  if (psortie==NULL){
    printf("erreur à l'ouverture");
    exit(1);
  }

  for (i=0;i<=NSITES;i++)
    for (j=0;j<=S_MAX;j++)
      for (k=0; k<=nhsc_max;k++)
        if (d[i][j][k] !=1.) n1++;
  fprintf(psortie,"%d\n",n1);
  printf(" \n nombre de lignes de densité %d\n", n1);


  for (i=0;i<=NSITES;i++)
    for (j=0;j<=S_MAX;j++)
      for (k=0; k<=nhsc_max;k++)
        if (d[i][j][k] !=1.){
          fprintf(psortie," %d %d %d %d %d %d %15.3f\n ",i,j*2-S_MAX, k, 2*i-NSITES, j*2-S_MAX,2*k-nhsc_max,d[i][j][k]);
        }

    printf("%d %d %15.3f\n ",NSITES,S_MAX,d[NSITES][S_MAX][nhsc_max]);
    fclose(psortie);
}

 void montecarlo (){
  int icount=1,i,j,k,is_v,snew,nhsnew,nhscnew,nbre;
  int nrun, irun;
  float r,proba;
  float h[NSITES+1][S_MAX+1][nhsc_max+1];
  double moycarr,carmoy,ety,etyrel,som_h,scar_h,nval,moy,som_d;

  printf("avant while nhs, s, nhsc \n %d %d %d ",nhs,s,nhsc);

  while (icount<= NBRE_D){
    /* RAZ de l'histogramme h */
    for (i=0;i<=NSITES;i++){
      for (j=0;j<=S_MAX;j++){
        for (k=0; k<=nhsc_max;k++){
          h[i][j][k]=1.;
        }
      }
    }

     /* NPAS_MC*NSITES pas pour construire h[][] en utilisant d[][] */
      /*if (icount==1)
         nrun=NPAS_MC*NSITES*10;
        else
        */
    nrun=NPAS_MC;

    for (irun=1;irun<=nrun;irun++){
      nbre=1+rand()/(RAND_MAX+1.) *(NSITES-1+1);/* nbre au hazard entre 0 et NSITES-1 */
      i=lig[nbre];
      j=col[nbre];
      /*printf (" nbre: %d  i: %d  j: %d \n",nbre, i, j);  */


      is_v=tab[i+1][j]+tab[i-1][j]+tab[i][j+1]+tab[i][j-1];
      snew=s-2*is_v*tab[i][j];
      nhsnew=nhs-tab[i][j];
      nhscnew=nhsc;
      if (i==1 || i==NL){
        nhscnew=nhsc-tab[i][j] ;
      }else if (j==1 || j==NC){
        nhscnew=nhsc-tab[i][j] ;
      }

      /*printf (" nhs s nhsc  %d %d %d \n",nhs, s, nhsc);
        printf (" nhsnew snew nhscnew  %d %d %d \n",nhsnew, snew, nhscnew);
        getchar(); */
      if (d[nhs][(s+S_MAX)/2][nhsc] >= d[nhsnew][(snew+S_MAX)/2][nhscnew]){
        tab[i][j]=-tab[i][j];
        nhs=nhsnew;
        s=snew;
        nhsc=nhscnew;
      }else{
        proba=d[nhs][(s+S_MAX)/2][nhsc]/d[nhsnew][(snew+S_MAX)/2][nhscnew] ;
        r=(float)rand()/(RAND_MAX+1.) ;
        if ( r < proba ){
          tab[i][j]=-tab[i][j];
          nhs=nhsnew;
          s=snew;
          nhsc=nhscnew;
        }
      }
      h[nhs][(s+S_MAX)/2][nhsc] +=1;
    }

         /* on regarde h */
    som_h=0.;
    scar_h=0.;
    nval=0.;
    for (i=0;i<=NSITES;i++)
      for (j=0;j<=S_MAX;j++)
        for (k=0; k<=nhsc_max;k++)
          if (h[i][j][k]!=1){
            nval++;
            som_h+=h[i][j][k];
            scar_h+= h[i][j][k]*h[i][j][k];
          }

    moy=som_h/nval;
    moycarr=scar_h/nval;
    carmoy=moy*moy;
    ety=moycarr-carmoy;
    etyrel=sqrt(fabs(ety))/moy;
    printf("\n etape numero %d ecart type relatif=%lf",icount,etyrel);
    som_d=0;

    /* on améliore d */
    for (i=0;i<=NSITES;i++)
      for (j=0;j<=S_MAX;j++)
        for (k=0; k<=nhsc_max;k++){
          d[i][j][k]=d[i][j][k]*h[i][j][k];
          if (d[i][j][k]!=1)
            som_d+=d[i][j][k];
        }

    for (i=0;i<=NSITES;i++)
      for (j=0;j<=S_MAX;j++)
        for (k=0; k<=nhsc_max;k++)
          if (d[i][j][k] !=1){
              d[i][j][k]=d[i][j][k]/som_d;
              d[i][j][k] = d[i][j][k]* pow(2,NSITES);
          }

    icount+=1;
  }
}

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

