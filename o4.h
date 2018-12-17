
/* O4.H  Header file */

/*#define HISTER*/ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#define ntrans     1
/* #define DEBUG */
#define twopi  6.28318
#define pi 3.14159
#define HEATBATH     /*     si no se define hace Metropolis */  

#if ntrans>1
#include <channel.h>
#include <process.h>
#include <misc.h>
#endif

#define L  8               /* LATTICE SIZE */

#define V  (L*L*L*L)

#define maxit 5000      /* maximo numero de iteraciones por bin */

#define n_obs  11        /* numero de observables                */



typedef struct{float a0,a1,a2,a3;} o4v;

 /* random number generator */

#define NormRAN (1.0F/( (float) RAND_MAX+300.0F))
#define  RAN() ( (float) rand() * NormRAN )


#define FNORM   (2.3283063671E-10F)
#define RANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )
#define FRANDOM (FNORM*RANDOM)


#ifdef MAIN
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
#else
extern unsigned char ip,ip1,ip2,ip3;
extern unsigned ira[];
#endif

#define randmax NormRAN
struct s_datos
{
   int itmax,           /* Numero de medidas por bloque                  */
       mfresh,           /* frecuencia de las medidas                    */
       nbin,            /* numero de bloque                              */
       itcut,           /* proximo bloque a calcular                     */
       flag,            /* conf de partida: 0(random), 1(fria),2(backup) */
       seed;            /* semilla random                                */
   float beta1,         /* acoplamiento a primeros vecinos               */
         beta2         /* acoplamiento a segundos vecinos                */ 
#ifdef HISTER
        ,dbeta1,        /* variacion de beta1 para la histeresis         */
         dbeta2         /* variacion de beta2 para la histeresis         */
#endif
             ;   
};


#define NDATINT   6     /* numero de campos int   en s_datos */
#ifdef HISTER
#define NDATFLOAT 4     /* numero de campos float en s_datos */
#else
#define NDATFLOAT 2     /* numero de campos float en s_datos */
#endif

#define Normaener   ( (float) (1.0/(double) V )  ) 
#define Normaener1  ( (float) (1.0/((double) V * 8)) )
#define Normaener2  ( (float) (1.0/((double) V * 24)) )

#define NOBS_HISTER 5
#define LPATH 100

FILE *Finput,*Foutput,*Fconfig;








