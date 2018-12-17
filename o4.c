

/* MAIN PROGRAM  */

#define MAIN
#include "o4.h"
#include "o4aritm.h"


/* Declaraciones de funciones externas al main */

extern void Init_Rand(int); 
extern void Overrelaxed(int);
extern void Direccionamientos(void);
extern void Lee_Datos(void);
extern void Genera_O4(int);
extern void Neigh(int);
extern void Inicializa(int);
extern void Energia(int);
extern void Canonic_Update(int);
extern void Order(void);
extern void Mide_Ener(void);

FILE *Fheat;

o4v f[V],SUMA,suma1,suma2,cero,celda[16];

int neigh_px,neigh_py,neigh_pz,neigh_pt,
    neigh_mx,neigh_my,neigh_mz,neigh_mt;

int x_p[L],y_p[L],z_p[L],t_p[L],
    x_m[L],y_m[L],z_m[L],t_m[L];

float E_1,E_2,ener1,ener2,Mf,Mhpaf,Mafp;

float beta1,beta2,resultados[NOBS_HISTER];

int jpp[6],jpm[6],jmp[6],jmm[6],it;

int jpx,jpy,jpz,jpt,
    jmx,jmy,jmz,jmt;

int seed,flag,veces,good;

struct s_datos datos;

float obs_dat[n_obs][maxit];


int main(void)
{

/* Declaracion de variables que solo se usan en `main` */

   int x,y,z,t,site,mfr,mfresh,ibin,j;  /* contadores */

   Fheat=fopen("heatb.plot","w");
  

   Direccionamientos();         /* define las condiciones de contorno */
 

   Lee_Datos();
                
   cero.a0=cero.a1=cero.a2=cero.a3=0.0F;
   

   Init_Rand(datos.seed);


   Inicializa(datos.flag);

   for(ibin=datos.itcut;ibin<datos.nbin;ibin++)    /* numero de bloques */
   {

      Init_Rand((unsigned) datos.seed);

#ifdef HISTER

   if(ibin<datos.nbin/2)
   {
     datos.beta1 +=datos.dbeta1;
     datos.beta2 +=datos.dbeta2;
   }
   else
   {   
     datos.beta1 -=datos.dbeta1;
     datos.beta2 -=datos.dbeta2;
   }
     
#endif

    beta1=datos.beta1;
    beta2=datos.beta2;
    mfresh=datos.mfresh;
  

      for(it=0;it<datos.itmax;it++)
      {
         for(mfr=0;mfr<mfresh;mfr++)              /* loop sin medidas */
         { 
            site=0;
     
            for(t=0;t<L;t++)
            {
               neigh_pt=t_p[t];
               neigh_mt=t_m[t];
      
               for(z=0;z<L;z++)
               {
                 neigh_pz=z_p[z];
                 neigh_mz=z_m[z];
         
                 for(y=0;y<L;y++)
                 {
                     neigh_py=y_p[y];
                     neigh_my=y_m[y];

                     for(x=0;x<L;x++)
                     {
                        neigh_px=x_p[x];
                        neigh_mx=x_m[x];
                     
                        if(mfr==0) Canonic_Update(site); 
                           Overrelaxed(site);
                                  
                        site++;
                        good=0;
                     }             /* coord x */
                 }                /* coord y */
              }                  /* coord z */
            }                   /* coord t */   
         }        /* fin de mfresh */

     
                              /*  MEDIDAS  */
      Mide_Ener();

      Order();


      obs_dat[0][it] = E_1 * Normaener1/beta1;
      obs_dat[1][it] = E_2 * Normaener2/beta2;
      obs_dat[2][it] = Mf;
      obs_dat[3][it] = Mhpaf;
      obs_dat[4][it] = Mafp;
      obs_dat[5][it] = Mf*Mf;
      obs_dat[6][it] = Mhpaf*Mhpaf;
      obs_dat[7][it] = Mafp*Mafp;
      obs_dat[8][it] = Mf*Mf*Mf*Mf;
      obs_dat[9][it] = Mhpaf*Mhpaf*Mhpaf*Mhpaf;
      obs_dat[10][it] = Mafp*Mafp*Mafp*Mafp;




      }                       /* fin de 'itmax' */

      Escribe_Hister(ibin);  /* fichero de histeresis */

#ifndef HISTER
      Escribe_Result(ibin);
#endif  
  
      /* Control */

     printf("Energia Ios vecinos = %f\n",obs_dat[0][datos.itmax-1]);
     printf("Energia 2os vecinos = %f\n",obs_dat[1][datos.itmax-1]);
     printf("M ferromag. = %f\n",obs_dat[2][datos.itmax-1]);
     printf("M Hiperpla. = %f\n",obs_dat[3][datos.itmax-1]);
     printf("M Planos    = %f\n",obs_dat[4][datos.itmax-1]);
      
     datos.seed=RANDOM;
     datos.itcut=ibin+1;


#ifndef HISTER
     escribe_conf(0);
#endif


   }      /*  final de NBIN  */   

   fclose(Fheat);
   
   return 0;


}  /*      END         */


























