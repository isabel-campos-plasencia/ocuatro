/* PROGRAMA DE ANALISIS PARA EL MODELO SIGMA ANTIFERROMAGNETICO*/
		     /*    HEADER	*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


# define L       6	 
# define maxit 		5000
# define maxbloque 	 20
# define nom_fich	 "OUT"
# define n_obs_medid      11  /* estos se miden directamente del output */    
# define n_obs_FS	 19   /* los otros 8 se calculan */
# define maxfiles	 5000
# define nbetas  	 10
# define n_inter	 100
# define n_correl	 100
# define maxmed		 150000
# define n_obs_plot      5


float v_dat[n_obs_FS][maxit],B[nbetas+1],B_err[nbetas+1];
long int frec[maxbloque][n_inter+1];
long int vol;
float xc[nbetas+1],xc_err[nbetas+1],kappa[nbetas+1],kappa_err[nbetas+1];
double xfrec[n_obs_FS][maxbloque][n_inter+1],SumO[n_obs_FS];
double coup_maxder[n_obs_FS][maxbloque+1], err_coup[n_obs_FS];
double maxder[n_obs_FS][maxbloque+1], err_maxder[n_obs_FS];
double O_v[n_obs_FS][nbetas+1], O_err[n_obs_FS][nbetas+1];
double derO_v[n_obs_FS][nbetas+1], derO_err[n_obs_FS][nbetas+1];
double ymed,x0,h,c1,c2,delta,derlO_v[n_obs_FS][nbetas+1],
        derlO_err[n_obs_FS][nbetas+1];
char   *ficheros[maxfiles],n_coup[10],n_deri[10], nombre[150],nom[100];
int    n_vol,permutar,d;

FILE *Foutplt;

int nblo,lblo,n1,n2;


struct s_datos
{
 int	itmax,
		mfresh,
		  nbin,
		itcut,
		flag,
		seed;
 float	beta1,
		beta2;
 };

 struct s_datos datos;
 char cadena[n_obs_FS][100]={
		      		 "Energia primeros vecinos",
	      			 "Energia segundos vecinos",
				 "Modulo de magnetizacion",
				 "Magnetizacion a hiperpl.",
		      		 "Magnetizacion a planos",
                                 "Cuadrado de M ",
                                 "Cuadrado de Mhpaf",
                                 "Cuadrado de Mafp ",
                                 " M to the 4th    ",
                                 " Mhpaf to the 4th ",
/* **** medidos *********** */   " Mafp  to the 4th ",
                                 " M0 * E1          ",
                                 " M1 * E1          ",
                                 " M2* E1        ",
                                 " M0 * E2        ",
                                 " M1 * E2       ",
                                 " M2 * E2       ",
                                 " E1^2          ",
                                 " E2^2          "
	      			};


