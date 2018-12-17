/* O4MED.C  MEDIDA DE ENERGIA Y PARAMETROS DE ORDEN */

#include "o4.h"
#include "o4aritm.h"

extern float ener1,ener2,E_1,E_2;

extern o4v f[],SUMA,suma1,suma2,cero,celda[];

extern float beta1,beta2,Mf,Mhpaf,Mafp;

extern int jpx,jpy,jpz,jpt,jmx,jmy,jmz,jmt,
           jpp[],jpm[],jmp[],jmm[];

extern int x_p[],x_m[],y_p[],y_m[],
           z_p[],z_m[],t_p[],t_m[];

extern int neigh_px,neigh_py,neigh_pz,neigh_pt,
           neigh_mx,neigh_my,neigh_mz,neigh_mt,it;

void Energia(int n)
{
  o4v jpp6,jpm6,jmp6,jmm6;
    
  suma2 = cero;
  suma1 = cero;
  SUMA = cero;

  Neigh(n);                         /* Calcula los vecinos */


                                   /*  ener de primeros vecinos */

  _suma8(suma1,f[jpx],f[jmx],f[jpy],f[jmy],f[jpz],
         f[jmz],f[jpt],f[jmt]); 

  _multesc(suma1,suma1,beta1);    


                                   /*  ener de segundos vecinos */

  _suma6(jpp6,f[jpp[0]],f[jpp[1]],f[jpp[2]],f[jpp[3]],f[jpp[4]],f[jpp[5]]); 
  _suma6(jpm6,f[jpm[0]],f[jpm[1]],f[jpm[2]],f[jpm[3]],f[jpm[4]],f[jpm[5]]);
  _suma6(jmp6,f[jmp[0]],f[jmp[1]],f[jmp[2]],f[jmp[3]],f[jmp[4]],f[jmp[5]]);
  _suma6(jmm6,f[jmm[0]],f[jmm[1]],f[jmm[2]],f[jmm[3]],f[jmm[4]],f[jmm[5]]);
 
  suma2.a0 = jpp6.a0 + jpm6.a0 + jmp6.a0 + jmm6.a0 ;
  suma2.a1 = jpp6.a1 + jpm6.a1 + jmp6.a1 + jmm6.a1 ;
  suma2.a2 = jpp6.a2 + jpm6.a2 + jmp6.a2 + jmm6.a2 ;
  suma2.a3 = jpp6.a3 + jpm6.a3 + jmp6.a3 + jmm6.a3 ;

  

  _multesc(suma2,suma2,beta2);

  _suma(SUMA,suma1,suma2);           /* necesario para el update  */

  ener1 = _prodescalar(f[n],suma1);
  ener2 = _prodescalar(f[n],suma2);
 
}

void Mide_Ener(void)
{
  
     int nn,x,y,z,t;

     nn=0;
     E_1=E_2=0.0;  

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

                 
                 Energia(nn);               
               
               
                 E_1 += ener1;
                 E_2 += ener2;

                 nn++;        
	       }
	    }
	 }
      }
  
     
#ifdef DEBUG
printf("E1 = %d\t%f\n",it,(E_1)*Normaener1);
printf("E2 = %d\t%f\n",it,(E_2)*Normaener2); 
#endif

}

void Order(void)
{

   int x,y,z,t,i,j,par[4],direccion,sp;
   int site;
   o4v MHPAF[4],MAFP[4][4],MF;
   float mhpaf[4],mafp[4][4],N;

   Mf=Mhpaf=Mafp=0.0;           /* Inicializacion necesaria */
   MF = cero;
   for(j=0;j<16;j++)
     {
       celda[j]=cero;
     }

   for(j=0;j<4;j++)   
       for(i=0;i<4;i++) 
       {
          MAFP[i][j] = cero;
          MHPAF[i] = cero;
       }

   site=0;

   for(t=0;t<L;t++)
      for(z=0;z<L;z++)     
         for(y=0;y<L;y++)
            for(x=0;x<L;x++)
	    {
                                       /* MF Ferro. */
	       MF.a0 += f[site].a0;
	       MF.a1 += f[site].a1;
               MF.a2 += f[site].a2;
               MF.a3 += f[site].a3;      

                                       /* MHPAF = AntiHiperplanos */

               par[0] = (x%2);         /* paridad comp. x */
               par[1] = (y%2);
               par[2] = (z%2);
               par[3] = (t%2);
               
               for(i=0;i<4;i++)
	       {
                  if ( par[i] == 0) par[i]=1;
                     else par[i] = -1;
                  MHPAF[i].a0 += f[site].a0*par[i];
                  MHPAF[i].a1 += f[site].a1*par[i];
                  MHPAF[i].a2 += f[site].a2*par[i];
                  MHPAF[i].a3 += f[site].a3*par[i];
	       }
          
                                         /* MAFP = Antiplanos  */
              for(j=3;j>=0;j--)
	      {
                 direccion=j;
                 for(i=0;i<direccion;i++)
		 {
                     sp = (par[i]*par[j]);
               
                     MAFP[j][i].a0 += f[site].a0*sp;
                     MAFP[j][i].a1 += f[site].a1*sp;
                     MAFP[j][i].a2 += f[site].a2*sp;
                     MAFP[j][i].a3 += f[site].a3*sp;
		   }

		}




                celda[(x%2)+2*(y%2)+4*(z%2)+8*(t%2)].a0 += f[site].a0;
                celda[(x%2)+2*(y%2)+4*(z%2)+8*(t%2)].a1 += f[site].a1;
                celda[(x%2)+2*(y%2)+4*(z%2)+8*(t%2)].a2 += f[site].a2;
                celda[(x%2)+2*(y%2)+4*(z%2)+8*(t%2)].a3 += f[site].a3;
               


                site++;

	     }        /* sweep en V */


                                      /* Medida de los parametros */
    N=Normaener*Normaener;

    Mf = (MF.a0*MF.a0 + MF.a1*MF.a1 
          + MF.a2*MF.a2 + MF.a3*MF.a3)*N;

    Mf = sqrt(Mf);
    for(i=0;i<4;i++)
    {
       mhpaf[i] = ( MHPAF[i].a0*MHPAF[i].a0 + 
                    MHPAF[i].a1*MHPAF[i].a1 +
                    MHPAF[i].a2*MHPAF[i].a2 + 
                    MHPAF[i].a3*MHPAF[i].a3)*N;

       Mhpaf += mhpaf[i];
    }

    Mhpaf = sqrt(Mhpaf);

     for(j=3;j>=0;j--)
     {
        direccion=j;
        for(i=0;i<direccion;i++)
	{
           mafp[j][i] = ( MAFP[j][i].a0*MAFP[j][i].a0 +
                           MAFP[j][i].a1*MAFP[j][i].a1 +
                            MAFP[j][i].a2*MAFP[j][i].a2 +
                           MAFP[j][i].a3*MAFP[j][i].a3)*N;
           Mafp += mafp[j][i];
	}
      }

      Mafp = sqrt(Mafp);
       
      for(i=0;i<16;i++)
	{
           _multesc(celda[i],celda[i],Normaener*16)
	 }


#ifdef DEBUG

printf("Parametro Ferromagnetico %f\n",Mf);
printf("Parametro de orden Hyperplanos %f\n",Mhpaf);
printf("Parametro de orden planos %f\n",Mafp);

#endif   



 }


