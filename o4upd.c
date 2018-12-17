
 /* O4UPD.C   FUNCIONES DE UPDATE  */

#include "o4.h"
#include "o4aritm.h"


extern o4v f[],SUMA,suma1,suma2;
extern int neigh_px,neigh_py,neigh_pz,neigh_pt,
           neigh_mx,neigh_my,neigh_mz,neigh_mt,it,good;
 
extern int jpp[],jpm[],jmp[],jmm[],jmp[],veces;

extern int jpx,jpy,jpz,jpt,
           jmx,jmy,jmz,jmt;

extern FILE *Fheat;
 
void Neigh(int j)
{
 
  jpx = j + neigh_px;        /* CALCULO DE LOS SITES VECINOS */
  jmx = j + neigh_mx;
  jpy = j + neigh_py;
  jmy = j + neigh_my;        /* 8 Primeros    vecinos   */
  jpz = j + neigh_pz;
  jmz = j + neigh_mz;
  jpt = j + neigh_pt;
  jmt = j + neigh_mt;


                             /* 24 segundos vecinos */
  jpp[0] = jpx + neigh_py;
  jpp[1] = jpx + neigh_pz;
  jpp[2] = jpx + neigh_pt;
  jpp[3] = jpy + neigh_pz;
  jpp[4] = jpy + neigh_pt;
  jpp[5] = jpz + neigh_pt;
  
  jpm[0] = jpx + neigh_my;
  jpm[1] = jpx + neigh_mz;
  jpm[2] = jpx + neigh_mt;
  jpm[3] = jpy + neigh_mz;
  jpm[4] = jpy + neigh_mt;
  jpm[5] = jpz + neigh_mt;

  jmp[0] = jmx + neigh_py;
  jmp[1] = jmx + neigh_pz;
  jmp[2] = jmx + neigh_pt;
  jmp[3] = jmy + neigh_pz;
  jmp[4] = jmy + neigh_pt;
  jmp[5] = jmz + neigh_pt;

  jmm[0] = jmx + neigh_my;
  jmm[1] = jmx + neigh_mz;
  jmm[2] = jmx + neigh_mt;
  jmm[3] = jmy + neigh_mz;
  jmm[4] = jmy + neigh_mt;
  jmm[5] = jmz + neigh_mt;

}



void Overrelaxed(int j)
{

     o4v fi_over,temp,proy,dif,dif_2;

     float factor,norm_temp,norm_SUMA,cosa; 

     Energia(j);

 
     temp=f[j]; 
    
    
     factor = 2.0F*_prodescalar(SUMA,temp)/_norm2(SUMA);

     fi_over.a0 = SUMA.a0 *factor - temp.a0;
     fi_over.a1 = SUMA.a1 *factor - temp.a1;
     fi_over.a2 = SUMA.a2 *factor - temp.a2;
     fi_over.a3 = SUMA.a3 *factor - temp.a3;  


 /* CAMBIO DEL VECTOR EN LA MATRIZ DE CAMPO */

       f[j] = fi_over;
     
     
}

  

void Canonic_Update(int j)
{

    float umag,umagin,damdum,umaginb,norma2,invnorma;
    float a_0,a_1,a_2,a_3,rad,red,x;
    o4v uint, anew,temp,temp2;
    float n,accion_old,accion_new,bolt;
    int nhit,delta2,i; 


    Energia(j);         /* calculo `stap'  */

    uint = SUMA;


#ifdef HEATBATH

      umag=sqrt(_norm2(uint)); 
   
      umagin=1.0F/umag;
    do
    {
        if ((x=FRANDOM))        
            a_0=1.0F+log(x)*umagin;
        else
            a_0=1.0F-100*umagin;        
 
        rad=1.0F-a_0*a_0;
       
        x=FRANDOM;
    }
    while (x*x > rad);
    
    /* fprintf(Fheat,"%g\n",rad);  PRINTS BOLTZMAN DISTRIBUTION */

    do
    {
        a_1=FRANDOM-0.5F;
        a_2=FRANDOM-0.5F;
        red=a_1*a_1+a_2*a_2;
    }
    while (red > 0.25F);

    x=FRANDOM;
    a_3=sqrt(rad)*(2.0F*x-1.0F);
    damdum=rad/red*4.0F*x*(1.0F-x); 
    umaginb=umagin*sqrt(damdum);

    anew.a0 = a_0*umagin;
    anew.a1 = a_3*umagin;
    anew.a2 = a_2*umaginb;
    anew.a3 = a_1*umaginb;


    _rotacion(temp,anew,uint);   /* Vuelve al espacio original */

    f[j] = temp;
    
#else

      /* SINO, HACE  METROPOLIS  */

    temp2=f[j];
    accion_old = _prodescalar(temp2,uint);

    nhit = 2;
    delta2 = 2;

    for(i=0;i<nhit;i++)
    {
       do
       {
         temp.a0 = FRANDOM - 0.5F;
         temp.a1 = FRANDOM - 0.5F;
         temp.a2 = FRANDOM - 0.5F;
         temp.a3 = FRANDOM - 0.5F;
       }
       while (_norm2(temp) > 0.25);

       _multesc(temp,temp,delta2);
                                            /* NORMALIZA */
       _suma(temp,temp,temp2);

       invnorma = 1.0F/sqrt(_norm2(temp));

       _multesc(temp,temp,invnorma);

       accion_new = _prodescalar(temp,uint);

       if(accion_new < accion_old)
       {
          x=FRANDOM;
          if (exp(accion_new-accion_old) > x)
	  {
             temp2 = temp;
             accion_old = accion_new;
             good++;
         
       
	  }
	}
        else
	{
           temp2=temp;
           accion_old=accion_new;
           good++;
      
	 }
     }
  
     f[j] = temp2;

#endif    
}
    




























 
