/*  PROGRAMA DE ANALISIS DEL MODELO O(4) Antiferro */


# include "ana.h"



void muestra_datos(struct s_datos *dat)
{

	printf("\n El input de los archivos es: \n");
	printf("itmax  %d \n", dat->itmax);
	printf("mfresh  %d \n", dat->mfresh);
	printf("nbin  %d \n", dat->nbin);
	printf("itcut  %d \n", dat->itcut);
	printf("flag  %d \n", dat->flag);
	printf("seed  %d \n", dat->seed);
	printf("\n");
	printf("beta1  %f \n", dat->beta1);
	printf("beta2  %f \n", dat->beta2);
	printf("\n");
}

void lee_argumentos( int argc, char *argv[],
		     struct s_datos *datos,
		     float *beta1, float *beta2)
{
	int n,nbl,lbl;
	char nombres[100], name[100];
	FILE *Finput, *Fnombres;
	static struct s_datos datos_old;
	

	switch(argc)
	{
	 case 4: sscanf(argv[3],"%d", &nbl);
	 case 3: sscanf(argv[2],"%d", &n2);
		 sscanf(argv[1],"%d", &n1); break;
	 default: printf("ERROR: para lanzar introduce\n ");
		  printf(" ana%d (1er fich) (ult. fich)",L);
		  printf(" (n_bloques)\n ");
		  exit(0);
	 }
        nblo=nbl;

	 if(n2-n1+1>maxfiles)
	 {
	    printf("MAXIMUM NUMBER OF FILES = %i \n", maxfiles);
	    exit(0);
	 }

	 printf("\n STARTING \n");
	 printf(" ANALYSIS OF OUT%03d.DAT - OUT%03d.DAT\n",n1,n2);
	 printf(" LATTICE SIZE %dx%dx%dx%d \n",L,L,L,L);
	 printf(" JACK-KNIFE BLOCKS %d\n ",nblo);

	 if((nblo)>maxbloque)
	 {
		printf("\n nblo  no puede superar %d \n",maxbloque);
		exit(0);
	 }

	 lbl=(n2-n1+1)/nbl;
	 if(lbl==0)
		exit(0);
	 lblo=lbl;
	 n1=n2-lblo*nblo+1;

	 for(n=n1;n<=n2;n++)
	 {
		sprintf(name,"%s%03d.DAT",nom_fich,n);
		if ((Finput=fopen(name,"rb"))==NULL)
		{
			printf("\n THIS FILE DOES NOT EXIST \n",name);
			exit(0);
		}

		fread(datos,sizeof(*datos),1,Finput);
	
		if(n==n1)
		{
		   memcpy(&datos_old,datos, sizeof(datos_old));
		}
		if(n>n1 && (datos->beta1 != datos_old.beta1 ||
			    datos->beta2 != datos_old.beta2 ||
			    datos->itmax != datos_old.itmax ||
			    datos->mfresh != datos_old.mfresh))
		{
			printf("NOT THE SAME INPUT FROM %d \n",n);
			exit (0);

		}
		fclose(Finput);

		if((ficheros[n-n1]= (char*)malloc(strlen(name)))==NULL)
		{
			printf("Fuera de memoria \n");
			exit(0);
		}
		strcpy(ficheros[n-n1],name);
	
          }

	 printf("\n %d \n",datos->itmax);
	 muestra_datos(datos);
	 *beta1=datos->beta1;
	 *beta2=datos->beta2;

}

void lee_datos(int n)
{
	int idat,idatnew,imag;
	FILE *Finput;
	struct s_datos dat;

	
	Finput=fopen(ficheros[n],"rb");
	fread(&dat,sizeof(dat),1,Finput);
         
  	   for(idat=0;idat<n_obs_medid;idat++)	     	     
	     fread(&v_dat[idat][0],4,(size_t)dat.itmax,Finput);
       
           for(idatnew=0;idat<dat.itmax;idat++)
	   {
              for(imag=0;imag < 3;imag++)
              {

/* M*E_1 */      v_dat[11+imag][idatnew] =
                       (v_dat[0][idatnew]*v_dat[2+imag][idatnew]);

/* M*E_2 */      v_dat[14+imag][idatnew] =
                       (v_dat[1][idatnew]*v_dat[2+imag][idatnew]);

	       }

/* E_1 ^2 */      v_dat[17][idatnew] = v_dat[0][idatnew]*v_dat[0][idatnew];
/* E_2 ^2 */      v_dat[18][idatnew] = v_dat[1][idatnew]*v_dat[1][idatnew];

              
	   }

           
	
  	fclose(Finput);
}


void Histograma(int acoplo)
{
          int ib,j;
	  float frecu[n_inter+1],y;
	  FILE *Fhistog;
	  char nomhist[150];

	  sprintf(nomhist,"histog%d_%5.4f_%5.4f_%d.plt",acoplo,datos.beta1,datos.beta2,L);
	  Fhistog=fopen(nomhist,"w");
	  fprintf(Fhistog,"title top fon T 'RETICULO L=%d, beta1=%f, beta2=%f, HISTOGRAMA de le energia a %dos vecinos'\n",L,datos.beta1,datos.beta2,(acoplo+1));
	  fprintf(Fhistog,"title left font T 'frecuencia'\n");
	  fprintf(Fhistog,"title bottom font T 'Espectro de energia a %dos vecinos'\n",(acoplo+1));
      
	  for(j=0;j<=n_inter;j++)
	    frecu[j]=0.;
	  for(j=0;j<=n_inter;j++)
	   for(ib=0;ib<nblo;ib++)
	    frecu[j]+=frec[ib][j];
	  for(j=0;j<=n_inter;j++)
	  {
	    y=c1*j+c2;
	    fprintf(Fhistog,"%f\t%f\n", y,frecu[j]);
	  }
	  fprintf(Fhistog,"bargraph solid fullwidth\n");
	  fclose(Fhistog);
}

void FS(void)
{
	int ib,ibb,i,j,k,q,p,volu;
	double SumO[n_obs_FS],SumO2[n_obs_FS],SumderO[n_obs_FS],
               SumderO2[n_obs_FS],SumderlO[n_obs_FS],SumderlO2[n_obs_FS];
	double sumO[n_obs_FS],sumOe[n_obs_FS],sumx[n_obs_FS];
	double expo,expo_frec,O[n_obs_FS],derO[n_obs_FS],derlO[n_obs_FS];
	long int sumf;
	double en, sum,sume,x,y,lO[n_obs_FS];
       

    
        for(k=0;k<n_obs_FS;k++)
         for(ib=0;ib<=nblo;maxder[k][ib++]=0);

 	volu=L*L*L*L;

	for(j=0;j<=nbetas;j++)
	{
	 for(k=0;k<n_obs_FS;k++)
	 {
		SumO[k]=0.;
		SumO2[k]=0.;
		SumderO[k]=0.;
		SumderO2[k]=0.;
               
	 }

	 for(ib=0;ib<=nblo;ib++)
	 {
		sum=sume=0.;

		for(k=0;k<n_obs_FS;k++)
		 {
		  sumO[k]=0.;
		  sumOe[k]=0.;

		 }

		 x= x0-delta + j*h;

		 for(i=0;i<n_inter;i++)
		 {
			sumf=0.;
			for(k=0;k<n_obs_FS;k++)
				sumx[k]=0.;
			y=c1*i+c2;
			for(ibb=0;ibb<nblo;ibb++)
				if(ib!=ibb)
				{
					sumf += frec[ibb][i];

					for(k=0;k<n_obs_FS;k++)
						sumx[k] += xfrec[k][ibb][i];
				}

			expo=exp((x-x0)*(y-ymed)*vol);
			expo_frec = expo*sumf;
			sum += expo_frec;
			sume += y*expo_frec;

			for(k=0;k<n_obs_FS;k++)
			{
				sumO[k] += sumx[k]*expo;
				sumOe[k] += sumx[k]*y*expo;
			}
		 }

		 en=sume/sum;
		 for(k=0;k<n_obs_FS;k++)
		 {
		     O[k]=sumO[k]/sum;
	             derO[k]=(sumOe[k]/sum - O[k]*en)*vol;
           	    
		 }

		
		 if(ib<nblo)
		  for(k=0;k<n_obs_FS;k++)
		  {
			SumO[k] += O[k];
			SumO2[k] += O[k]*O[k];
                        SumderO[k] += derO[k];
                        SumderO2[k] += derO[k]*derO[k];
                        
		   }
			
                 for(k=0;k<n_obs_FS;k++)
		    if(fabs(derO[k]) > fabs(maxder[k][ib]))
		    { 
                        
                       maxder[k][ib]=derO[k];
                       coup_maxder[k][ib]=x;
		    }

		 
	     } /* nblo */
		for(k=0;k<n_obs_FS;k++)
		{
			SumO[k]/=nblo;
		        SumderO[k]/=nblo;
			SumderlO[k]/=nblo;
                        O_v[k][j]=O[k];
		        derO_v[k][j]=derO[k];
			O_err[k][j]=
                                   sqrt(fabs(SumO2[k]/nblo-SumO[k]*SumO[k])*
                                   (double)(nblo-1));

                        derO_err[k][j]=
                                   sqrt(fabs(SumderO2[k]/nblo -
                                   SumderO[k]*SumderO[k])*(double)(nblo-1));

		}

             
	} /* bucle en betas */

  
}
void Maxder(void)
{
	int ib,ibb,i,j,k,q,p,volu;
	double SumO[n_obs_FS],SumO2[n_obs_FS],SumderO[n_obs_FS],SumderO2[n_obs_FS];
	double sumO[n_obs_FS],sumOe[n_obs_FS],sumx[n_obs_FS];
	double expo,expo_frec,O[n_obs_FS],derO[n_obs_FS];
        long double Sum, Sum2,Sumder,Sumder2;
	long int sumf;
	double en, sum,sume,x,y;


	volu=L*L*L*L;
	 for(ib=0;ib<=nblo;ib++)
	 {

	  for(k=0;k<n_obs_FS;k++)
	  {

		  x=coup_maxder[k][ib]-h/2-h/20;
		  for(j=-10;j<=10;j++)
		  {
			x+=h/20; 

			if(k<n_obs_FS)
			{

			  sumO[k]=0.;
			  sumOe[k]=0.;
			  sum=sume=0.;


			 for(i=0;i<n_inter;i++)
			 {
				sumf=0.;
			      
				sumx[k]=0.;

				for(ibb=0;ibb<nblo;ibb++)
					if(ib!=ibb)
					{
					sumf += frec[ibb][i];
					sumx[k] += xfrec[k][ibb][i];
					}

				expo=exp((x-x0)*(y-ymed)*vol);
				expo_frec = expo*sumf;
				sum += expo_frec;
				sume += y*expo_frec;

	
	
					sumO[k] += sumx[k]*expo;
					sumOe[k] += sumx[k]*y*expo;
				
			 }

			 en=sume/sum;
		      
			 
			 O[k]=sumO[k]/sum;

/* esto es la derivada */

		         derO[k]=(sumOe[k]/sum-O[k]*en)*vol;
			 
			}
			
		   
			if((derO[k])>(maxder[k][ib]))
			{
				maxder[k][ib]=derO[k];
				coup_maxder[k][ib]=x;
			}
	  }/* bucle en betas */
  

	 } /* bucle en observable  */
    
	} /* bucle en bloques */

	for(k=0;k<n_obs_FS;k++)
	{
		Sum=Sum2=Sumder=Sumder2=0.;
		for(ib=0;ib<nblo;ib++)
		{
			Sum += coup_maxder[k][ib];
			Sum2 += coup_maxder[k][ib]*coup_maxder[k][ib];
			Sumder += maxder[k][ib];
			Sumder2 += maxder[k][ib]*maxder[k][ib];
		}
		Sum/=nblo;
		err_coup[k]=sqrt(fabs(Sum2/nblo-Sum*Sum)*(double)(nblo-1));
		Sumder/=nblo;
		err_maxder[k]=sqrt(fabs(Sumder2/nblo-Sumder*Sumder)*(double)(nblo-1));
	}
}


Binder_Par(int nE)     /* Calcula el parametro de Binder de las Energias */
{

  double binder,x,cumul;
  int j,imag,vol;
  FILE *FBINDER,*FXC,*FKAPPA;
  char name1[100],name2[100],name3[100];
  double sm2,sm4,avm2_2,avm2_4,avm2_6,avm4_2;
  double sm,binder_error,jen,jen_err;


  vol=L*L*L*L;

  for(imag=0;imag<3;imag++)        /* para los tres parametros de orden */
  {

      sprintf(name1,"E%dbinder%d_%d.plt",nE,imag,L);
      
      sprintf(name3,"E%dkappa%d_%d.plt",nE,imag,L);
      
      sprintf(name2,"E%dxc%d_%d.plt",nE,imag,L);
     
      FXC=fopen(name2,"w");
      FBINDER=fopen(name1,"w");
      FKAPPA=fopen(name3,"w");

      for(j=0;j<=nbetas;j++)
      {
                jen = O_v[nE][j];

                cumul = O_v[8+imag][j]/(O_v[5+imag][j]*
                                        O_v[5+imag][j]);
                
                B[j] = (1.0 - (cumul/3.0));

/* KAPPA es un desastre */

                kappa[j] = vol*((O_v[11+imag+3*nE][j]/O_v[2+imag][j]) -
                                 O_v[nE][j]);
                
                /* calculo del error en el cumulante */

/************************************************************/

               sm2 = O_err[5+imag][j]*O_err[5+imag][j];
               sm4 = O_err[8+imag][j]*O_err[8+imag][j];
               avm2_2 = O_v[5+imag][j]*O_v[5+imag][j];
               avm2_4 = avm2_2 * avm2_2;
               avm2_6 = avm2_4 * avm2_2;
               avm4_2 = O_v[8+imag][j]*O_v[8+imag][j];

               binder_error = (sm4/avm2_4) + 0.2222*(sm2/avm2_6)*avm4_2;

/************************************************************/  

                B_err[j] = sqrt(binder_error);	
                
                x = x0 - delta + j*h;

/***************** SUSCEPTIBILIDAD conexa ************************/

               xc[j] = vol*x*(O_v[5+imag][j] -
                       (O_v[2+imag][j]*O_v[2+imag][j])); 
               
               sm = O_err[2+imag][j]*O_err[2+imag][j];

               xc_err[j]= sqrt(vol*vol*x*x*
                                ((O_err[5+imag][j]*O_err[5+imag][j]) +
                                 2*O_v[2+imag][j]*O_v[2+imag][j]*sm));

/*******************************************************************/

                fprintf(FBINDER,"%f\t%f\t%f\n",x,B[j],B_err[j]);  
                fprintf(FXC,"%f\t%f\t%f\n",x,xc[j],xc_err[j]);  
                fprintf(FKAPPA,"%f\t%f\n",x,kappa[j]);
                 
                       
      }
      fprintf(FBINDER, "join\n");
      fclose(FBINDER);
      fprintf(FXC,"join\n");
      fclose(FXC);
      fprintf(FKAPPA,"join\n");
      fclose(FKAPPA);
      
   }

  
}


void titulo (void)
{
	int i;
	sprintf(nombre,"title top font T 'L=%d, bb1e=%f, bb2e=%f, 
                   #MC= %dx%d Nbl=%d, Lblo=%d, %s'",
		   L,datos.beta1,datos.beta2,datos.itmax,datos.mfresh,
                   nblo, lblo, nom);
}

void Dibuja_resultados (int v)
{
	int j;
	double x;
	char formato[]= "new plot \n"
	      		"%s \n"
		        "title left font T '%s' \n"
			"title right font T '%s' \n"
			"case               '%s' \n"
			"set order x y dy \n"
			"set labels left font T \n"
		        "set labels bottom font T \n"
			"set symbol %s \n";


	/* OBSERVABLE  */

	fprintf(Foutplt,formato,nombre,"O",n_coup,"G- -","9N");

	for(j=0;j<=nbetas;j++)
	{
		x=x0-delta + j*h;
		fprintf(Foutplt,"%10f  %10f\n",x,O_v[v][j]);
	}
	fprintf(Foutplt,"join \n");

	for(j=0;j<=nbetas;j+=(nbetas/10))
	{

		x=x0-delta + j*h;
		fprintf(Foutplt,"%10f  %10f  %10f\n",x,O_v[v][j], O_err[v][j]);
	}
	fprintf(Foutplt,"plot\n");
	fprintf(Foutplt,"set symbol 5N \n");
	fprintf(Foutplt, " %10f %10f %10f\n", x0,O_v[v][nbetas/2],O_err[v][nbetas/2]);
	fprintf(Foutplt, "plot \n");

	/* DERIVADA  */
	fprintf(Foutplt, formato, nombre,n_deri,n_coup,"G- -","9N");
    
	for(j=0;j<=nbetas;j++)
	{
		x=x0-delta + j*h;
		fprintf(Foutplt,"%10f  %10f\n",x,derO_v[v][j]);
	}
	fprintf(Foutplt,"join \n");

	for(j=0;j<=nbetas;j+=(nbetas/10))
	{

		x=x0-delta + j*h;
       	fprintf(Foutplt,"%10f  %10f  %10f\n",x,derO_v[v][j], derO_err[v][j]);
	}
	
       
        fprintf(Foutplt,"plot\n");
	fprintf(Foutplt,"set symbol 3N \n");
	fprintf(Foutplt,"set order x dx y dy \n");
	fprintf(Foutplt,"%10f %10f %10f %10f\n",coup_maxder[v][nblo],err_coup[v],maxder[v][nblo],err_maxder[v]);
	fprintf(Foutplt,"set order x y dy \n");
	fprintf(Foutplt, " %10f %10f %10f \n", x0,derO_v[v][nbetas/2],derO_err[v][nbetas/2]);
	fprintf(Foutplt, "plot \n");

	fprintf(Foutplt,"title bottom font T' C=%8.6f+/-%8.6f; D=%8f+/-%8f'\n",
					coup_maxder[v][nblo],err_coup[v],maxder[v][nblo],err_maxder[v]);


}






void main(int argc,char *argv[])
{

	char nombreresult[150],n_der[2][10],n_cou[2][5],nombrefichero[150];
	int  n_vo[2],j,n,it,ib,t,i,p,Num_Ener;
	float datos_beta[2],sum,sumen,sumen2,ymin,ymax,beta1,beta2,coup;
        float en,cv,sigma,y;
	FILE *Foutput;




     
	strcpy(n_cou[0],"b1");
      
	strcpy(n_cou[1],"b2");

	lee_argumentos(argc,argv,&datos,&beta1,&beta2);

	n_vo[0]=1;
	n_vo[1]=1;
	datos_beta[0]=datos.beta1;
	datos_beta[1]=datos.beta2;

	printf("Los observables analizados son:\n");
	for(j=0;j<n_obs_FS;j++)
		printf("          %d\t%s\n",j,cadena[j]);
       
  

	for(d=1;d>=0;d--)
	{

               sprintf(nombreresult,"%dresult_%5.4f_%5.4f_%d.plt",
                       d,beta1,beta2,L);
        	if((Foutplt=fopen(nombreresult,"w"))==NULL)
	        {
         		puts("No puedo abrir el archivo de resultados\n");
		        exit(1);
	        }

	      
		strcpy(n_coup,n_cou[d]);
		n_vol=n_vo[d];
		coup=datos_beta[d];
		vol=(long int)n_vol*L*L*L*L;
		

		sum=sumen=sumen2=0.;
		ymin=10.;               /* Valores maximo y minimo de E */
		ymax=-10.0;

		for(n=n1;n<=n2;n++)
		{
			lee_datos(n-n1);
			for(it=0;it<datos.itmax;it++)
			{
				y=v_dat[d][it];
 				sum++;
				sumen += y;
				sumen2 += y*y;
				ymin=(y<ymin)?y:ymin;
				ymax=(y>ymax)?y:ymax;				
			}
		}
		en=sumen/sum;                  /* Dispersion de la gaussiana */
		cv=fabs(sumen2/sum-en*en);     /* para hacer FS              */
		sigma=sqrt(cv);

		for(n=0;n<n_obs_FS;n++)
			for(ib=0;ib<nblo;ib++)
				for(j=0;j<=n_inter;j++)
					xfrec[n][ib][j]=frec[ib][j]=0.;

		h=1.0/(ymax-ymin);
		c1=1.0/(float)n_inter/h;
		c2=-0.5/(float)n_inter/h+ymin;

		for(n=n1;n<=n2;n++)

		{                            /* Clasificacion en intervalos */
			ib=(n-n1)/lblo;      /* de los Observ. */
			lee_datos(n-n1);
			for(it=0;it<datos.itmax;it++)
			{
			        y=v_dat[d][it];
				i=(int)((y-ymin)*h*n_inter);

				frec[ib][i]++;
				for(t=0;t<n_obs_FS;t++)
					xfrec[t][ib][i]+=v_dat[t][it];
		
			}     
		}
		delta=3/(sigma*vol);  /* validez de la aproximacion de */
		h=2.0*delta/nbetas;     /* FS = n / (sigma*V)            */
		ymed=en;
		x0=coup;
		Histograma(d);
                FS();	      
                Maxder();                            
	

        for(n=0;n<n_obs_plot;n++)
	{

           sprintf(nom," O= %s ",cadena[n]);
           titulo();
           Dibuja_resultados(n);
	}
        
        Binder_Par(d);
 

    fclose(Foutplt);
    }

    

    for(n=0;n<n_obs_FS;n++)
    {
       printf("\n %d. %s= %8.6f +/- %8.6f \n",n,cadena[n],
               O_v[n][nbetas/2],O_err[n][nbetas/2]);
       sprintf(nombrefichero,"%dvsb_%d.plt",n,L);

       if((Foutput=fopen(nombrefichero,"a"))==NULL)
       {
          puts("\n No puedo abrir el fichero vsb.plt ");
          exit(1);
       }

       fprintf(Foutput,"%f %f\t%f\t%f\n",beta1,beta2,O_v[n][nbetas/2],
               O_err[n][nbetas/2]);
       fclose(Foutput);
     }

   



}          /*  END */
	
	














