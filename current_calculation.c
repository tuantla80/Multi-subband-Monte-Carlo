/* *****************************************************************************
This functin is to calculate the current based on the charges entering and exiting each
   terminal contact -> noisy: due to discrete nature of the electrons -> fitting data

Starting date: April 6, 2010
Latest update: April 6, 2010
****************************************************************************** */
#include <stdio.h>
#include "constants.h"
#include "nrutil.h"

void current_calculation(double *cur_av){
     // Goi cac ham
     void linreg(double x[],double y[],double sigma[],int N,double a_fit[],double sig_a[],double yy[],double *chisqr );
     double Get_Vd_start(),Get_Vd_end(),Get_Vd_step();
     double Get_Vg_start(),Get_Vg_end(),Get_Vg_step();
     double Get_dt(),Get_tot_time(),Get_transient_time(); 
     // Bien local
     double Vd_start = Get_Vd_start();
     double Vd_end   = Get_Vd_end();
     double Vd_step  = Get_Vd_step(); 
     double Vg_start = Get_Vg_start();
     double Vg_end   = Get_Vg_end();
     double Vg_step  = Get_Vg_step();
     double dt = Get_dt();                         printf("\n dt=%le",dt);
     double tot_time = Get_tot_time();             printf("\n total_time =%le",tot_time);
     double transient_time = Get_transient_time(); printf("\n transient_time=%le",transient_time); 
    
// Thuc hien          
     double *rtime,*rtime2,*source,*source2,*drain,*drain2,*sigma,*yy,*a_fit,*sig_a;
     int i, N=0, N_process=0;//N: is the dimension of the matrix -> to avoid RUNTIME ERROR
     N_process=((int)((Vg_end-Vg_start)/Vg_step) + 1)*((int)((Vd_end-Vd_start)/Vd_step) + 1); 
     int t1 = (int)(tot_time/1.0e-16);// vi khi doc vao cac gia tri thi neu DO CHINH XAC o input file khac nhau thi gia tri thuc te se khac nhau 1 ti
     int t2 = (int)(transient_time/1.0e-16);
     int t3 = (int)(dt/1.0e-16); //printf("\n t1=%d, t2=%d, t3=%d", t1,t2,t3);
     N = (int)((t1-t2)/t3)*N_process; printf("\n N=%d,N_process=%d",N,N_process);//getchar();

     rtime   = dvector(1,N);
     rtime2  = dvector(1,N);
     source  = dvector(1,N);
     source2 = dvector(1,N);
     drain   = dvector(1,N);
     drain2  = dvector(1,N);
     sigma   = dvector(1,N);
     yy      = dvector(1,N);
     a_fit   = dvector(1,2);
     sig_a   = dvector(1,2);
     for(i=1; i<=N; i++){
         rtime[i] = 0.0; rtime2[i] = 0.0; source[i] = 0.0; source2[i] = 0.0; 
         drain[i] = 0.0; drain2[i] = 0.0; sigma[i] = 0.0; yy[i] = 0.0;
     }
     for(i=1; i<=2; i++){ a_fit[i] = 0.0; sig_a[i] = 0.0; }
       
     //double width = device_width*1.0e6; // [micro m] 
    
     FILE *ff_transient;
     ff_transient=fopen("charge_source_drain_transient.dat","r"); 
     if(ff_transient==NULL){printf("can't open charge_source_drain_transient.dat\n");return;} 
     int j=1;
     while (!feof(ff_transient)){
           fscanf(ff_transient,"%le %le %le ",&rtime[j],&source[j],&drain[j]); 
           //NOTE: rtime[]: la thoi gian chay -> don vi se la [s]
           //      source[]:la tong so hat vao ra o cuc Source -> don vi se la [C] culong
           //      drain[]: la tong so hat vao ra o cuc Drain -> don vi se la [C] culong  
           j = j + 1;// Ket qua cuoi cung j la tong so dong trong file charge_source_drain_transient.dat
      }
     fclose(ff_transient);  
     int j_fix = j - 1; // j_fix la tong so dong trong file charge_source_drain_transient.dat
     int jskip = (int)(j_fix/2); // la 1/2 so dong trong file charge_source_drain_transient.dat
     
     FILE *f;
     f=fopen("calculated_current_2_txt.dat","a"); // Chua thay co ung dung gi ca
     if(f==NULL){printf("Can't open calculated_current_2_txt.dat.\n");return;}
     fprintf(f,"\n # rtime[s]  source   drain  \n");
                   
     double sum_source = 0.0;
     double sum_drain = 0.0;
      
     int jj=0;
     j=0;
     for(j=1;j<=j_fix;j++){ // di tu dong dau tien den dong cuoi cung trong file charge_source_drain_transient.dat
        sum_source += source[j]; // tong so hat (o cuc Source) theo tat ca cac dong
        sum_drain  += drain[j];   // tong so hat (o cuc Drain) theo tat ca cac dong
        if(j >jskip) { // Ta chi dung cac lenh sau cho nhung dong o nua sau cua file charge_source_drain_transient.dat
           jj = jj + 1; // tinh so dong o nua sau, chu y ban dau jj=0 nhe
           rtime2[jj] = rtime[j];
           sigma[jj] = 1.0;
           source2[jj] = -sum_source;
           drain2[jj] = sum_drain;
           fprintf(f,"%le %le %le \n ",rtime[j],-sum_source,sum_drain); 
           // Xem hinh 15 quyen sach cua GS do
        }// End of if
     }// End of for
     fclose(f);
        
     double chisqr = 0.0;
     linreg(rtime2,source2,sigma,jj,a_fit,sig_a,yy,&chisqr); // min cho Source
        /* Nhu vay co the thay ta chi lam min o nua sau cua file 
           charge_source_drain_transient.dat thoi */
     //double cur_source = 1.602e-19*a_fit[2]/width*1.e6;
     double cur_source = q*a_fit[2];// Tu ta them vao

     linreg(rtime2,drain2,sigma,jj,a_fit,sig_a,yy,&chisqr); // min cho Drain
     //double cur_drain = 1.602e-19*a_fit[2]/width*1.e6;
     double cur_drain = q*a_fit[2];
       
     //printf("\n Source current = %le",cur_source);
     //printf("\n Drain_current = %le",cur_drain);
     *cur_av = 0.5*(cur_source + cur_drain);
     printf("\n Tai current_calculation():  Average current = %le",0.5*(cur_source + cur_drain)); //getchar();

      // Free for local vectors
      //*  Khong hieu tai sao neu dung free_dvector o dat thi lai co loi segmentation fault (29-Jan-2009)
       free_dvector(rtime,1,N);
       free_dvector(rtime2,1,N);
       free_dvector(source,1,N);
       free_dvector(source2,1,N);
       free_dvector(drain,1,N);
       free_dvector(drain2,1,N);
       free_dvector(sigma,1,N);
       free_dvector(yy,1,N);
       free_dvector(a_fit,1,2);
       free_dvector(sig_a,1,2);
       // */
      
      return;
  }// End of current_calculation()
/*********************************************************************************
  Function to perform linear regression (fit a line)
  Inputs
   x       Independent variable
   y       Dependent variable
   sigma   Estimated error in y
   N       Number of data points
  Outputs
   a_fit   Fit parameters; a(1) is intercept (phan chan), a(2) is slope
   sig_a   Estimated error in the parameters a()
   yy      Curve fit to the data
   chisqr  Chi squared statistic
   
   Reference: fit.c - Numerical Recipes for C - page 665
              "Fitting data to a Straight Line"
 
 Latest update January 7, 2009 
***********************************************************************************/
 // Vi x,y, sigma la cac mang dau vao va khong thay doi nen ta co the them constant dang truoc
 // vi du: void linreg(constant double x[],constant double y[],constant double sigma[]
  
 #include<math.h>
 void linreg( double x[],double y[],double sigma[],int N ,double a_fit[],double sig_a[],double yy[],double *chisqr )
{
  int i = 0;
  double sigmaTerm = 0.0, s = 0.0, sx = 0.0, sy = 0.0;
  double sxy = 0.0, sxx = 0.0, denom = 0.0, delta = 0.0;

  // Evaluate various sigma sums
      s = 0.0;
      sx = 0.0;
      sy = 0.0;
      sxy = 0.0;
      sxx = 0.0;
      
      for(i=1;i<=N;i++){
        sigmaTerm = 1.0/sigma[i]*sigma[i];
        s = s + sigmaTerm;
        sx = sx + x[i]*sigmaTerm;
        sy = sy + y[i]*sigmaTerm;
        sxy = sxy + x[i]*y[i]*sigmaTerm;
        sxx = sxx + x[i]*x[i]*sigmaTerm;
      }
      denom = s*sxx - sx*sx;

      // Compute intercept a_fit[1] and slope a_fit[2]
      a_fit[1] = (sxx*sy - sx*sxy)/denom;
      a_fit[2] = (s*sxy - sx*sy)/denom;

      // Compute error bars for intercept and slope
      sig_a[1] = sqrt(sxx/denom);
      sig_a[2] = sqrt(s/denom);

      // Evaluate curve fit at each data point and compute Chi^2
      *chisqr = 0.0;
      for(i=1;i<=N;i++){
        yy[i] = a_fit[1]+a_fit[2]*x[i]; // Curve fit to the data
        delta = (y[i]-yy[i])/sigma[i];
        *chisqr = *chisqr + delta*delta; // Chi square
      }
      return;
}// End of linreg() function
    
      

