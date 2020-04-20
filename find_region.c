/* ****************************************************************
Tu 1 diem (i,j,k) la toa do cua x,y va z tuong ung ->
     ta biet duoc diem nay nam o region nao (S, D or Channel)
     To find a doping region based on node location
     Region = 1 -> Source
     Region = 2 -> Drain
     Region = 3 -> Channel
     
Starting date: Feb 11, 2010
Latest update: Feb 11, 2010
****************************************************************** */
#include <stdio.h>

int find_region(int i,int j,int k){
    // Goi cac ham
    int Get_nx0(),Get_nxa(),Get_nxb(),Get_nx1();
    int Get_ny0(),Get_nya(),Get_nyb(),Get_ny1();
    int Get_nz0(),Get_nza(),Get_nzb(),Get_nz1();
    
    // Cac bien local
    int nx0 = Get_nx0();
    int nxa = Get_nxa();
    int nxb = Get_nxb();
    int nx1 = Get_nx1(); 
    int ny0 = Get_ny0();
    int nya = Get_nya();
    int nyb = Get_nyb();
    int ny1 = Get_ny1();
    int nz0 = Get_nz0();
    int nza = Get_nza();
    int nzb = Get_nzb();
    int nz1 = Get_nz1();
    
    // Thuc hien
    int i_reg = 0; 
    if((i>=nx0)&&(i<=nxa)&&(j>=nya)&&(j<=nyb)&&(k>=nza)&&(k<=nzb))
         { i_reg = 1; } // Source 
    else if((i>nxa)&&(i<nxb)&&(j>=nya)&&(j<=nyb)&&(k>=nza)&&(k<=nzb))
         { i_reg = 3; } // Channel. RAT CHU Y
    else if((i>=nxb)&&(i<=nx1)&&(j>=nya)&&(j<=nyb)&&(k>=nza)&&(k<=nzb))
         { i_reg = 2; } // Drain
    else { i_reg = 0; } // Outside the regions we consider                                                             
                                                                  
    return i_reg ; //Tuc la ham nay se dua ra duoc o vi tri [i,j] thi la region so may
}


















