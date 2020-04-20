/* ***********************************************************************
To check the charge neutrality of the source and drain ohmic contacts
   + Neu so hat LON HON thi Eliminate no
   + Neu so hat NHO HON thi Create  
   + Considering for SIDE CONTACT case
Co nghia chi xac dinh cai phan cua Source va Drain ap vao

Starting date: March 27, 2010
Latest update: March 29, 2010 
************************************************************************* */
#include <stdio.h>

void check_source_drain_contacts(){
     // Goi cac ham
     void init_kspace(int ne,int i),init_realspace(int ne,int i); 
     int Get_nsource_side_carriers(),Get_ndrain_side_carriers();
     int Get_n_used(); // Tinh SO HAT dang SU DUNG
     void Set_n_used( int n);
     double **Get_p();
     int *Get_valley();
     int Get_iss_eli(),Get_idd_eli(),Get_iss_cre(),Get_idd_cre();
     void Set_iss_eli(int eli_p),Set_idd_eli(int eli_p),Set_iss_cre(int cre_p),Set_idd_cre(int cre_p);
     double Get_mesh_size_x();
              
     // Cac bien local
     int nsource_side_carriers = Get_nsource_side_carriers();
     int ndrain_side_carriers  = Get_ndrain_side_carriers();
     double **p     = Get_p();
     int *valley = Get_valley();
     double mesh_size_x = Get_mesh_size_x();
     
     int Get_nx0(),Get_nx1(); 
     int nx0 = Get_nx0(); int nx1 = Get_nx1();
     int nx_max = nx1-nx0; // int nx_max = Get_nx_max();
    
     
// Thuc hien
//******************************************************************************
   // Buoc 1: Delete extra carriers at the source and drain contacts
     int iss_eli_update = 0, idd_eli_update = 0;
     double  x_position = 0.0;                
     int n,ix;
     //int count =0; // for checking
     int npt_source=0, npt_drain=0;
     int n_used = Get_n_used(); //printf("\n Truoc Check source and Drain contact n_used = %d",n_used);
     for(n = 1; n<=n_used; n++){ // Chay cho tung hat
        if(valley[n] !=9){// Chi chay cho cac hat co iv KHAC 9, Neu =9: xem ham delete_particle()
           //count +=1;
           x_position = p[n][2]; // lay vi tri hat thu n
           ix = (int)(x_position/mesh_size_x + 0.5);
           if(ix > nx_max) { 
               ix = nx_max; 
            }
           // Tai Source contact    
           if(ix == 0){ // Phan mep ben tay trai (cuc Source)
              if(nsource_side_carriers > npt_source )   
                   {
                      npt_source += 1; // tang them 1 hat
                   } // cho den khi nsource_side_carriers = npt_source
              else { // Qua so hat de tao charge neutrality. Can GHI NHO se Loai hat nay
                      valley[n] = 9;
                      iss_eli_update = Get_iss_eli();// 3 lenh tuong duong 1 lenh iss_eli += 1;
                      iss_eli_update += 1;  // Loai 1 hat ra
                      Set_iss_eli(iss_eli_update);
                   }
           }// End of if(ix ==0)
        
           // Tai Drain contact     
           if(ix ==nx_max){// Phan mep ben tay phai (cuc Drain)
              if(ndrain_side_carriers > npt_drain )
                   {
                     npt_drain += 1;
                   }
              else {
                     valley[n] = 9;
                     idd_eli_update = Get_idd_eli();
                     idd_eli_update += 1;
                     Set_idd_eli(idd_eli_update);
                
                   }
          }// End of if(ix ==nx_max)
       } // End of if(valley[n] !=9) 
     }// End of for(int n = 1; n<=n_used; n++)
     //printf("\n So hat co valley KHAC 9 la count= %d ",count);
     //printf("\n iss_eli = %d,Kieu 2 iss_eli = %d ", iss_eli_update,Get_iss_eli()); // Only for checking
     //printf("\n idd_eli = %d,Kieu 2 idd_eli = %d ", idd_eli_update,Get_idd_eli());
     // End of Buoc 1: Delete extra carriers at the source and drain contacts
//******************************************************************************

     // Buoc 2. Create carriers at the source and drain contacts
     int iss_cre_update=0, idd_cre_update=0; 
     int nele_diff = 0; // XAY RA khi sau khi quet het so hat ma nsource_side_carriers > npt_source
                        // So particle KHAC NHAU giua (npt_source va nsource_side_carriers) 
                        // cho Source contact, tuong tu cho Drain contact. ->
                        // need to create particles to make charge neutrality at the Souce
                        // and Drain ohmic contacts */
     int ne = 0; // Bien trung gian  
     
     // Tai Source contact
     ne = 0;
     ix = 0;
     if(nsource_side_carriers > npt_source){
        nele_diff = nsource_side_carriers - npt_source;
        ne = Get_n_used(); 
        while(nele_diff >= 1){ 
              ne = ne + 1;  // Add vao cac vi tri tiep theo va lan luot
              init_kspace(ne,ix);
              init_realspace(ne,ix);
              nele_diff = nele_diff - 1;
              iss_cre_update = Get_iss_cre();
              iss_cre_update += 1;
              Set_iss_cre(iss_cre_update);
         } // End of while
         
     Set_n_used(ne);// vi trong qua trinh tren "ne" da tang len, ta can update gia tri
    } // End of  if(nsource_side_carriers > npt_source)

    // Tai Drain contact
    ix = nx_max; 
    if(ndrain_side_carriers > npt_drain){
       nele_diff = ndrain_side_carriers - npt_drain;
       ne = Get_n_used();
       while(nele_diff >= 1){
             ne = ne + 1; // Add vao cac vi tri tiep theo va lan luot
             init_kspace(ne,ix);
             init_realspace(ne,ix);
             nele_diff = nele_diff - 1;
             idd_cre_update = Get_idd_cre();
             idd_cre_update += 1;
             Set_idd_cre(idd_cre_update);
       } // End of while
     
     Set_n_used(ne);        
    }// End of if(ndrain_side_carriers > npt_drain)
     //printf("\n iss_cre = %d,Kieu 2 iss_cre = %d ", iss_cre_update,Get_iss_cre()); // Only for checking
     //printf("\n idd_cre = %d,Kieu 2 idd_cre = %d ", idd_cre_update,Get_idd_cre());
     //printf("\n Sau Check source and Drain contact n_used = %d",Get_n_used());
     // Hien nhien so hat se tang len neu co su Creation vi cac hat bi delete do out hay eliminate
     // se thuc hien o ham detele_particle() goi sau ham check_source_drain_contacts()
return;
}
/* ************ End of check_source_drain_contacts() function ******************/
