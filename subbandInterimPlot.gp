#Tu gnuplot dung lenh load 'subbandPlot.gp' la chay. Muon ve cai nao thi BO DAU # o cai do

#Plot TAT CA Valley va Subbands. DUNG cho SO SUBBAND =4
#plot 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:2, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:3, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:4,\
     'subbandInterim-Vd0.50-Vg0.50.dat' using 1:5, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:6, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:7, \
     'subbandInterim-Vd0.50-Vg0.50.dat' using 1:8, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:9, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:10,\
     'subbandInterim-Vd0.50-Vg0.50.dat' using 1:11,'subbandInterim-Vd0.50-Vg0.50.dat' using 1:12,'subbandInterim-Vd0.50-Vg0.50.dat' using 1:13 \
     ,'potential_x_at_midpoint_yz_in2DSchrodingerInterim-Vd0.50-Vg0.50.dat'
 
#Plot Valley 2 va 3 cho Tat ca Subbands. DUNG cho SO SUBBAND =4
plot  'subbandInterim-Vd0.50-Vg0.50.dat' using 1:6, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:7, \
     'subbandInterim-Vd0.50-Vg0.50.dat' using 1:8, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:9, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:10,\
     'subbandInterim-Vd0.50-Vg0.50.dat' using 1:11,'subbandInterim-Vd0.50-Vg0.50.dat' using 1:12,'subbandInterim-Vd0.50-Vg0.50.dat' using 1:13 \
     ,'potential_x_at_midpoint_yz_in2DSchrodingerInterim-Vd0.50-Vg0.50.dat'

#Plot for all subbands in Valley 1
plot 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:2, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:3, 'subbandInterim-Vd0.50-Vg0.50.dat' using 1:4,\
     'subbandInterim-Vd0.50-Vg0.50.dat' using 1:5, 'potential_x_at_midpoint_yz_in2DSchrodingerInterim-Vd0.50-Vg0.50.dat'

