# Set range  truong hop le EletronInitSubband thi range 1 den NSELECT
# Ve Electron Distribution theo Subband
# Ten file co the can duoc them vao neu Vd va Vg la khac nhau

#Cho Valley 1
#plot 'ElectronInitDistributionSubbands-Vd0.50-Vg0.50.dat' u 1:2
#with imp ls 1
#set title "Initial Electron Distribution in Valley 1"

#Cho Valley 2
#plot 'ElectronInitDistributionSubbands-Vd0.50-Vg0.50.dat'u 1:3
#with imp ls 1
#set title "Initial Electron Distribution in Valley 2"

#Cho Valley 3
#plot 'ElectronInitDistributionSubbands-Vd0.50-Vg0.50.dat'u 1:4
#with imp ls 1
#set title "Initial Electron Distribution in Valley 3"

#Cho All Valleys
plot 'ElectronInterimDistributionEnergy-Vd0.50-Vg0.50.dat' u 1:2,'ElectronInterimDistributionEnergy-Vd0.50-Vg0.50.dat'u 1:3, 'ElectronInterimDistributionEnergy-Vd0.50-Vg0.50.dat'u 1:4
#set title "Initial Electron Distribution in All Valleys"

set xlabel "Subband Index"
set ylabel "Number_of_Electrons"


