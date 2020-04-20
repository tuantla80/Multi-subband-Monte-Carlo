# Set style line la 20 thi no nhin bar se ro, neu lon hon vi du 25 tro len thi nhin khong ro
# Cho truong hop mo full man hinh
set style line  1 lt 1 lw 20

# Set range  truong hop le EletronInitEnergy thi range cua energy nhu the la vua dep
#set xr [0.0:0.12]

#set multiplot plot de len nhau

# Ve Electron Distribution theo Energy
# Ten file co the can duoc them vao neu Vd va Vg la khac nhau

#Cho Valley 1
#plot 'ElectronInitDistributionEnergy-Vd0.50-Vg0.50.dat'u 1:2 with imp ls 1
#set title "Initial Electron Distribution in Valley 1"

#Cho Valley 2
#plot 'ElectronInitDistributionEnergy-Vd0.50-Vg0.50.dat'u 1:3 with imp ls 1
#set title "Initial Electron Distribution in Valley 2"

#Cho Valley 3
#plot 'ElectronInitDistributionEnergy-Vd0.50-Vg0.50.dat'u 1:4 with imp ls 1
#set title "Initial Electron Distribution in Valley 3"

#Cho All Valleys
plot 'ElectronInitDistributionEnergy-Vd0.50-Vg0.50.dat'u 1:5 with imp ls 1
#set title "Initial Electron Distribution in All Valleys"



set xlabel "Energy [eV]"
set ylabel "Number of Electrons"


