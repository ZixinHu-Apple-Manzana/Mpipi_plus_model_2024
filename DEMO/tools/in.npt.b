# 5ns is sufficient already
variable temperature equal 150
variable randomSeed equal 123456789

# Initialise parameters for units and box type
units       real
dimension   3
boundary    p p p
atom_style  full

# Initital pair style and bond parameters
bond_style  harmonic

dielectric  80.0

# Create new simulation box or load existing box

fix		CG_CMAP all cmap ../../table/cg.cmap
read_data       data.name fix CG_CMAP crossterm CMAP
fix_modify	CG_CMAP energy yes

# Set pair parameters
# include         set_pair_param-Mpipi.dat
include         ./set_param_Mpipi_bonded_rs.dat
#include		./set_param_Mpipi_nobonded.dat

#reset_atom_ids
velocity        all create ${temperature} ${randomSeed}

special_bonds	lj/coul 0 1 1

# neighbour list settings
neighbor  3.5 multi

# Timestep and computational parameters
comm_style      tiled
timestep        10

# Relax sim by performing minimisation (this is helpful if we change the potential or sequence)
min_style       fire
timestep        0.0000001
minimize 0.0 1.0e-8 1000 100000
timestep        0.00001
minimize 0.0 1.0e-8 1000 100000
timestep        0.1
minimize 0.0 1.0e-8 1000 100000
timestep        10
minimize 0.0 1.0e-8 1000 100000


neigh_modify    every 10 delay 0

fix             fxnve   all nve
fix             fxlange all langevin ${temperature} ${temperature} 100000.0 ${randomSeed}
#fix             fxnvt   all temp/rescale 10000 ${temperature} ${temperature} 1 1
#fix             fxnvt   all nvt temp ${temperature} ${temperature} 1000.0
#fix             fxt     all temp/berendsen ${temperature} ${temperature} 100000.0
#fix		fxp 	all press/berendsen z 1 1 10000
fix		fxp 	all press/berendsen iso 1 1 10000

fix             fxbal  all balance 1000 1.1 rcb

fix             com_motion_remove all momentum 10000 linear 1 1 1 angular rescale

restart         10000000 ./tmp.restart

#compute         myChunk2 all chunk/atom bin/1d z lower 0.005 units reduced
#fix             densityRegions2 all ave/chunk 1000 80 100000 myChunk2 density/mass file densities_chunked2.dat

#fix fixCOM all recenter INIT INIT INIT  # this helps to prevent drift of the interfaces and makes it easier to average over many measurements

# Thermo output settings
thermo          1000
#thermo_style    custom step pe ecoul ke temp press density
#thermo_style	custom step lz pzz fmax ebond eangle edihed f_CG_CMAP evdwl ecoul pe ke etotal temp press
thermo_style	custom step lz pxx pyy pzz fmax ebond eangle edihed f_CG_CMAP evdwl ecoul pe ke etotal temp press
#thermo_style	custom step lz pzz fmax ebond eangle edihed evdwl ecoul pe ke etotal temp press

#thermo_style	custom step fmax ebond eangle edihed evdwl ecoul pe ke etotal temp press
thermo_modify   flush yes

#dump            1 all custom 200000 result.lammpstrj id mol type q x y z
dump		dumptraj all xtc 1000 ./md.xtc 

log     ./log.lammps.cg

# 10ns
# Equilibration
run		200000
run		200000
run		200000
run		200000
run		200000

#run         600000

# 20ns
# Run simulation
run		2000000
#run		2000000
#run		2000000
#run		2000000
#run		2000000

# write final configuration
write_data      final-structure.dat nocoeff

