# temperature allowing float
from math import ceil
import sys

print("python write_in_file.py <temperature> <T_damp> <step_length> <bonded(y)>")

T = float(sys.argv[1])
T_damp = sys.argv[2]
step_length = float(sys.argv[3])
with_bonded = sys.argv[4]
# T=600
# step_length=50
rds = 123456789
start_t = 50
n = ceil((T - start_t) / step_length)
version_name = "topaz London Blue"

with open("in.llps.raw", "w+") as f_w:
    f_w.writelines("# version " + version_name + "\n")
    f_w.writelines("variable temperature equal " + str(T) + "\n")
    f_w.writelines("variable randomSeed equal " + str(rds) + "\n")
    f_w.writelines("\n")

    f_w.writelines("# Initialize parameters from units and box type\n")
    f_w.writelines("units\t\t\treal\ndimension\t\t3\nboundary\t\tp p p\natom_style\t\tfull\n\n")

    f_w.writelines("# Initial pair style and bond parameters\n")
    f_w.writelines("bond_style\t\tharmonic\ndielectric\t\t80.0\n\n")

    f_w.writelines("# Create new simulation box or load existing box\n")
    if with_bonded == "y":
        f_w.writelines("fix\t\t\t\tCG_CMAP all cmap ../../table/cg.cmap\n")
        f_w.writelines("read_data\t\tdata.name fix CG_CMAP crossterm CMAP\n")
        f_w.writelines("fix_modify\t\tCG_CMAP energy yes\n\n")
        f_w.writelines("# Set pair parameters\n")
        f_w.writelines("include\t\t\t./set_param_Mpipi_bonded_rs.dat\n\n")
        f_w.writelines("special_bonds\tlj/coul 0 0 0\n\n")
    else:
        f_w.writelines("read_data\t\tdata.name\n\n")
        f_w.writelines("# Set pair parameters\n")
        f_w.writelines("include\t\t./set_param_Mpipi_nobonded.dat\n\n")
        f_w.writelines("special_bonds\tlj/coul 0 1 1\n\n")

    f_w.writelines("# neighbour list settings\n ")
    f_w.writelines("neighbor\t\t3.5 multi\n\n")

    f_w.writelines("# Timestep and computational parameters\n")
    f_w.writelines("comm_style\t\ttiled\n")
    f_w.writelines("timestep\t\t10\n\n")

    f_w.writelines("# relaxation by energy minimization\n")

    f_w.writelines("min_style\t\tfire\n")
    for i in [0.0000001, 0.00001, 0.1, 10]:
        f_w.writelines("timestep\t\t" + str(i) + "\nminimize\t\t0.0 1.0e-8 1000 100000\n")
    f_w.writelines("\n")

    f_w.writelines("neigh_modify\tevery 10 delay 0\n\n")

    f_w.writelines("velocity\t\tall create " + str(start_t) + " ${randomSeed}\n\n")

    f_w.writelines("fix\t\t\t\tfxvne all nve\n\n")

    f_w.writelines("fix\t\t\t\tfxbal all balance 1000 1.1 rcb\n\n")

    f_w.writelines("fix\t\t\t\tcom_motion_remove all momentum 1000 linear 1 1 1 angular rescale\n\n")

    f_w.writelines("restart\t\t\t10000000 ./tmp.restart\n\n")

    f_w.writelines("compute\t\t\tmyChunk2 all chunk/atom bin/1d z lower 0.0005 units reduced\n")
    f_w.writelines("fix\t\t\t\tdensityRegions2 all ave/chunk 1000 10 20000 myChunk2 density/mass file "
                   "densities_chunked2.dat\n\n")

    f_w.writelines("thermo\t\t\t1000\n")
    if with_bonded == "y":
        f_w.writelines(
            "thermo_style\tcustom step fmax ebond eangle edihed f_CG_CMAP evdwl ecoul pe ke etotal temp press\n")
    else:
        f_w.writelines("thermo_style\tcustom step fmax ebond eangle edihed evdwl ecoul pe ke etotal temp press\n")
    f_w.writelines("thermo_modify\tflush yes\n\n")

    f_w.writelines("dump\t\t\tdumptraj all xtc 1000 ./md.xtc\n")
    f_w.writelines("log\t\t\t\t./log.lammps.cg\n\n")

    f_w.writelines("# warming up\n# " + str(2 * n) + " ns.\n")
    for i in range(int(n - 1)):
        f_w.writelines(
            "fix\t\t\t\tfxl" + str(2 * i + 1) + " all langevin " + str(start_t + i * step_length) + " " + str(
                start_t + i * step_length) + " " + T_damp + " ${randomSeed}\n\n")
        f_w.writelines("run\t\t\t\t100000\n\n")
        f_w.writelines("unfix\t\t\tfxl" + str(2 * i + 1) + "\n")
        f_w.writelines(
            "fix\t\t\t\tfxl" + str(2 * i + 2) + " all langevin " + str(start_t + i * step_length) + " " + str(
                start_t + (i + 1) * step_length) + " " + T_damp + " ${randomSeed}\n\n")
        f_w.writelines("run\t\t\t\t100000\n\n")
        f_w.writelines("unfix\t\t\tfxl" + str(2 * i + 2) + "\n")

    f_w.writelines(
        "fix\t\t\t\tfxl" + str(2 * n - 1) + " all langevin " + str(start_t + (n - 1) * step_length) + " " + str(
            start_t + (n - 1) * step_length) + " " + T_damp + " ${randomSeed}\n\n")
    f_w.writelines("run\t\t\t\t100000\n\n")
    f_w.writelines("unfix\t\t\tfxl" + str(2 * n - 1) + "\n")
    f_w.writelines("fix\t\t\t\tfxl" + str(2 * n) + " all langevin " + str(
        start_t + (n - 1) * step_length) + " ${temperature} " + T_damp + " ${randomSeed}\n\n")
    f_w.writelines("run\t\t\t\t100000\n\n")
    f_w.writelines("unfix\t\t\tfxl" + str(2 * n) + "\n")
    f_w.writelines(
        "fix\t\t\t\tfxl" + str(2 * n + 1) + " all langevin " + str(
            T) + " ${temperature} " + T_damp + " ${randomSeed}\n\n")

    f_w.writelines("# 1ns\n# eq \n")
    for i in range(5):
        f_w.writelines("run\t\t\t\t20000\n")
    f_w.writelines("\n")
    f_w.writelines("# 100ns\n# Run simulation\n")
    for i in range(5):
        f_w.writelines("run\t\t\t\t2000000\n")
    f_w.writelines("\n")
    f_w.writelines("# 3000ns. i.e. 3us\n")
    for i in range(3):
        for j in range(5):
            f_w.writelines("run\t\t\t\t20000000\n")
            f_w.writelines("write_data\t\texam_"+str(i*5+j)+".dat nocoeff\n\n")
        f_w.writelines("\n")

    f_w.writelines("write_data\t\tfinal-structure.dat nocoeff\n")
