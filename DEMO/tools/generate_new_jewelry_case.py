# modified the llps_file: add n_node
# added get_last_frame. n_??? 
# add new tpr for jewelry_case
import os
import sys


def generate_prep_file(protein_name):
    os.chdir(protein_name)
    os.chdir("prep")
    os.system("grep '^ATOM' " + protein_name + ".pdb > fix1.pdb")
    os.system("python ../../tools/change_ter_name.py")
    os.system("cp fix2.pdb fixed.pdb")
    os.system("smog2 -i fixed.pdb -CA")
    os.system("gmx editconf -f smog.gro -o chain_big_box.gro -box 100 100 100 -c")
    os.system("cp chain_big_box.gro chain_huge_box.gro")

    os.chdir("..")

    os.system("mkdir n_prep")
    os.system("mkdir b_prep")

    os.system("rm n_prep/*")
    os.system("rm b_prep/*")

    os.system("cp prep/* n_prep/")
    os.system("cp prep/* b_prep/")

    # n_prep
    os.chdir("n_prep")
    os.system("python ../../tools/gro2data_huge.py")
    os.system("python ../../tools/write_prep.py 300 100000 50 n")
    os.system("cp in.llps.raw in.llps")
    #os.system("python ../../tools/write_pbs.py 1 4 4 " + protein_name + "_np")
    os.system("python ../../tools/write_pbs.py 1 4 24 " + protein_name + "_np")
    os.system("cp ../../tools/set_param_Mpipi_nobonded.dat .")
    os.chdir("..")

    # b_prep
    os.chdir("b_prep")
    os.system("python ../../tools/gro2data_angular_dihedral_huge_rs.py")
    os.system("python ../../tools/write_prep.py 300 100000 50 y")
    os.system("cp in.llps.raw in.llps")
    #os.system("python ../../tools/write_pbs.py 1 4 4 " + protein_name + "_bp")
    os.system("python ../../tools/write_pbs.py 1 4 24 " + protein_name + "_bp")
    os.system("cp ../../tools/set_param_Mpipi_bonded_rs.dat .")
    os.chdir("..")

    os.chdir("..")

    return 0


def generate_npt_file(protein_name, n_chain, xy_axis, z_axis):
    os.chdir(protein_name)

    os.chdir("n_prep")
    #os.system("gmx grompp -f ../../tools/fake.mdp -c smog.gro -p smog.top -o fake.tpr -maxwarn 1")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p smog.top -o fake.tpr -maxwarn 1")
    #os.system("echo '0' | gmx trjconv -f md.xtc -s fake.tpr -o prep.gro -tu ns -b 103 -e 103 -pbc mol")
    os.system("echo '0' | gmx trjconv -f md.xtc -s fake.tpr -o prep.gro -tu ns -b 2500 -e 2500 -pbc mol")
    os.system("gmx insert-molecules -ci prep.gro -box " + xy_axis + " " + xy_axis + " " + z_axis + " -o disperse_chain_huge_box.gro -rot xyz -nmol " + n_chain + " -seed 1999 -radius 1.0 -try 50")
    os.chdir("..")

    os.chdir("b_prep")
    #os.system("gmx grompp -f ../../tools/fake.mdp -c smog.gro -p smog.top -o fake.tpr -maxwarn 1")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p smog.top -o fake.tpr -maxwarn 1")
    #os.system("echo '0' | gmx trjconv -f md.xtc -s fake.tpr -o prep.gro -tu ns -b 103 -e 103 -pbc mol")
    os.system("echo '0' | gmx trjconv -f md.xtc -s fake.tpr -o prep.gro -tu ns -b 2500 -e 2500 -pbc mol")
    os.system("gmx insert-molecules -ci prep.gro -box " + xy_axis + " " + xy_axis + " " + z_axis + " -o disperse_chain_huge_box.gro -rot xyz -nmol " + n_chain + " -seed 1999 -radius 1.0 -try 50")
    os.chdir("..")

    os.system("mkdir n_npt")
    os.system("mkdir b_npt")

    os.chdir("n_npt")
    os.system("cp ../n_prep/disperse_chain_huge_box.gro ./chain_huge_box.gro")
    os.system("sed '$ d' ../n_prep/smog.top > topol.top")
    os.system("echo '  Macromolecule   " + n_chain + "' >> topol.top")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p topol.top -o llps_fake.tpr -maxwarn 1")
    os.system("cp ../../tools/set_param_Mpipi_nobonded.dat .")
    os.system("python ../../tools/gro2data_huge.py")
    os.system("cp ../../tools/in.npt.n ./in.llps")
    os.system("python ../../tools/write_pbs.py 1 128 24 " + protein_name + "_n_npt")
    os.chdir("..")

    os.chdir("b_npt")
    os.system("cp ../b_prep/disperse_chain_huge_box.gro ./chain_huge_box.gro")
    os.system("sed '$ d' ../b_prep/smog.top > topol.top")
    os.system("echo '  Macromolecule   " + n_chain + "' >> topol.top")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p topol.top -o llps_fake.tpr -maxwarn 1")
    os.system("cp ../../tools/set_param_Mpipi_bonded_rs.dat .")
    os.system("python ../../tools/gro2data_angular_dihedral_huge_rs.py")
    os.system("cp ../../tools/in.npt.b ./in.llps")
    os.system("python ../../tools/write_pbs.py 1 128 24 " + protein_name + "_b_npt")
    os.chdir("..")

    os.chdir("..")

    return 0


def jewelry_case(protein_name, xy_axis, z_axis):
    os.chdir(protein_name)

    os.chdir("n_npt")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p topol.top -o llps_fake.tpr -maxwarn 1")
    os.system("echo '0' | gmx trjconv -f md.xtc -s llps_fake.tpr -o after_npt.gro -tu ns -b 30 -e 30 -pbc mol")
    #os.system("echo '0' | gmx trjconv -f md.xtc -s llps_fake.tpr -o after_npt.gro -tu ns -b 100 -e 100 -pbc mol")
    os.system("gmx editconf -f after_npt.gro -c -box " + xy_axis + " " + xy_axis + " " + z_axis + " -o after_npt_chain_huge_box.gro")
    os.system("mkdir ../n_jewelry_case/")
    os.system("cp after_npt_chain_huge_box.gro ../n_jewelry_case/chain_huge_box.gro")
    os.system("cp llps_fake.tpr ../n_jewelry_case/")
    os.system("cp set_param_Mpipi_nobonded.dat ../n_jewelry_case/")
    os.chdir("../n_jewelry_case/")
    os.system("python ../../tools/gro2data_huge.py")
    os.chdir("..")

    os.chdir("b_npt")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p topol.top -o llps_fake.tpr -maxwarn 1")
    os.system("echo '0' | gmx trjconv -f md.xtc -s llps_fake.tpr -o after_npt.gro -tu ns -b 30 -e 30 -pbc mol")
    #os.system("echo '0' | gmx trjconv -f md.xtc -s llps_fake.tpr -o after_npt.gro -tu ns -b 100 -e 100 -pbc mol")
    os.system("gmx editconf -f after_npt.gro -c -box " + xy_axis + " " + xy_axis + " " + z_axis + " -o after_npt_chain_huge_box.gro")
    os.system("mkdir ../b_jewelry_case/")
    os.system("cp after_npt_chain_huge_box.gro ../b_jewelry_case/chain_huge_box.gro")
    os.system("cp llps_fake.tpr ../b_jewelry_case/")
    os.system("cp set_param_Mpipi_bonded_rs.dat ../b_jewelry_case/")
    os.chdir("../b_jewelry_case/")
    os.system("python ../../tools/gro2data_angular_dihedral_huge_rs.py")
    os.chdir("..")

    os.chdir("..")

    return 0


def self_assembly_case(protein_name, n_chain, xy_axis, z_axis):
    os.chdir(protein_name)

    os.chdir("n_prep")
    os.system("gmx grompp -f ../../tools/fake.mdp -c smog.gro -p smog.top -o fake.tpr -maxwarn 1")
    os.system("echo '0' | gmx trjconv -f md.xtc -s fake.tpr -o prep.gro -tu ns -b 103 -e 103 -pbc mol")
    os.system("gmx insert-molecules -ci prep.gro -box " + xy_axis + " " + xy_axis + " " + z_axis + " -o disperse_chain_huge_box.gro -rot xyz -nmol " + n_chain + " -seed 1999 -radius 1.0 -try 50")
    os.chdir("..")

    os.chdir("b_prep")
    os.system("gmx grompp -f ../../tools/fake.mdp -c smog.gro -p smog.top -o fake.tpr -maxwarn 1")
    os.system("echo '0' | gmx trjconv -f md.xtc -s fake.tpr -o prep.gro -tu ns -b 103 -e 103 -pbc mol")
    os.system("gmx insert-molecules -ci prep.gro -box " + xy_axis + " " + xy_axis + " " + z_axis + " -o disperse_chain_huge_box.gro -rot xyz -nmol " + n_chain + " -seed 1999 -radius 1.0 -try 50")
    os.chdir("..")

    os.system("mkdir n_self_assembly_case/")
    os.system("mkdir b_self_assembly_case/")

    os.chdir("n_self_assembly_case")
    os.system("cp ../n_prep/disperse_chain_huge_box.gro chain_huge_box.gro")
    os.system("sed '$ d' ../n_prep/smog.top > topol.top")
    os.system("echo '  Macromolecule   " + n_chain + "' >> topol.top")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p topol.top -o llps_fake.tpr -maxwarn 1")
    os.system("cp ../../tools/set_param_Mpipi_nobonded.dat .")
    os.system("python ../../tools/gro2data_huge.py")

    os.chdir("..")

    os.chdir("b_self_assembly_case")
    os.system("cp ../b_prep/disperse_chain_huge_box.gro chain_huge_box.gro")
    os.system("sed '$ d' ../b_prep/smog.top > topol.top")
    os.system("echo '  Macromolecule   " + n_chain + "' >> topol.top")
    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p topol.top -o llps_fake.tpr -maxwarn 1")
    os.system("cp ../../tools/set_param_Mpipi_bonded_rs.dat .")
    os.system("python ../../tools/gro2data_angular_dihedral_huge_rs.py")
    os.chdir("..")

    os.chdir("..")

    return 0


def llps_file(protein_name, simulation_temperature, n_node): 
    os.chdir(protein_name)
    os.system("mkdir n_" + simulation_temperature)
    os.system("mkdir b_" + simulation_temperature)

    os.system("cp n_jewelry_case/* n_" + simulation_temperature) 
    os.chdir("n_" + simulation_temperature)
    os.system("python ../../tools/write_in_file.py " + simulation_temperature + " 100000 50 n")
    os.system("cp in.llps.raw in.llps")
    os.system("python ../../tools/write_pbs.py " + n_node +" 128 24 " + protein_name + "n" + simulation_temperature)
    os.chdir("..")

    os.system("cp b_jewelry_case/* b_" + simulation_temperature)
    os.chdir("b_" + simulation_temperature)
    os.system("python ../../tools/write_in_file.py " + simulation_temperature + " 100000 50 y")
    os.system("cp in.llps.raw in.llps")
    os.system("python ../../tools/write_pbs.py " + n_node + " 128 24 " + protein_name + "b" + simulation_temperature)
    os.chdir("..")

    os.chdir("..")

    return 0


def self_assembly_file(protein_name, simulation_temperature): 
    os.chdir(protein_name)
    os.system("mkdir n_" + simulation_temperature)
    os.system("mkdir b_" + simulation_temperature)

    os.system("cp n_self_assembly_case/* n_" + simulation_temperature) 
    os.chdir("n_" + simulation_temperature)
    os.system("python ../../tools/write_self_assembly_in_file.py " + simulation_temperature + " 100000 50 n")
    os.system("cp in.llps.raw in.llps")
    os.system("python ../../tools/write_pbs.py 1 128 24 " + protein_name + "n" + simulation_temperature)
    os.chdir("..")

    os.system("cp b_self_assembly_case/* b_" + simulation_temperature)
    os.chdir("b_" + simulation_temperature)
    os.system("python ../../tools/write_self_assembly_in_file.py " + simulation_temperature + " 100000 50 y")
    os.system("cp in.llps.raw in.llps")
    os.system("python ../../tools/write_pbs.py 1 128 24 " + protein_name + "b" + simulation_temperature)
    os.chdir("..")

    os.chdir("..")

    return 0


def ex_file(protein_name, simulation_temperature, index_ex, n_node):
    index_last = int(index_ex) - 1

    os.chdir(protein_name)

    os.system("mkdir " + "ne" + index_ex + "_" + simulation_temperature)
    os.system("mkdir " + "be" + index_ex + "_" + simulation_temperature)

    os.system("rm ne" + index_ex + "_" + simulation_temperature + "/*")
    os.system("rm be" + index_ex + "_" + simulation_temperature + "/*")
    os.system("cp n_jewelry_case/* ne" + index_ex + "_" + simulation_temperature)
    os.system("cp b_jewelry_case/* be" + index_ex + "_" + simulation_temperature)

    # ne1_300
    os.chdir("ne" + index_ex + "_" + simulation_temperature)
    if index_last == 0:
        os.system("cp ../n_" + simulation_temperature + "/tmp.restart.* .")
    else:
        os.system("cp ../ne" + str(index_last) + "_" + simulation_temperature + "/tmp.restart.* .")

    os.system("python ../../tools/write_ex_in_file.py " + simulation_temperature + " 100000 50 n")
    os.system("cp in.llps.ex.raw in.llps.ex")
    os.system("python ../../tools/write_ex_pbs.py " + n_node + " 128 24 " + protein_name + "n" + index_ex + simulation_temperature)
    os.chdir("..")

    # be1_300
    os.chdir("be" + index_ex + "_" + simulation_temperature)
    if index_last == 0:
        os.system("cp ../b_" + simulation_temperature + "/tmp.restart.* .")
    else:
        os.system("cp ../be" + str(index_last) + "_" + simulation_temperature + "/tmp.restart.* .")

    os.system("python ../../tools/write_ex_in_file.py " + simulation_temperature + " 100000 50 y")
    os.system("cp in.llps.ex.raw in.llps.ex")
    os.system("python ../../tools/write_ex_pbs.py " + n_node + " 128 24 " + protein_name + "b" + index_ex + simulation_temperature)
    os.chdir("..")

    os.chdir("..")

    return 0


def self_assembly_ex_file(protein_name, simulation_temperature, index_ex, n_node):
    index_last = int(index_ex) - 1

    os.chdir(protein_name)

    os.system("mkdir " + "ne" + index_ex + "_" + simulation_temperature)
    os.system("mkdir " + "be" + index_ex + "_" + simulation_temperature)

    os.system("rm ne" + index_ex + "_" + simulation_temperature + "/*")
    os.system("rm be" + index_ex + "_" + simulation_temperature + "/*")
    os.system("cp n_self_assembly_case/* ne" + index_ex + "_" + simulation_temperature)
    os.system("cp b_self_assembly_case/* be" + index_ex + "_" + simulation_temperature)

    # ne1_300
    os.chdir("ne" + index_ex + "_" + simulation_temperature)
    if index_last == 0:
        os.system("cp ../n_" + simulation_temperature + "/tmp.restart.* .")
    else:
        os.system("cp ../ne" + str(index_last) + "_" + simulation_temperature + "/tmp.restart.* .")

    os.system("python ../../tools/write_self_assembly_ex_in_file.py " + simulation_temperature + " 100000 50 n")
    os.system("cp in.llps.ex.raw in.llps.ex")
    os.system("python ../../tools/write_ex_pbs.py " + n_node + " 128 24 " + protein_name + "n" + index_ex + simulation_temperature)
    os.chdir("..")

    # be1_300
    os.chdir("be" + index_ex + "_" + simulation_temperature)
    if index_last == 0:
        os.system("cp ../b_" + simulation_temperature + "/tmp.restart.* .")
    else:
        os.system("cp ../be" + str(index_last) + "_" + simulation_temperature + "/tmp.restart.* .")

    os.system("python ../../tools/write_self_assembly_ex_in_file.py " + simulation_temperature + " 100000 50 y")
    os.system("cp in.llps.ex.raw in.llps.ex")
    os.system("python ../../tools/write_ex_pbs.py " + n_node + " 128 24 " + protein_name + "b" + index_ex + simulation_temperature)
    os.chdir("..")

    os.chdir("..")

    return 0


def run_prep(protein_name):
    os.chdir(protein_name)

    os.chdir("n_prep")
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("b_prep")
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("..")

    return 0


def run_npt(protein_name):
    os.chdir(protein_name)

    os.chdir("n_npt")
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("b_npt")
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("..")

    return 0


def run_llps(protein_name, simulation_temperature):
    os.chdir(protein_name)

    os.chdir("n_" + simulation_temperature)
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("b_" + simulation_temperature)
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("..")
    return 0


def run_ex(protein_name, simulation_temperature, index_ex):
    os.chdir(protein_name)

    os.chdir("ne" + index_ex + "_" + simulation_temperature)
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("be" + index_ex + "_" + simulation_temperature)
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("..")

    return 0


def write_single_file(protein_name, simulation_temperature):
    # os.chdir("sasbdb_" + protein_name)
    os.chdir(protein_name)

    os.system("mkdir b_" + simulation_temperature)

    os.chdir("b_" + simulation_temperature)
    os.system("cp ../prep/chain_huge_box.gro .")
    os.system("cp ../prep/fake.tpr .")
    os.system("cp ../../tools/set_param_Mpipi_bonded_rs.dat .")
    os.system("python ../../tools/gro2data_angular_dihedral_huge_rs.py")
    os.system("python ../../tools/write_in_single.py " + simulation_temperature + " 100000 y")
    os.system("cp in.llps.raw in.llps")
    os.system("python ../../tools/write_pbs.py 1 24 24 " + protein_name + "_b_" + simulation_temperature)
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.system("mkdir n_" + simulation_temperature)

    os.chdir("n_" + simulation_temperature)
    os.system("cp ../prep/chain_huge_box.gro .")
    os.system("cp ../prep/fake.tpr .")
    os.system("cp ../../tools/set_param_Mpipi_nobonded.dat .")
    os.system("python ../../tools/gro2data_huge.py")
    os.system("python ../../tools/write_in_single.py " + simulation_temperature + " 100000 n")
    os.system("cp in.llps.raw in.llps")
    os.system("python ../../tools/write_pbs.py 1 24 24 " + protein_name + "_n_" + simulation_temperature)
    os.system("qsub sub.pbs")
    os.chdir("..")

    os.chdir("..")

    return 0


def last_frame(): 
    raw_out=os.popen("ls tmp.restart.*")
    t_r=raw_out.read()
    tt_r=t_r.replace("tmp.restart.","")
    ttt_r=tt_r.split()
    l_num=[]
    for i in ttt_r: 
        l_num.append(int(i))

    tmp_max=max(l_num)
    gmx_num = tmp_max//100

    with open ("frame_info.txt","w+") as f_w: 
        f_w.writelines(str(gmx_num))

    os.system("gmx grompp -f ../../tools/fake.mdp -c chain_huge_box.gro -p ../n_npt/topol.top -o llps_fake.tpr -maxwarn 1")
    os.system("echo '0' | gmx trjconv -f md.xtc -s llps_fake.tpr -o last.gro -b " + str(gmx_num) + " -e " + str(gmx_num) + " -pbc mol &")

    return 0


def get_last_frame(protein_name, simulation_temperature, ex_index): 
    os.chdir(protein_name)

    #os.chdir("ne" + ex_index + "_" + simulation_temperature)
    os.chdir("n_" + simulation_temperature)
    last_frame()
    os.chdir("..")

    #os.chdir("be" + ex_index + "_" + simulation_temperature)
    os.chdir("b_" + simulation_temperature)
    last_frame()
    os.chdir("..")
    
    os.chdir("..")
    return 0


if __name__ == "__main__":
    #     test()
    fn_protein_temperature = sys.argv[1]
    command = sys.argv[2]

    protein = []
    temperature = []
    n_chain = []
    xy_axis = []
    z_axis = []
    index_ex = []
    n_node = []
    with open(fn_protein_temperature, "r") as f_r1:
        t = f_r1.readlines()

    #    list_length = len(t)

    for i in t:
        ti = i.split()
        print(ti)
        if not i.startswith("#"):
            if command == "prep_file" or command == "run_prep" or command == "run_npt":
                protein.append(ti[0])
            elif command == "npt_file" or command == "jewelry_case" or command == "self_assembly_case":
                protein.append(ti[0])
                n_chain.append(ti[1])
                xy_axis.append(ti[2])
                z_axis.append(ti[3])
            #elif command == "llps_file" or command == "self_assembly_file": # or command == "run_llps":
            elif command == command == "self_assembly_file": # or command == "run_llps":
                protein.append(ti[0])
                temperature.append(ti[1])
            elif command == "llps_file" or "run_llps" or "ex_file" or "self_assembly_ex_file" or "last_frame":
                protein.append(ti[0])
                temperature.append(ti[1])
                index_ex.append(ti[2])
                n_node.append(ti[3])
            elif command == command == "run_ex":
                protein.append(ti[0])
                temperature.append(ti[1])
                index_ex.append(ti[2])
            else:
                print("No such command! Please python generate_new_jewelry_case.py file_name "
                      "prep_file/run_prep/npt_file/run_npt/jewelry_case/self_assembly_case/llps_file/self_assembly_file/run_llps/ex_file/self_assembly_ex_file/run_ex/last_frame")
        else:
            pass

    list_length = len(protein)

    for i in range(list_length):
        if command == "prep_file":
            generate_prep_file(protein[i])
        elif command == "run_prep":
            run_prep(protein[i])
        elif command == "run_npt":
            run_npt(protein[i])
        elif command == "npt_file":
            generate_npt_file(protein[i], n_chain[i], xy_axis[i], z_axis[i])
        elif command == "jewelry_case": 
            jewelry_case(protein[i], xy_axis[i], z_axis[i])
        elif command == "self_assembly_case": 
            self_assembly_case(protein[i], n_chain[i], xy_axis[i], z_axis[i])
        elif command == "llps_file":
            llps_file(protein[i], temperature[i], n_node[i])
            print(temperature[i], n_node[i])
        elif command == "self_assembly_file":
            self_assembly_file(protein[i], temperature[i])
            print(temperature[i])
        elif command == "run_llps": 
            run_llps(protein[i], temperature[i])
            print(temperature[i])
        elif command == "ex_file":
            ex_file(protein[i], temperature[i], index_ex[i], n_node[i])
            print(protein[i], temperature[i], index_ex[i], n_node[i])
        elif command == "self_assembly_ex_file":
            self_assembly_ex_file(protein[i], temperature[i], index_ex[i], n_node[i])
            print(protein[i], temperature[i], index_ex[i], n_node[i])
        elif command == "run_ex":
            run_ex(protein[i], temperature[i], index_ex[i])
            print(protein[i], temperature[i], index_ex[i])
        elif command == "last_frame": 
            get_last_frame(protein[i],temperature[i],index_ex[i])
        else:
            pass
        print(protein[i])
