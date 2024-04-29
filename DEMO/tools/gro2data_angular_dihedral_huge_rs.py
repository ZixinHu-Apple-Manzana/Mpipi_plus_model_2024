from datetime import datetime

gro_fn = "chain_huge_box.gro"
protein_name = "name"


# pair_coeff_fn = "pair_coeff.txt"
# cg_model = "HPS_M3"

def hello():
    s = "# This script converts .gro file to data. file and pair coefficient file that can be read by lammps. \n"
    return s


class Residue:  # number, aaType, str_x, str_y, str_z
    residue_number = 0
    aaType = 0
    # in A
    x = "0.00"
    y = "0.00"
    z = "0.00"

    def __init__(self, residue_number, aaType, str_x, str_y, str_z):
        self.residue_number = residue_number
        self.aaType = aaType
        self.x = str_x
        self.y = str_y
        self.z = str_z


# class AminoAcidType: # name, fullname, hydrophobicity, radius, index, mass, charge
#     name = ""
#     fullname = ""
#     hydrophobicity = 0.0
#     radius = 0.0
#     index = 0
#     mass = 0.00
#     charge = 0.0
#
#     def __init__(self, name, fullname, hydrophobicity, radius, index, mass, charge):
#         self.name = name
#         self.fullname = fullname
#         self.hydrophobicity = hydrophobicity
#         self.radius = radius
#         self.index = index
#         self.mass = mass
#         self.charge = charge

class AminoAcidType:  # name, fullname, hydrophobicity, radius, index, mass, charge
    name = ""
    fullname = ""
    hydrophobicity = 0.0
    radius = 0.0
    index = 0
    mass = 0.00
    charge = 0.0
    angular_pool1 = 0;
    angular_pool2 = 0;
    angular_pool3 = 0;

    def __init__(self, name, fullname, index, mass, charge, angular_pool1, angular_pool2, angular_pool3):
        self.name = name
        self.fullname = fullname
        # self.hydrophobicity = hydrophobicity
        # self.radius = radius
        self.index = index
        self.mass = mass
        self.charge = charge
        self.angular_pool1 = angular_pool1;
        self.angular_pool2 = angular_pool2;
        self.angular_pool3 = angular_pool3;


# M = AminoAcidType("NB_M", "MET", 0, 131.19, 0.0)
# G = AminoAcidType("NB_G", "GLY", 1, 57.05, 0.0)
# K = AminoAcidType("NB_K", "LYS", 2, 128.17, 1.0)
# T = AminoAcidType("NB_T", "THR", 3, 101.10, 0.0)
# R = AminoAcidType("NB_R", "ARG", 4, 156.18, 1.0)
# A = AminoAcidType("NB_A", "ALA", 5, 71.07, 0.0)
# D = AminoAcidType("NB_D", "ASP", 6, 115.08, -1.0)
# E = AminoAcidType("NB_E", "GLU", 7, 129.10, -1.0)
# Y = AminoAcidType("NB_Y", "TYR", 8, 163.17, 0.0)
# V = AminoAcidType("NB_V", "VAL", 9, 99.13, 0.0)
# L = AminoAcidType("NB_L", "LEU", 10, 113.15, 0.0)
# Q = AminoAcidType("NB_Q", "GLN", 11, 128.13, 0.0)
# W = AminoAcidType("NB_W", "TRP", 12, 186.07, 0.0)
# F = AminoAcidType("NB_F", "PHE", 13, 147.17, 0.0)
# S = AminoAcidType("NB_S", "SER", 14, 87.07, 0.0)
# H = AminoAcidType("NB_H", "HIS", 15, 137.14, 0.5)
# N = AminoAcidType("NB_N", "ASN", 16, 114.10, 0.0)
# P = AminoAcidType("NB_P", "PRO", 17, 97.11, 0.0)
# C = AminoAcidType("NB_C", "CYS", 18, 103.13, 0.0)
# I = AminoAcidType("NB_I", "ILE", 19, 113.15, 0.0)
M = AminoAcidType("NB_M", "MET", 0,  131.199997,    0.0,    1,  8,  1)
G = AminoAcidType("NB_G", "GLY", 1,  57.049999,     0.0,    4,  13, 1)
K = AminoAcidType("NB_K", "LYS", 2,  128.199997,    0.75,   1,  7,  1)
T = AminoAcidType("NB_T", "THR", 3,  101.099998,    0.0,    3,  4,  1)
R = AminoAcidType("NB_R", "ARG", 4,  156.199997,    0.75,   1,  7,  1)
A = AminoAcidType("NB_A", "ALA", 5,  71.080002,     0.0,    2,  5,  1)
D = AminoAcidType("NB_D", "ASP", 6,  115.099998,    -0.75,  6,  11, 1)
E = AminoAcidType("NB_E", "GLU", 7,  129.100006,    -0.75,  2,  6,  1)
Y = AminoAcidType("NB_Y", "TYR", 8,  163.199997,    0.0,    1,  1,  1)
V = AminoAcidType("NB_V", "VAL", 9,  99.070000,     0.0,    5,  12, 2)
L = AminoAcidType("NB_L", "LEU", 10, 113.199997,    0.0,    1,  9,  1)
Q = AminoAcidType("NB_Q", "GLN", 11, 128.100006,    0.0,    1,  7,  1)
W = AminoAcidType("NB_W", "TRP", 12, 186.199997,    0.0,    1,  8,  1)
F = AminoAcidType("NB_F", "PHE", 13, 147.199997,    0.0,    1,  1,  1)
S = AminoAcidType("NB_S", "SER", 14, 87.080002,     0.0,    7,  10, 1)
H = AminoAcidType("NB_H", "HIS", 15, 137.100006,    0.375,  2,  2,  1)
N = AminoAcidType("NB_N", "ASN", 16, 114.099998,    0.0,    2,  11, 1)
P = AminoAcidType("NB_P", "PRO", 17, 97.120003,     0.0,    2,  14, 3)
C = AminoAcidType("NB_C", "CYS", 18, 103.099998,    0.0,    2,  3,  1)
I = AminoAcidType("NB_I", "ILE", 19, 113.199997,    0.0,    5,  12, 2)


# C = AminoAcidType("NB_C", "CYS", 0.594595, 5.48, 0, 103.13, 0.0)
# M = AminoAcidType("NB_M", "MET", 0.837838, 6.18, 1, 131.19, 0.0)
# F = AminoAcidType("NB_F", "PHE", 1.0, 6.36, 2, 147.17, 0.0)
# I = AminoAcidType("NB_I", "ILE", 0.972973, 6.18, 3, 113.15, 0.0)
# L = AminoAcidType("NB_L", "LEU", 0.972973, 6.18, 4, 113.15, 0.0)
# V = AminoAcidType("NB_V", "VAL", 0.891892, 5.86, 5, 99.13, 0.0)
# W = AminoAcidType("NB_W", "TRP", 0.945946, 6.78, 6, 186.07, 0.0)
# Y = AminoAcidType("NB_Y", "TYR", 0.864865, 6.46, 7, 163.17, 0.0)
# A = AminoAcidType("NB_A", "ALA", 0.72973, 5.04, 8, 71.07, 0.0)
# G = AminoAcidType("NB_G", "GLY", 0.648649, 4.50, 9, 57.05, 0.0)
# T = AminoAcidType("NB_T", "THR", 0.675676, 5.62, 10, 101.10, 0.0)
# S = AminoAcidType("NB_S", "SER", 0.594595, 5.18, 11, 87.07, 0.0)
# N = AminoAcidType("NB_N", "ASN", 0.432432, 5.68, 12, 114.10, 0.0)
# Q = AminoAcidType("NB_Q", "GLN", 0.513514, 6.02, 13, 128.13, 0.0)
# D = AminoAcidType("NB_D", "ASP", 0.378378, 5.58, 14, 115.08, -1.0)
# E = AminoAcidType("NB_E", "GLU", 0.459459, 5.92, 15, 129.10, -1.0)
# H = AminoAcidType("NB_H", "HIS", 0.513514, 6.08, 16, 137.14, 0.5)
# R = AminoAcidType("NB_R", "ARG", 0.0, 6.56, 17, 156.18, 1.0)
# K = AminoAcidType("NB_K", "LYS", 0.513514, 6.36, 18, 128.17, 1.0)
# P = AminoAcidType("NB_P", "PRO", 1.0, 5.56, 19, 97.11, 0.0)

# CN = AminoAcidType("NB_CN", "CYSN", 0.594595, 5.48, 0, 104.14, 1.0)
# MN = AminoAcidType("NB_MN", "METN", 0.837838, 6.18, 1, 132.20, 1.0)
# FN = AminoAcidType("NB_FN", "PHEN", 1.0, 6.36, 2, 148.18, 1.0)
# IN = AminoAcidType("NB_IN", "ILEN", 0.972973, 6.18, 3, 114.16, 1.0)
# LN = AminoAcidType("NB_LN", "LEUN", 0.972973, 6.18, 4, 114.16, 1.0)
# VN = AminoAcidType("NB_VN", "VALN", 0.891892, 5.86, 5, 100.14, 1.0)
# WN = AminoAcidType("NB_WN", "TRPN", 0.945946, 6.78, 6, 187.08, 1.0)
# YN = AminoAcidType("NB_YN", "TYRN", 0.864865, 6.46, 7, 164.18, 1.0)
# AN = AminoAcidType("NB_AN", "ALAN", 0.72973, 5.04, 8, 72.08, 1.0)
# GN = AminoAcidType("NB_GN", "GLYN", 0.648649, 4.50, 9, 58.06, 1.0)
# TN = AminoAcidType("NB_TN", "THRN", 0.675676, 5.62, 10, 102.11, 1.0)
# SN = AminoAcidType("NB_SN", "SERN", 0.594595, 5.18, 11, 88.08, 1.0)
# NN = AminoAcidType("NB_NN", "ASNN", 0.432432, 5.68, 12, 115.11, 1.0)
# QN = AminoAcidType("NB_QN", "GLNN", 0.513514, 6.02, 13, 129.14, 1.0)
# DN = AminoAcidType("NB_DN", "ASPN", 0.378378, 5.58, 14, 116.09, 0.0)
# EN = AminoAcidType("NB_EN", "GLUN", 0.459459, 5.92, 15, 130.11, 0.0)
# HN = AminoAcidType("NB_HN", "HISN", 0.513514, 6.08, 16, 138.15, 1.5)
# RN = AminoAcidType("NB_RN", "ARGN", 0.0, 6.56, 17, 157.19, 2.0)
# KN = AminoAcidType("NB_KN", "LYSN", 0.513514, 6.36, 18, 129.18, 2.0)
# PN = AminoAcidType("NB_PN", "PRON", 1.0, 5.56, 19, 98.12, 1.0)
#
# CC = AminoAcidType("NB_CC", "CYSC", 0.594595, 5.48, 0, 120.14, -1.0)
# MC = AminoAcidType("NB_MC", "METC", 0.837838, 6.18, 1, 148.20, -1.0)
# FC = AminoAcidType("NB_FC", "PHEC", 1.0, 6.36, 2, 164.18, -1.0)
# IC = AminoAcidType("NB_IC", "ILEC", 0.972973, 6.18, 3, 130.16, -1.0)
# LC = AminoAcidType("NB_LC", "LEUC", 0.972973, 6.18, 4, 130.16, -1.0)
# VC = AminoAcidType("NB_VC", "VALC", 0.891892, 5.86, 5, 116.14, -1.0)
# WC = AminoAcidType("NB_WC", "TRPC", 0.945946, 6.78, 6, 203.08, -1.0)
# YC = AminoAcidType("NB_YC", "TYRC", 0.864865, 6.46, 7, 180.18, -1.0)
# AC = AminoAcidType("NB_AC", "ALAC", 0.72973, 5.04, 8, 88.08, -1.0)
# GC = AminoAcidType("NB_GC", "GLYC", 0.648649, 4.50, 9, 74.06, -1.0)
# TC = AminoAcidType("NB_TC", "THRC", 0.675676, 5.62, 10, 118.11, -1.0)
# SC = AminoAcidType("NB_SC", "SERC", 0.594595, 5.18, 11, 104.08, -1.0)
# NC = AminoAcidType("NB_NC", "ASNC", 0.432432, 5.68, 12, 131.11, -1.0)
# QC = AminoAcidType("NB_QC", "GLNC", 0.513514, 6.02, 13, 145.14, -1.0)
# DC = AminoAcidType("NB_DC", "ASPC", 0.378378, 5.58, 14, 132.09, -2.0)
# EC = AminoAcidType("NB_EC", "GLUC", 0.459459, 5.92, 15, 146.11, -2.0)
# HC = AminoAcidType("NB_HC", "HISC", 0.513514, 6.08, 16, 154.15, -0.5)
# RC = AminoAcidType("NB_RC", "ARGC", 0.0, 6.56, 17, 173.19, 0.0)
# KC = AminoAcidType("NB_KC", "LYSC", 0.513514, 6.36, 18, 145.18, 0.0)
# PC = AminoAcidType("NB_PC", "PROC", 1.0, 5.56, 19, 114.12, 0.0)

# aaList = [C, M, F, I, L,
#           V, W, Y, A, G,
#           T, S, N, Q, D,
#           E, H, R, K, P,
#           CN, MN, FN, IN, LN,
#           VN, WN, YN, AN, GN,
#           TN, SN, NN, QN, DN,
#           EN, HN, RN, KN, PN,
#           CC, MC, FC, IC, LC,
#           VC, WC, YC, AC, GC,
#           TC, SC, NC, QC, DC,
#           EC, HC, RC, KC, PC]


aaList = [M, G, K, T, R,
          A, D, E, Y, V,
          L, Q, W, F, S,
          H, N, P, C, I]


def name_to_aa(triname):
    for i in aaList:
        if i.fullname == triname:
            return i
        else:
            pass
    print("No such aa: " + triname + "!\n")
    return 0

# default residue number: res_index=0
# customized residue number: res_index=res_index
def line_to_residue(line, res_index=0):
    t_line_list = [line[0:5], line[5:10], line[10:15], line[15:20], line[20:28], line[28:36],
                   line[36:44]]  # residue number; residue name; atom name; atom number; x; y; z
    #    t_line_list = [line[0:5]]
    #    res_name = "".join([i for i in t_line_list[0] if not i.isdigit()])
    res_name = t_line_list[1].strip()
    index_x = t_line_list[4].index(".")
    index_y = t_line_list[5].index(".")
    index_z = t_line_list[6].index(".")
    
    if res_index == 0: 
        res_index = line[15:20]
    else: 
        pass
    
    t_residue = Residue(int(res_index),
                        name_to_aa(res_name[0:3]),
                        t_line_list[4][0:index_x] + t_line_list[4][index_x + 1] + "." + t_line_list[4][index_x + 2:],
                        t_line_list[5][0:index_y] + t_line_list[5][index_y + 1] + "." + t_line_list[5][index_y + 2:],
                        t_line_list[6][0:index_z] + t_line_list[6][index_z + 1] + "." + t_line_list[6][index_z + 2:])
    return t_residue


def fake_mdp():
    with open("fake.mdp", "w+") as f_w:
        f_w.write(
            "integrator = steep; \ndt = 0.0005; \nnsteps = 500000;\nemtol = 0.1; \nemstep = 0.01; \n\nnstxout = 200;\nnstvout = 200;\nnstfout = 200;\nnstlog = 5000;\nnstenergy = 5000;\n\nenergygrps = system;\ncutoff-scheme = group; \nnstlist = 1; \nns-type = grid; \npbc = xyz;\nrlist = 3.0; \nrcoulomb = 3.0; \nrvdw = 3.0; \ncoulombtype = User; \nvdwtype = User; \nDispCorr = no; \ntable-extension = 0; \ntc-grps = system; \ntau-t = 0.02; \nref-t = 300;\npcoupl = no; \ngen-vel = yes;\ngen-temp = 300; ")


if __name__ == "__main__":
    print(hello())
    chain_length = 0

    with open(gro_fn, "r") as f_r:
        t = f_r.readlines()

    residue_list = []
    # aa_type_list = []
    charge_type_list = []

    index_N = 0
    index_C = 0
    n_chain = 0

    box_line = t[len(t) - 1].split()
    box_x = float(box_line[0]) * 10
    box_y = float(box_line[1]) * 10
    box_z = float(box_line[2]) * 10

    for i in range(2, len(t)):
        if t[i][8] == "N":  # XXXN 5 6 7 8
            index_N = i
        elif t[i][8] == "C":  # XXXC 5 6 7 8
            index_C = i
            break
        else:
            pass

    chain_length = index_C - index_N + 1
    n_chain = int((len(t) - 3) / chain_length)

    for i in range(2, len(t) - 1):
        #    for i in range(2, 2 + chain_length):
        t_res = line_to_residue(t[i], i - 1)
        residue_list.append(t_res)
        # if t_res.aaType not in aa_type_list:
        #     aa_type_list.append(t_res.aaType)
        # else:
        #     pass

    # for i in aa_type_list:
    #     if i.charge != 0.0:
    #         charge_type_list.append(i)
    #     else:
    #         pass

    with open("data." + protein_name, "w+") as f_w:
        f_w.write("\n\n")
        f_w.write("\t" + str(len(residue_list)) + " atoms\n")
        f_w.write("\t" + str(n_chain * (chain_length - 1)) + " bonds\n")
        f_w.write("\t" + str(n_chain * (chain_length - 2)) + " angles\n")
        f_w.write("\t" + str(n_chain * (chain_length - 3)) + " dihedrals\n")
        f_w.write("\t" + str(n_chain * (chain_length - 4)) + " crossterms\n\n")

        f_w.write("\t" + str(len(aaList)*2) + "\tatom types\n")
        f_w.write("\t1\tbond types\n")
        #f_w.write("\t1\tangle types\n")
        f_w.write("\t294\tangle types\n")
        f_w.write("\t1\tdihedral types\n\n")

        f_w.write("\t0.00\t" + "{:10.2f}".format(box_x) + "\txlo xhi\n")
        f_w.write("\t0.00\t" + "{:10.2f}".format(box_y) + "\tylo yhi\n")
        f_w.write("\t0.00\t" + "{:10.2f}".format(box_z) + "\tzlo zhi\n\n")

        f_w.write("\tMasses\n\n")

        for i in range(len(aaList)):
            f_w.write("\t" + str(i + 1) + "\t" + "{:>10.2f}".format(aaList[i].mass) + "\t\t# " + aaList[i].name + "\n")

        for i in range(len(aaList)):
            f_w.write("\t" + str(i + 21) + "\t" + "{:>10.2f}".format(aaList[i].mass) + "\t\t# " + aaList[i].name + "\n")

        f_w.write("\n")
        f_w.write("Atoms\n\n")

        for i in residue_list:
            f_w.write("\t" + "{:>10}".format(i.residue_number) + "\t"
                    + "{:>10}".format(i.residue_number) + "\t"
                    + "{:>7}".format(aaList.index(i.aaType) + 1) + "\t"
                      + "{:>7}".format(i.aaType.charge) + "\t"
                      + i.x + "\t"
                      + i.y + "\t"
                      + i.z + "\t\t"
                      + "# " + i.aaType.name + "\n")

        f_w.write("\n")

        # below for Bond section

        f_w.write("Bond Coeffs\n\n")
        # f_w.write("\t1\t100\t3.82\n\n")
        f_w.write("\t1\t9.600\t3.81\n\n")
        f_w.write("Bonds\n\n")

        #        for i in range(1, len(residue_list)):
        #            f_w.write("\t" + str(i) + "\t\t1\t" + str(i) + "\t\t" + str(i+1) + "\n")
        index = 1
        for i in range(n_chain):
            for j in range(1, chain_length):
                f_w.write("\t" + "{:>10}".format(index) + "\t1\t" + "{:>10}".format(i * chain_length + j) + "\t" + "{:>7}".format(i * chain_length + j + 1) + "\n")
                index = index + 1
        f_w.write("\n")

        # below for Angle section
        #f_w.write("Angles\n\n")
        #index = 1
        #for i in range(n_chain):
        #    for j in range(1, chain_length - 1):
        #        f_w.write("\t" + str(index) + "\t\t1\t" +
        #                  str(i * chain_length + j) + "\t\t" + str(i * chain_length + j + 1) + "\t\t" +
        #                  str(i * chain_length + j + 2) + "\n")
        #        index = index + 1
        #f_w.write("\n")
        ## below for residue-specific section
        f_w.write("Angles\n\n")
        index = 1
        for i in range(n_chain): 
            for j in range(1, chain_length - 1): 
                angle_type = (residue_list[j-1].aaType.angular_pool1 - 1) * 14 * 3 + (residue_list[j].aaType.angular_pool2 - 1) * 3 + residue_list[j+1].aaType.angular_pool3;
                f_w.write("\t" + "{:>10}".format(index) + "\t" + "{:>10}".format(angle_type) + "\t" + 
                        "{:>10}".format(i * chain_length + j) + "\t" + "{:>10}".format(i * chain_length + j + 1) + "\t" +
                        "{:>10}".format(i * chain_length + j + 2) + "\n")
                index = index + 1
        f_w.write("\n")



        # below for Dihedral section
        f_w.write("Dihedrals\n\n")
        index = 1
        for i in range(n_chain):
            for j in range(1, chain_length - 2):
                f_w.write("\t" + "{:>10}".format(index) + "\t1\t" +
                        "{:>10}".format(i * chain_length + j) + "\t" + "{:>10}".format(i * chain_length + j + 1) + "\t" +
                        "{:>10}".format(i * chain_length + j + 2) + "\t" + "{:>10}".format(i * chain_length + j + 3) + "\n")
                index = index + 1
        f_w.write("\n")

        f_w.write("CMAP\n\n")

        #        for i in range(1, len(residue_list) - 3):
        #            f_w.write("\t" + str(i) + "\t\t1\t" + str(i) + "\t\t" + str(i+1) + "\t\t" + str(i+2) + "\t\t" + str(i+3) + "\t\t" + str(i+4) + "\n")
        index = 1
        for i in range(n_chain):
            for j in range(1, chain_length - 3):
                f_w.write("\t" + "{:>10}".format(index) + "\t1\t" + "{:>10}".format(i * chain_length + j) + "\t" + "{:>10}".format(i * chain_length + j + 1) + "\t" + "{:>10}".format(i * chain_length + j + 2) + "\t" + "{:>10}".format(i * chain_length + j + 3) + "\t" + "{:>10}".format(i * chain_length + j + 4) + "\n")
                index = index + 1

    # with open(pair_coeff_fn, "w+") as f_w2:
    #     for i in range(len(aa_type_list)):
    #         for j in range(i, len(aa_type_list)):
    #             f_w2.write("pair_coeff\t" + str(i + 1) + "\t" + str(j + 1) + "\ttable\t../../table/HPS_tables/HPS_" + aa_type_list[i].name + ".table\tHPS_" + aa_type_list[i].name + "_" +aa_type_list[j].name + "\t30.0\n")
    #     for i in range(len(charge_type_list)):
    #         for j in range(i, len(charge_type_list)):
    #             f_w2.write("pair_coeff\t" + str(aa_type_list.index(charge_type_list[i]) + 1) + "\t" + str(aa_type_list.index(charge_type_list[j]) + 1) + "\tcoul/debye\t\t30.0\n")

    fake_mdp()
