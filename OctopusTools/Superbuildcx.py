import numpy as np
import os
import OctopusTools.Basic as cx

def str2float(inputlist):
    # Convert list[str] to np.array [float]
    outputlist = []
    for i in range(len(inputlist)):
        outputlist.append(float(inputlist[i]))
    return (np.array(outputlist))

def SupercellInput(inp1, atom1, inp2, atom2, theta, max, vaccum_height, layer_distance):

    ###  Produce the VASP2 file!!!! ######
    with open(inp1, "r") as lattice1:
        vectors1 = lattice1.readlines()
        vector_a1 = str2float(vectors1[2].split())
        vector_b1 = str2float(vectors1[3].split())
        vector_c1 = str2float(vectors1[4].split())
        Distance_W_S = (str2float(vectors1[9].split())[2] - str2float(vectors1[8].split())[2]) \
                    * np.linalg.norm(vector_c1)
    with open(inp2, "r") as lattice2:
        vectors2 = lattice2.readlines()
        vector_a2 = str2float(vectors2[2].split())
        vector_b2 = str2float(vectors2[3].split())
        vector_c2 = str2float(vectors2[4].split())
        Distance_W_Se = (str2float(vectors2[9].split())[2] - str2float(vectors2[8].split())[2]) \
                    * np.linalg.norm(vector_c2)

    len_a1 = np.linalg.norm(vector_a1)
    len_a2 = np.linalg.norm(vector_a2)
    len_a3 = (len_a1 + len_a2) / 2
    ####  new lattice vector of new supercell ##########
    vector_a3 =  len_a3 * np.array([1,0,0])
    vector_b3 =  len_a3 * np.array([-0.5, np.sqrt(3)/2,0])
    vector_c3 = vector_c1 + vector_c2
    len_a3 = np.linalg.norm(vector_a3)
    len_b3 = np.linalg.norm(vector_b3)
    len_c3 = np.linalg.norm(vector_c3)

    #### Direct coordinates ####
    Mo1_position  = [0, 0, 0]
    Mo2_position = [0, 0, 0]
    S1_position = [0.333333, 0.666667, 0 ]
    S2_position = [0.333333, 0.666667, 0]
    Se1_position = [0.333333, 0.666667, 0]
    Se2_position = [0.333333, 0.666667, 0]
    ### Real Coordinates ####
    Mo1_position_C = Mo1_position[0] * vector_a3 + Mo1_position[1] * vector_b3 + Mo1_position[2] * vector_c3
    Mo2_position_C = Mo2_position[0] * vector_a3 + Mo2_position[1] * vector_b3 + Mo2_position[2] * vector_c3
    S1_position_C = S1_position[0] * vector_a3 + S1_position[1] * vector_b3 + S1_position[2] * vector_c3
    S2_position_C = S2_position[0] * vector_a3 + S2_position[1] * vector_b3 + S2_position[2] * vector_c3
    Se1_position_C = Se1_position[0] * vector_a3 + Se1_position[1] * vector_b3 + Se1_position[2] * vector_c3
    Se2_position_C = Se2_position[0] * vector_a3 + Se2_position[1] * vector_b3 + Se2_position[2] * vector_c3

    ##################################################################################################
    height_S1 = vaccum_height
    height_W1 = height_S1 + Distance_W_S
    height_S2 = height_W1 + Distance_W_S
    height_Se1 = height_S2 + layer_distance
    height_W2 = height_Se1  + Distance_W_Se
    height_Se2 = height_W2 + Distance_W_Se
    Total_height = height_Se2 + vaccum_height

    Rewrite_Z_Coordinate = [height_W1, height_W2, height_S1,  height_S2, height_Se1,  height_Se2]

    with open("POSCAR", "w") as f:
        f.write("Bilayerd POSCAR set Z = 0 \n")
        f.write("1.0\n")
        f.write("%.6f    %.6f    %.6f \n" % (vector_a3[0],vector_a3[1],vector_a3[2] ))
        f.write("%.6f    %.6f    %.6f \n" % (vector_b3[0], vector_b3[1], vector_b3[2]))
        f.write("%.6f    %.6f    %.6f \n" % (vector_c3[0], vector_c3[1], Total_height))
        f.write("%s %s %s\n" % (atom1[0], atom1[1], atom2[1]) )
        f.write("2 2 2 \nCartesian \n")
        f.write("%.6f    %.6f    %.6f \n" % (Mo1_position_C[0], Mo1_position_C[1], Rewrite_Z_Coordinate[0]))
        f.write("%.6f    %.6f    %.6f \n" % (Mo2_position_C[0], Mo2_position_C[1], Rewrite_Z_Coordinate[1]))
        f.write("%.6f    %.6f    %.6f \n" % (S1_position_C[0], S1_position_C[1], Rewrite_Z_Coordinate[2]))
        f.write("%.6f    %.6f    %.6f \n" % (S2_position_C[0],S2_position_C[1], Rewrite_Z_Coordinate[3]))
        f.write("%.6f    %.6f    %.6f \n" % (Se1_position_C[0], Se1_position_C[1], Rewrite_Z_Coordinate[4]))
        f.write("%.6f    %.6f    %.6f \n" % (Se2_position_C[0], Se2_position_C[1], Rewrite_Z_Coordinate[5]))

    ###### Create Octopus input file TMD.xyz  original point(0,0,0) at center  ##############
    with open("TMD.xyz", "w") as tmd:
        tmd.write("%d \n" % 6 )
        tmd.write("##Bilayer TMD Cartesian Ccoordinates (0,0,0) at center \n")
        octopus_origin = (vector_a3 + vector_b3)/2 + [0,0, Total_height / 2]
        all_atom = [atom1[0], atom1[0], atom1[1], atom1[1], atom2[1], atom2[1]]

        with open("POSCAR", "r") as file:
            pos = file.readlines()
            for j in np.arange(8,14,1):
                pos1 = str2float(pos[j].split())
                Cartesian_coordinate = pos1 -  octopus_origin
                tmd.write(" %s    %.6f      %.6f      %.6f\n" % (all_atom[j-8],  \
                    Cartesian_coordinate[0], Cartesian_coordinate[1], Cartesian_coordinate[2]))
    return()

#################################################################################################
def test():

    SupercellInput(inp1 = "./MoS2", \
                   atom1 = ["Mo", "S"], \
                   inp2 = "./MoSe2", \
                   atom2 = ["Mo", "Se"], \
                   theta = 0, \
                   max = 1,\
                   vaccum_height = 10 ,   \
                   layer_distance = 5)
    ########## Input variables ################################
    # inp1:  VASP input file1 for substract layer
    # atom1:  atom configuration of inp1
    # inp2:  VASP input file2 for top layer
    # atom2: atom configuration of inp2
    # theta: twisted angle of two layer
    # max:  variable max_el of Supercore_cell
    # vaccum_height: angstrom
    # layer_distance : angstrom
    # OutputFormat = "vasp" or "octopus"
    # CoordinateSystem = "xyz" or "reduced"
    ##########  Output files ##########################################
    #  VASP1: simply add two layers with Direct coordinates

    return()

if __name__ == "__main__" :
    test()
