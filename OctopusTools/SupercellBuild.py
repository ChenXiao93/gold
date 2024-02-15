import supercell_core as sc
import numpy as np
import os

def str2float(inputlist):
    # Convert list[str] to np.array [float]
    outputlist = []
    for i in range(len(inputlist)):
        outputlist.append(float(inputlist[i]))
    return (np.array(outputlist))

def SupercellInput(inp1, atom1, inp2, atom2, theta, max, vaccum_height, layer_distance):
    # substract layer a;  top layer b
    layer_a = sc.read_POSCAR(inp1, atomic_species=atom1 )
    layer_b = sc.read_POSCAR(inp2, atomic_species=atom2)
    # bilayer
    bilayer = (sc.heterostructure().set_substrate(layer_a).add_layer(layer_b))
    # optimal algorithm
    res = bilayer.opt(max_el=max, thetas=[theta * sc.DEGREE], algorithm='fast', log=True)
    #######  Save VASP1 simply add two layers with Direct coordinates #####################
    res.superlattice().save_POSCAR("POSCAR")
    if os.path.exists("VASP1"):
        os.remove("VASP1")
    os.rename("POSCAR", "VASP1")
    ######   Output Part
    print(" Max Strain is %.2f [percent]" % (res.max_strain() * 100 ))  # max strain  [percent]
    print("Atom number is %s" % res.atom_count())
    ######### Read lattice vectors and the length from VASP1 #####
    with open("VASP1", "r") as lattice:
        vectors = lattice.readlines()
        vector_a = str2float(vectors[2].split())
        vector_b = str2float(vectors[3].split())
        vector_c = str2float(vectors[4].split())
        len_a = np.linalg.norm(vector_a)
        len_b = np.linalg.norm(vector_b)
        len_c = np.linalg.norm(vector_c)

        atomSpeices = vectors[5].split()  # ['Mo', 'S', 'Se']
        atomnumbers = str2float(vectors[6].split())  # [2. 2. 2.]

        ###############
        all_atom = []  # ['Mo', 'Mo', 'S', 'S', 'Se', 'Se']
        for i in range(len(atomSpeices)):
            j = 0
            while j < atomnumbers[i]:
                all_atom.append(atomSpeices[i])
                j += 1
        # print(vector_a, vector_b, vector_c, len_a, len_b, len_c, atomnumbers, atomSpeices, all_atom)
    ######### Read Atom coordinates from VASP1 #####
    f = np.loadtxt("VASP1", skiprows = 8)
    x = f[:,0]  # x reduced coordinates of all  atoms
    y = f[:,1]  # y reduced coordinates of all  atoms
    z = f[:,2]  # z reduced coordinates of all  atoms

    Atom_z_position = np.array(sorted(set(z), key = list(z).index)) * len_c # delete the repeated number
    Position_W1 = Atom_z_position[0]
    Position_W2 = Atom_z_position[1]
    Position_S1 = Atom_z_position[2]
    Position_S2 = Atom_z_position[3]
    Position_Se1 = Atom_z_position[4]
    Position_Se2 = Atom_z_position[5]

    Distance_W_S = np.abs(Position_W1 - Position_S1)
    Distance_W_Se = np.abs(Position_W2 - Position_Se2)

    Number_W1 = int(atomnumbers[1] / 2)  # W of WS2
    Number_W2 = int(atomnumbers[2] / 2)  # W of WSe2

    height_S1 = vaccum_height
    height_W1 = height_S1 + Distance_W_S
    height_S2 = height_W1 + Distance_W_S
    height_Se1 = height_S2 + layer_distance
    height_W2 = height_Se1  + Distance_W_Se
    height_Se2 = height_W2 + Distance_W_Se
    Total_height = height_Se2 + vaccum_height

    Rewrite_Z_Coordinate = np.concatenate([np.ones(Number_W1) * height_W1, np.ones(Number_W2) * height_W2,
                                           np.ones(Number_W1) * height_S1, np.ones(Number_W1) * height_S2,
                                           np.ones(Number_W2) * height_Se1, np.ones(Number_W2) * height_Se2])
    #######  Save POSCAR simply add two layers with Direct coordinates #####################
    with open("POSCAR", "w") as newf:
        for i in [0,1,2,3]:
            newf.write(vectors[i])
        newf.write("0 0 %.3f \n" % Total_height)
        newf.write(vectors[5])  #  Mo S Se
        newf.write(vectors[6])  #  2 2 2
        newf.write("Cartesian\n")

        for j in range( int( sum(atomnumbers) ) ):
            Cartesian_coordinate = (x[j]) * vector_a + (y[j]) * vector_b  + (z[j]) * vector_c
            newf.write("%.6f      %.6f      %.6f\n" % (Cartesian_coordinate[0],Cartesian_coordinate[1], Rewrite_Z_Coordinate[j]))

    os.remove("VASP1")
    # ###### Create Octopus input file TMD1.xyz  original point(0,0,0) at corner  ##############
    with open("TMD1.xyz", "w") as tmd1:
        tmd1.write("%d \n" % int( sum(atomnumbers)) )
        tmd1.write("##Bilayer TMD Cartesian Ccoordinates (0,0,0) at corner \n")

        for j in range( int( sum(atomnumbers) ) ):
            Cartesian_coordinate = (x[j]) * vector_a + (y[j]) * vector_b  + (z[j]) * vector_c
            tmd1.write(" %s    %.6f      %.6f      %.6f\n" % (all_atom[j], Cartesian_coordinate[0], \
                                                              Cartesian_coordinate[1], Rewrite_Z_Coordinate[j]))

    ###### Create Octopus input file TMD2.xyz  original point(0,0,0) at center  ##############
    with open("TMD2.xyz", "w") as tmd2:
        tmd2.write("%d \n" % int( sum(atomnumbers)) )
        tmd2.write("##Bilayer TMD Cartesian Ccoordinates (0,0,0) at center \n")
        octopus_origin = (vector_a + vector_b)/2

        for j in range( int( sum(atomnumbers) ) ):
            Cartesian_coordinate = ((x[j]) * vector_a + (y[j]) * vector_b  + (z[j]) * vector_c) \
                                   - [octopus_origin[0], octopus_origin[1],0]
            tmd2.write(" %s    %.6f      %.6f      %.6f\n" % (all_atom[j],  \
            Cartesian_coordinate[0], Cartesian_coordinate[1], Rewrite_Z_Coordinate[j] -  Total_height/2 ))
    return()
#################################################################################################
def test():

    SupercellInput(inp1 = "../MoS2.dat", \
                   atom1 = ["Mo", "S"], \
                   inp2 = "../MoSe2.dat", \
                   atom2 = ["Mo", "Se"], \
                   theta = 0, \
                   max = 1,\
                   vaccum_height = 30 ,   \
                   layer_distance = 20)
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
