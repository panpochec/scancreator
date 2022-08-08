from pathlib import Path
import numpy as np
import pandas as pd



### Specifing keys for searching ###
key1 = "Center     Atomic      Atomic             Coordinates (Angstroms)"
key2 = "symmetry adapted cartesian basis functions of"
key3 = "Symbolic Z-matrix:"
key4 = "Initialization pass."
key5 = "Trust Radius"


def finding_rows(what, where, start, finish):
    """Function allowing scanning of a given matrix's in search for a given parameter.

    Parameters
    ----------
    where : list
        list to scan through
    what : str
        search for a given string
    start: int
        denominates begining of 'where'
    start: int
        denominates end of 'where'
    Returns
    -------
    list
        Returns list of where indexes numbers.
    """
    temp_list = []              ### Temporary list, returned as output
    for x in range(0, len(where)):
        line_part = where[x]
        if str(what) in line_part[int(start): int(finish)]:
            temp_list.append(str(x + 1))
        else:
            pass
    return temp_list


def finding_in_rows(what, where, start, finish, str_start, str_finish):
    """Function allowing scanning of a given matrix's in search for a given string.
    'str_start' and 'str_finish' allow searching only a given part of the line
    (useful when given letter or string is non-specific, like D denominating dihedrals
    but also denominating change of energy at the end of the line)

    Parameters
    ----------
    what : list
        search for a string from a given list (0 of list due to convention)
    where : list
        Scanning of a given matrix
    start: int
        denominates begining of 'where'
    start: int
        denominates end of 'where'
    str_start : int
        allow searching only from a given point in a line of 'where'
    str_finish : int
        allow searching only to a given point in a line of 'where'

    Returns
    -------
    list
        Returns list of string contents.
    """
    temp_list = []              ### Temporary list, returned as output
    for x in range(int(start), int(finish)):
        line_part = where[x]
        if str(what[0]) in line_part[str_start: str_finish]:
            temp_list.append(str(where[x]))
        else:
            pass
    return temp_list

########### START OF MAIN FILE################

### Acquiring  log filelist ###
zrob = 0
dirPath = Path.cwd()
result =  [f.name for f in dirPath.glob('*.log') if f.is_file()]
ile = len(result)

### Iterating through log filelist ###
for y in range(0,ile):
    xyz = ''
    funkc = "wB97xD"
    baza = "6-311+G(2d,2p)"
    argu = result[y]
    argu2 = argu.replace("log", "xyz")
    argu3 = argu.replace(".log", "")
    argu4 = argu3 +"_BridgeOne.xyz"
    argu5 = argu3 +"_BridgeTwo.xyz"
    file = open(result[y], "r")
    lines = file.readlines()
    file.close()
    with_key1 = []
    with_key2 = []


    ### Appending new file content with header ###
    xyz = xyz + "#" + funkc + "/" + baza + " opt=modredundant\n \n"
    xyz = xyz + argu2 + "\n" + "\n" + "0 1 \n"

    ######### Atom + coordinates matrix for xyz creator #######

    for number, line in enumerate(lines, 0):        ### Finding key phrases in file
            if key1 in line:
                with_key1.append(number)
    for number, line in enumerate(lines, 0):
            if key2 in line:
                with_key2.append(number)
    last_key1 = with_key1[-1]
    last_key2 = with_key2[-1]
    start = int(last_key1) + 3
    finish = int(last_key2) -3
    for x in range(start, finish):                  ### Replacing atomic number with atom symbol
        str1 = lines[x]
        strp2 = str1 [30:70]
        strp = str1[17]
        if strp == "6":
            strp = "C"
        elif strp == "1":
            strp = "H"
        elif strp == "7":
            strp = "N"
        elif strp == "8":
            strp = "O"
        elif strp == "5":
            strp = "Br"
            strp2 = str1[31:70]
        else:
            print("ERROR")
        strf= strp + strp2 + "\n"
        xyz= xyz + strf

                            ###### Finding bridge atoms parameters #####
    for number, line in enumerate(lines, 0):        ### Line before atom-type matrix
        if key3 in line:
            a_matrix_start = int(number)+2
    for number, line in enumerate(lines, 0):        ### Line after atom-type matrix
        if key4 in line:
            a_matrix_end = int(number)-4
            p_matrix_start = int(number) + 7        ### Line begining parameter matrix
            break
    atom_matrix = lines[a_matrix_start:a_matrix_end]                       ### Creating atoms matrix

                                 ### Finding oxygens in rows
    rows_for_O = finding_rows("O", atom_matrix, 0, -1)

                                 ### Finding nitrogens in rows
    rows_for_N = finding_rows("N", atom_matrix, 0, -1)

                                 ### Finding hydrogens in rows
    rows_for_H = finding_rows("H", atom_matrix, 0, -1)




    for number, line in enumerate(lines, 0):        ### Finding line after parameter matrix
        if key5 in line:
            p_matrix_end = int(number)-1
            break

    parameter_matrix = lines[p_matrix_start:p_matrix_end]                      ### Parameter Matrix

    ### Dihedral matrix ###

    dihedral_matrix = []        ### Finding dihedral rows  in parameter matrix
    dihedral_matrix = finding_rows("D", parameter_matrix, 0, 20)

    col1 =[]
    col2 =[]
    col3 =[]
    col4 =[]
    col5 =[]

    ### Due to the construction of lines in parameter matrix,
    ### searching for atom numbers and building dataframe is done
    ### using a trick:
    ### We iterate thorugh each line of data matrix searching for a specific
    ### denominator, then appending respective column list and deleting the appended
    ### part with the denominator to make further iterations easier.
    ### In the case of dihedrals :
    ### ! D44   D(16,17,39,40)        179.9269         estimate D2E/DX2
    ### First we delete the part "! D44   D("
    ### Then we take count to the first ",", and append column list with 16
    ### and use the newly created entry '16' and denominator ',' and
    ### delete them from line. The resulting line:
    ### 17,39,40)        179.9269         estimate D2E/DX2
    ### We repeat the process until we acquire all of the parameters.
    ### Same process is used for angles and distances

    for x in range(int(dihedral_matrix[0])-1 ,int(dihedral_matrix[-1])):
        help_string = str(parameter_matrix[x])
        help_string = help_string.replace(str(parameter_matrix[x])[0:11], "")
        col1.append(int(str(help_string)[0: str(help_string).index(",", 0,50)]))
        help_string = help_string.replace(str(col1[x-(int(dihedral_matrix[0])-1)]) + ",", "", 1)
        col2.append(int(str((help_string))[0: str((help_string)).index(",", 0)]))
        help_string = help_string.replace(str(col2[x-(int(dihedral_matrix[0])-1)]) + ",", "", 1)
        col3.append(int(str((help_string))[0: str((help_string)).index(",", 0)]))
        help_string = help_string.replace(str(col3[x-(int(dihedral_matrix[0])-1)]) + ",", "", 1)
        col4.append(int(str((help_string))[0: str((help_string)).index(")", 0)]))
        help_string = help_string.replace(str(col4[x-(int(dihedral_matrix[0])-1)]) + ")", "", 1)
        help_string = help_string.replace(" ", "")
        col5.append(float(str((help_string))[0: str((help_string)).index("e", 0)]))

    df_dihedral = pd.DataFrame({
        'Atom 1': col1,
        'Atom 2': col2,
        'Atom 3': col3,
        'Atom 4': col4,
        'Value': col5
    })

    ### Angle matrix ###

    angle_matrix = []  ### Finding angle rows in parameter matrix
    angle_matrix = finding_rows("A", parameter_matrix, 0, 20)
    print("Rows for dihedral matrix:")
    print(angle_matrix)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    for x in range(int(angle_matrix[0]) - 1, int(angle_matrix[-1])):
        help_string = str(parameter_matrix[x])
        help_string = help_string.replace(str(parameter_matrix[x])[0:11], "")
        col1.append(int(str(help_string)[0: str(help_string).index(",", 0, 50)]))
        help_string = help_string.replace(str(col1[x - (int(angle_matrix[0]) - 1)]) + ",", "", 1)
        col2.append(int(str((help_string))[0: str((help_string)).index(",", 0)]))
        help_string = help_string.replace(str(col2[x - (int(angle_matrix[0]) - 1)]) + ",", "", 1)
        col3.append(int(str((help_string))[0: str((help_string)).index(")", 0)]))
        help_string = help_string.replace(str(col3[x - (int(angle_matrix[0]) - 1)]) + ")", "", 1)
        help_string = help_string.replace(" ", "")
        col4.append(float(str((help_string))[0: str((help_string)).index("e", 0)]))

    df_angle = pd.DataFrame({
        'Atom 1': col1,
        'Atom 2': col2,
        'Atom 3': col3,
        'Value': col4
    })

    ### Distance matrix ###

    distance_matrix = []  ### Finding angle rows in parameter matrix
    distance_matrix = finding_rows("R", parameter_matrix, 0, 20)
    print("Rows for dihedral matrix:")
    print(distance_matrix)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    for x in range(int(distance_matrix[0]) - 1, int(distance_matrix[-1])):
        help_string = str(parameter_matrix[x])
        help_string = help_string.replace(str(parameter_matrix[x])[0:11], "")
        col1.append(int(str(help_string)[0: str(help_string).index(",", 0, 50)]))
        help_string = help_string.replace(str(col1[x]) + ",", "", 1)
        col2.append(int(str((help_string))[0: str((help_string)).index(")", 0)]))
        help_string = help_string.replace(str(col2[x]) + ")", "", 1)
        help_string = help_string.replace(" ", "")
        col3.append(float(str((help_string))[0: str((help_string)).index("e", 0)]))

    df_distance = pd.DataFrame({
        'Atom 1': col1,
        'Atom 2': col2,
        'Value': col3
    })


    ### Assigning bridge atoms###
    bridge_atoms1 = []
    bridge_atoms2 = []
    for x in range(0, len(df_distance.iloc[:, 0])):         ### Oxygens and hydrogens in bridges
        if str(df_distance.iloc[x,0]) in rows_for_O and str(df_distance.iloc[x,1]) in rows_for_H:
            bridge_atoms1.append(df_distance.iloc[x,1])
            bridge_atoms2.append(df_distance.iloc[x,0])
        elif str(df_distance.iloc[x,1]) in rows_for_O and str(df_distance.iloc[x,0]) in rows_for_H:
            bridge_atoms1.append(df_distance.iloc[x,1])
            bridge_atoms2.append(df_distance.iloc[x,0])

    bridge_oxygens = []
    bridge_hydrogens = []

    for x in range(0, len(bridge_atoms1)):   ### Discrimination between oxygens and hydrogens in bridges
        if bridge_atoms1[x] in rows_for_O:
            bridge_oxygens.append(bridge_atoms1[x])
            bridge_hydrogens.append(bridge_atoms2[x])
        else:
            bridge_oxygens.append(bridge_atoms2[x])
            bridge_hydrogens.append(bridge_atoms1[x])

    for row in df_distance.iterrows():      ### Assigining bridge atoms to respective bridges
        if " " + str(bridge_oxygens[0]) +".00" in str(row) and " " + str(bridge_hydrogens[0]) +".00" in str(row):
            bridge_a_oxygen = bridge_oxygens[0]
            bridge_a_hydrogen = bridge_hydrogens[0]
            bridge_b_oxygen = bridge_oxygens[1]
            bridge_b_hydrogen = bridge_hydrogens[1]
        elif " " + str(bridge_oxygens[0]) +".00" in str(row) and " " + str(bridge_hydrogens[1]) +".00" in str(row):
            bridge_a_oxygen = bridge_oxygens[0]
            bridge_a_hydrogen = bridge_hydrogens[1]
            bridge_b_oxygen = bridge_oxygens[1]
            bridge_b_hydrogen = bridge_hydrogens[0]


    #Selecting nitrogene for each bridge

    nitrogens_nitro = []
    nitrogens_amine = []

    for row in df_angle.iterrows():         ### Assigning nitrogens to potential substituents
        print(row)
        for n in range(0, len(rows_for_N)):
            if " " + str(rows_for_N[n])+".00" in str(row):
                for z in range(0, (len(rows_for_H)-1)):
                    if " " + str(rows_for_H[z])+".00" in str(row) and " " + str(rows_for_H[z+1])+".00" in str(row):
                        nitrogens_amine.append(rows_for_N[n])
                for t in range(0, (len(rows_for_O)-1)):
                    if " " + str(rows_for_O[t])+".00" in str(row) and " " + str(rows_for_O[t+1])+".00" in str(row):
                        nitrogens_nitro.append(rows_for_N[n])


    bridge_nitrogens =[]                ### Assigning rest to the bridge nitrogens group
    for x in range(0, len(rows_for_N)):
        if rows_for_N[x] not in nitrogens_nitro and rows_for_N[x] not in nitrogens_amine:
            bridge_nitrogens.append(int(rows_for_N[x]))

    ### Assigining bridge nitrogens to respective bridges using  distances
    ### containing nitrogen and dihedrals containing bridge oxygens

    for row in df_distance.iterrows():
        for t in range(0, len(bridge_nitrogens)):
            if " " + str(bridge_nitrogens[t]) + ".00" in str(row):
                dist_atom_a = df_distance.iloc[(int(str(row)[1: int(str(row).index(","))])),0]
                dist_atom_b = df_distance.iloc[(int(str(row)[1: int(str(row).index(","))])),1]
                for drows in df_dihedral.iterrows():
                    if " " + str(dist_atom_a) + ".0" in str(drows) and " " + str(bridge_a_oxygen) + ".0" in str(
                            drows):
                        bridge_a_nitrogen = dist_atom_b
                    elif " " + str(dist_atom_b) + ".0" in str(drows) and " " + str(bridge_a_oxygen) + ".0" in str(
                                drows):
                        bridge_a_nitrogen = dist_atom_a
                    else:
                        if " " + str(dist_atom_a) + ".0" in str(drows) and " " + str(bridge_b_oxygen) + ".0" in str(
                                drows):
                            bridge_b_nitrogen = dist_atom_b
                        elif " " + str(dist_atom_b) + ".0" in str(drows) and " " + str(bridge_b_oxygen) + ".0" in str(
                                    drows):
                            bridge_b_nitrogen = dist_atom_a

    print("Bridge a is "+ str(bridge_a_nitrogen)+" , "+ str(bridge_a_hydrogen)+" , "+ str(bridge_a_oxygen)+"\n" )
    print("Bridge b is "+ str(bridge_b_nitrogen)+" , "+ str(bridge_b_hydrogen)+" , "+ str(bridge_b_oxygen) )

    bridge_one_csv = ''
    bridge_one_csv = bridge_one_csv + xyz
    bridge_one_csv = bridge_one_csv + "\n" + str(bridge_a_oxygen) +" "+ str(bridge_a_hydrogen)+" "+ str(
        bridge_a_nitrogen) +" F\n" + str(bridge_a_oxygen)+ " " + str(bridge_a_hydrogen) + " S 20 0.05 \n\n"
    bridge_two_csv = ''
    bridge_two_csv = bridge_two_csv + xyz
    bridge_two_csv = bridge_two_csv + "\n" + str(bridge_b_oxygen) + " " + str(
        bridge_b_hydrogen) + " " + str(bridge_b_nitrogen) + " F\n" + str(bridge_b_oxygen) + " " + str(
        bridge_b_hydrogen) + " S 20 0.05 \n\n"

    LOGFILE = open(argu4, "a")
    LOGFILE.write(bridge_one_csv)
    LOGFILE.close()

    LOGFILE2 = open(argu5, "a")
    LOGFILE2.write(bridge_two_csv)
    LOGFILE2.close()
    print(bridge_atoms1)
