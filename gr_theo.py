#
# Program to find the lines in the theoretical calculations that have lifetimes and 
# print them out with a key that matches the levels file and linelist file
#

lev=[]

levels_filename = "fe_ii.lev"
lines_filename = "FeII_waveno_tot.E1"

# Extract the list of levels with lifetimes from the levels file

with open (levels_filename) as levfile:
    levs = levfile.readlines()
    for i in levs:
        j=i.split()
        if j[5] != "_":                     # i.e. level has a lifetime
            lev.append((j[0]+"*",j[1],j[2]))

for i in lev:
#    print(i[0],i[1],i[2])
    enlev = float(i[2])
    jval = i[1]

    with open(lines_filename) as fl:
        lines = fl.readlines()
#
# Match the levels based on the J-value and the energy. Difference should be less than 1 cm-1.
#
        mtlines = [x for x in lines if jval in x[74:79] and -1. < enlev - float(x[81:89 ]) < 1. ] 

    for j in mtlines:
        if jval in j[74:79] and -1. < enlev - float(j[81:89 ]) < 1. :
            print("%100s %s" % (j.strip(), i[0]))