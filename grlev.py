import sys

lifet = {}
lev=[]

# Extract the list of levels with lifetimes from the levels file

with open ("fe_ii.lev") as levfile:
    levs = levfile.readlines()
    for i in levs:
        j=i.split()
        if j[5] != "_":
            lev.append(j[0])

for j in lev:
    j=j[-4:]

# Get the list of transitions from each level from the identified linelist
    with open("feall.new") as fn:
        lines = fn.readlines()
        mtlines = [x for x in lines if j in x[69:80]]
# Print out all the lines in the corrected format.
        for i in mtlines:
            lower = i[56:62] + "_" + i[63:68]
            upper = i[69:74] + "_" + i[75:80]
            if i[1] == " "  :
                if float(i[31:40]) > 0.:
                    i = " G" + i[2:123]
                else: 
                    i = " V" + i[2:133]

            print(i[1],i[2:6],i[11:31],i[31:51],lower," ",upper,i[81:123])

