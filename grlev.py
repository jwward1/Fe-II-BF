#
# Program to extract the lines with lifetimes from the identified linelist feall.new and write in correct format
#
lev=[]

levels_filename = "fe_ii.lev"
lines_filename = "feall.new"


# Extract the list of levels with lifetimes from the levels file

with open (levels_filename) as levfile:
    levs = levfile.readlines()
    for i in levs:
        j=i.split()
        if j[5] != "_":
            lev.append((j[0],j[1],j[2]))

for j in lev:
    levid=j[0][-4:]
#    print(levid,j)

# Get the list of transitions from each level from the identified linelist
    with open(lines_filename) as fn:
        lines = fn.readlines()
        mtlines = [x for x in lines if levid in x[69:80]]
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
