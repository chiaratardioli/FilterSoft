#!/usr/bin/python

# Insert second line in 'results_kepl.dat'

file=open("results_kepl.dat","r")
lines = file.readlines()
file.close()

fout=open("tmp_results","w")

fout.write(lines[0])

# linea2=lines[0].split() 
#
# fout.write(' 0.400000000000E-02'+' '+linea2[1]+' '+linea2[2]+'
# '+linea2[3]+' '+linea2[4]+' '+linea2[5]+' '+linea2[6]+'
# '+linea2[7]+' '+linea2[8]+' '+linea2[9]+'\n')

linea2=lines[0].split()
linea3=lines[1].split()
linea4=lines[2].split()
var=float(linea4[0])-float(linea3[0])
# print linea2[0],linea3[0],linea4[0]
# print '{:.12E}'.format(var)

step = "%14.12E" % var
# print step

#print '     '+ step +' '+linea2[1]+' '+linea2[2]+' '+linea2[3]+'
#'+linea2[4]+' '+linea2[5]+' '+linea2[6]+' '+linea2[7]+' '+linea2[8]+'
#'+linea2[9]+'\n'


fout.write('     '+step+'      '+linea2[1]+'      '+linea2[2]+'      '+linea2[3]+'      '+linea2[4]+'      '+linea2[5]+'      '+linea2[6]+'      '+linea2[7]+'      '+linea2[8]+'      '+linea2[9]+'\n')

del lines[0]

fout.writelines(lines)

fout.close()
