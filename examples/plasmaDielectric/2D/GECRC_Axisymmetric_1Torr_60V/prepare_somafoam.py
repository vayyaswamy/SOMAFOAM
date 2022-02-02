# Script to convert BOLSIG+ output to SOMAFOAM input file format
# BOLSIG+ output should be converted to MS-DOS format

from math import *

inputfilename = 'Ar_BOLSIG.txt'
mobilitymodel = 'EON'
diffusionmodel = 'EON'

infile = open(inputfilename)



read = True
while (read == True):
	line = infile.readline()
	a = line.split()
	print a
	if (len(a) > 0):
		if (a[0] == 'R#'):
			read = False

# Writing mobility and diffusion data from BOLSIG+ output
mobility = open('constant/mu_electron','w')
diffusion = open('constant/D_electron','w')
mobility.write('('+'\n')
diffusion.write('('+'\n')
read = True
while (read == True):
	line = infile.readline()
	a = line.split()
	print a 
	
	if (len(a) > 0):
		if (mobilitymodel == 'EON'):
			mobility.write(a[1]+' '+a[3]+'\n')
		elif (mobilitymodel == 'Te'):
			mobility.write(a[2]+' '+a[3]+'\n')
		if (diffusionmodel == 'EON'):
			diffusion.write(a[1]+' '+a[4]+'\n')
		elif (diffusionmodel == 'Te'):
			mobility.write(a[2]+' '+a[4]+'\n')
	else:
		read = False

mobility.write(')'+'\n')
diffusion.write(')'+'\n')

mobility.close()
diffusion.close()

# Writing reaction rate data from BOLSIG+ output
read = True

while (read == True):
	line = infile.readline()
	a = line.split()
	print a
	if (len(a) > 0):
		if (a[0] == 'R#'):
			read = False
			Nr = len(a)-5

for i in range(1,Nr+1):
	outfilename = 'constant/reaction_'+str(i)
	outfile = open(outfilename,'w')
	outfile.write('('+'\n')
	outfile.close()


read = True
while (read == True):
	line = infile.readline()
	a = line.split()
	print a 

	if (len(a) > 0):
		Nr = len(a) - 3
		Te = float(a[2])/1.5
		
		for i in range(1,Nr+1):
			outfilename = 'constant/reaction_'+str(i)
			outfile = open(outfilename,'a')
			if (float(a[i+2]) < 1e-50):
				k = 1e-50
			else:
				k = float(a[i+2])

			outfile.write('('+str(Te*1.602e-19/1.38e-23)+' '+str(log(k*6.023e26))+')\n')
			outfile.close()
		
	else:
		read = False


for i in range(1,Nr+1):
	outfilename = 'constant/reaction_'+str(i)
	outfile = open(outfilename,'a')
	outfile.write(')'+'\n')
	outfile.close()



	