#!/bin/python3
#transfer the gauge field in 'laser' file to Electric field
#unit is eV/nm
Int = open("./e.txt","w")  
laser=open("./laser", "r")
Amplitude=laser.readlines()

Int.write("time[fs]"+"         "+"electricfield[V/nm]"+"\n")
Int.write("0"+"         "+"0"+"\n")
dt=0.016   #  1 [hbar/Ha] = 0.02419 [fs]
c=137.131  # [a.u.]
for i in range(7,130006):
    t1=float(Amplitude[i].split()[1])*0.02419
    a1=float(Amplitude[i].split()[2])
    a0=float(Amplitude[i-1].split()[2])   # 1[Ha]=27.2114[eV], 1 [Bohr]=0.529 [Angstrom]
    e=((-1/c)*(a1-a0)/dt)*51.44*10   #Ha/Bohr eV/Angstrom*e ,  e is a.u. = 1 
    Int.write(str(t1)+"         "+str(e)+"\n")

laser.close()
Int.close()
