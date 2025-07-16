#/bin/python

#This script reads in a textfile containing the cell parameter in 3x3
#a11 a12 a13
#a21 a22 a23
#a31 a32 a33
#Then build a kmesh based on user specified input nk1 nk2 nk3 offset1 offset2 offset3.
#The mesh in crystal coordinate is the same way how OCEAN would do.
#Next it converts the mesh from crystal to cartesian based on the cell parameters and print out
#The kpoint in cryst (which should be the same as OCEAN/QE) and its corresponding magnitude

import os,re,sys
import numpy as np

def create_mesh(NK1,NK2,NK3,offset1,offset2,offset3):
    NK1=int(NK1)
    NK2=int(NK2)
    NK3=int(NK3)
    #print(NK1,NK2,NK3)
    ntot=(NK1*NK2*NK3+7)
    ntot=int(ntot)
   #print(ntot)
    k1=float(offset1)
    k2=float(offset2)
    k3=float(offset3)
    q=np.zeros([ntot,3])
    for i in range(int(NK1)):
            for j in range(int(NK2)):
                for k in range(int(NK3)):
                 iglobal = k + j*int(NK3)+ i*int(NK2)*int(NK3)
                 q[iglobal,0] = float(i)/int(NK1)+k1/int(NK1)
                 q[iglobal,1] = float(j)/int(NK2)+k2/int(NK2)
                 q[iglobal,2] = float(k)/int(NK3)+k3/int(NK3)
                 #print(q[iglobal,:])
    q[iglobal+1,:]=np.array([1.0,0.0,0.0])
    q[iglobal+2,:]=np.array([0.0,1.0,0.0])
    q[iglobal+3,:]=np.array([0.0,0.0,1.0])
    q[iglobal+4,:]=np.array([1.0,1.0,0.0])
    q[iglobal+5,:]=np.array([1.0,0.0,1.0])
    q[iglobal+6,:]=np.array([0.0,1.0,1.0])
    q[iglobal+7,:]=np.array([1.0,1.0,1.0])
    #for i in range(ntot):
        #print(q[i,:])
    return q

def at2bg(File):
    at=np.loadtxt(File)
    print("Read in cell parameters")
    print(at)
    #print(at[:,0])
    print("")
    a1=at[:,0]
    a2=at[:,1]
    a3=at[:,2]
    alat=np.sqrt(np.sum(np.square(a1)))
    #print(a1)
    print("alat=",alat)
    print("")
    a1=a1/alat
    a2=a2/alat
    a3=a3/alat
    print("Crystal axes")
    print(a1)
    print(a2)
    print(a3)
    V=np.dot(np.cross(a1,a2),a3)
    print(V)
    b1=np.cross(a2,a3)/V
    b2=np.cross(a3,a1)/V
    b3=np.cross(a1,a2)/V
    #print(b1)
    #print(b2)
    #print(b3)
    #print(at.shape)
    bg=np.zeros([3,3])
    bg[:,0]=b1
    bg[:,1]=b2
    bg[:,2]=b3
    print("")
    print("Reciprocal axes")
    print(bg)
    return bg

def computeQcart(bg,kmesh):
    ntot=len(kmesh)
    qcart=np.zeros([ntot,7])
    for i in range(ntot):
        qcart[i,0:3]=kmesh[i,:]
        qcart[i,3:6]=np.matmul(bg,kmesh[i,:])
        qcart[i,-1]=np.sqrt(np.sum(np.square(qcart[i,3:6])))
        print("%d %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f" %(i+1,qcart[i,0],qcart[i,1],qcart[i,2],qcart[i,3],qcart[i,4],qcart[i,5],qcart[i,6]))
    #print(kmesh[1,:])
    #print(np.matmul(bg,kmesh[1,:]))

                    



if __name__== '__main__':
        if len(sys.argv) != 8:
         print("Incorrect input")
         sys.exit("USAGE: python kcrys2cart.py CELL.txt NK1 NK2 NK3 offset1 offset2 offset3")
        file1= sys.argv[1]
        if  os.path.isfile(file1) == False :
            sys.exit(file1,"not found. Exitig")
        kmesh=create_mesh(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])
        bg=at2bg(file1)
        #print(len(kmesh))
        computeQcart(bg,kmesh)

        #main()    
