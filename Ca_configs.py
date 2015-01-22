import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from copy import *
plt.close("all")

data = open("./gromos_out.txt").readlines()
out = open("./config.xyz", "w")

NCa = 1
NCl = 2
NWa = 116
N = NCa+NCl+3*NWa # Number of atoms    
r_cut = 0.3 # Radial cutoff  
r_ncut = 0.33 # Radial bond-cutoff
steps = 50
plot = True
save = True

x = np.zeros( (N,3) )
lc = 0
time = 0

for step in range(steps):
    
    atom = 0
    boxlen = 0

    print " +--------------------------------------+ \n\n Reading in atoms."
    for i in range(lc,lc+1000):
        
        line = data[i]
        
        if( "time" in line ):
            command = line[39:57]
            exec command
            print " "+command
            continue
        elif( "box[    0]" in line ):
            boxlen= float(line[19:30])
            continue
        elif( "ox[" in line ):
            continue
       
        if( "x[" in line ):
            command = "x[%d] = %s" % (atom,line[17:55])
            exec command
            atom += 1
            
        if( atom == N ):
            print " All atoms read in."
            lc = i+1
            break
        
    x = (x+boxlen) % boxlen       
    
    xCa = x[0]
    olist = []; hlist = [] # lists of "neighbors"
    otherlist = []
    
    print " \n Ca position: ", xCa
    print " Finding neighbors within |r| < ", r_cut, " nm."
    for i in range(1,atom):         
        if( i%3 == 0 ):
        
            r = x[i] - xCa
            for j in [0,1,2]:
                dist = [abs(r[j]-boxlen), abs(r[j]), abs(r[j]+boxlen)]
                r[j] = min(dist)
            
            if( np.linalg.norm(r) < r_cut ):
                    
                    olist.append(deepcopy(x[i]))              
                    
                    ohlist = []; rlist = []
#                    hlist.append(deepcopy(x[i+1]))
#                    hlist.append(deepcopy(x[i+2]))
                    for k in range(1,atom):          
                        if( k%3 != 0 ):
                            r = x[i] - x[k]
                            for j in [0,1,2]:
                                dist = [abs(r[j]-boxlen), abs(r[j]), abs(r[j]+boxlen)]
                                r[j] = min(dist)
                            rlist.append(np.linalg.norm(r))
                            ohlist.append(deepcopy(x[k]))
                    for k in range(2):
                        index = rlist.index( min(rlist) )
                        hlist.append(deepcopy(ohlist.pop(index)))
                        dummy = rlist.pop(index)
    
    no = len(olist); nh = len(hlist)
    if( 2*no != nh ):
        print "WARNING INCORRECT WATER MOLECULES!"
        print "%d Oxygen and %d Hydrogen" % (no,nh)
    else:
        print " %d found.\n" % no
    
    # Center the Ca
    for i in range(no):
        olist[i] = olist[i]-xCa
        for j in [0,1,2]:
            olist[i][j] = (olist[i][j]+boxlen/2.) % boxlen
    for i in range(nh):
        hlist[i] = hlist[i] -xCa
        for j in [0,1,2]:
            hlist[i][j] = (hlist[i][j]+boxlen/2.) % boxlen
    xCa = np.array([boxlen/2., boxlen/2., boxlen/2.])
    
    if(save):
        
#        xCa = np.array([0.,0.,0.])
#        for i in range(no):
#            r = olist.pop(i)
#            for j in [0,1,2]:
#                dist = [abs(r[j]-boxlen), r[j], (r[j]+boxlen)%boxlen]
#                r[j] = min(dist)
#            olist.insert(i,r)
#        for i in range(nh):
#            r = hlist.pop(i)
#            for j in [0,1,2]:
#                dist = [abs(r[j]-boxlen), r[j], (r[j]+boxlen)%boxlen]
#                r[j] = min(dist)
#            hlist.insert(i,r)
        
        out.write("%d \n" % (1+3*no) )
        out.write("time = %f ps, neighbors %d \n" %(time,no) )
        out.write("F %f %f %f \n" % (10*xCa[0],10*xCa[1],10*xCa[2]) )  ### CHANGE ATOM HERE
        
        for i in range(no):
            out.write("O %f %f %f \n" % (10*olist[i][0],10*olist[i][1],10*olist[i][2]) )
        for i in range(nh):
            out.write("H %f %f %f \n" % (10*hlist[i][0],10*hlist[i][1],10*hlist[i][2]) )
    
    
    olist = np.array(olist); hlist = np.array(hlist)
    if(plot):
        
        fig = plt.figure("Scatter 3d", figsize=(14,12) )
        ax = fig.add_subplot(111, projection='3d')
        
#        ax.scatter(x[:,0],x[:,1],x[:,2], c="k", s=10.)
        
        ax.scatter(xCa[0],xCa[1],xCa[2], c="g", s=150.)
        if( no > 0 ):    
            ax.scatter(olist[:,0],olist[:,1],olist[:,2], c="r", s=150.)
            ax.scatter(hlist[:,0],hlist[:,1],hlist[:,2], c="b", s=150.)
            
        ax.set_xlim([0,boxlen]); ax.set_ylim([0,boxlen]); ax.set_zlim([0,boxlen])
        ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
        ax.set_title("time = %.2f ps  #O = %d\n  r = %.3f nm  #H = %d" %(time,no,r_cut,nh) )
        
        plt.pause(0.2)
        plt.savefig("movie/%dps.png" %time )
        
        
out.close()
      

    
        