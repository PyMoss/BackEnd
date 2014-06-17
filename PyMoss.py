print('Hello World')
import numpy as np
import os,random
import matplotlib.pyplot as plt
import pylab as pylab
from numpy import linalg as LA
import scipy.optimize

def dataReducer(data, reductionlevel):
    if len(data) % reductionlevel != 0:
        print("This won't work. The reduction level doesn't match the amount of points.")
        return
    p = 0
    q = reductionlevel
    new = []
    for i in range(0,len(data)//reductionlevel):
        new.append(np.sum(data[p:q]))
        p+=reductionlevel
        q+=reductionlevel
    out=np.array(new)
    return out
def dataFolder(data, foldingpoint):
    newfold = []
    for i in range(0,len(data)//2):
        newfold.append(data[i]+data[len(data)-i-1])
    return newfold
def dataVelocity(data,zeropoint,rate):
    pos = []
    neg = []
    for i in range(zeropoint,len(data)):
        pos.append((i-zeropoint)*rate)
    for i in range(-zeropoint,0):
        neg.append((i)*rate)
    newvel = np.append(neg,pos)
    return newvel 
def mossplot(datain,xval,title="Mossbauer Spectrum",xaxis="Velocity (mm/s)",yaxis="Counts"):
    fig = pylab.plot(xval,datain)
    pylab.grid(True)
    pylab.legend()
    pylab.title(title)
    pylab.xlabel(xaxis)
    pylab.ylabel(yaxis)
    pylab.xlim(min(xval),max(xval))
    pylab.show()
    return
def dopplerEffect (velocity,Esource):
    c=2.998E8
    Egamma= (1+velocity)*Esource/c
def effectiveThickness (LMfactor,cross,abundance,density,thickness,molarmass):
    N=6.022E23
    result= LMfactor*N*abundance*density*thickness/molarmass
def lorentzian(x,hwhm,cent,intense,back=0):
    numerator =  (hwhm**2 )
    denominator = ( x - (cent) )**2 + hwhm**2
    y = intense*(numerator/denominator)+back
    return y
def residuals(p,y,x):
    err = y - lorentzian(x,p)
    return err
def multipleResiduals(p,x,yval):
    parin= np.zeros(len(x))
    for i in range(0,len(p),3):
        p0=p[i]
        p1=p[i+1]
        p2=p[i+2]
        if p0>0.5 or p0<0: return 1E8
        if p2>-10: return 1E8
        parin=np.add(parin,lorentzian(x,p0,p1,p2))
    err = yval - parin
    return err

filename = "H:\Google\Research\Software\wmoss\TPNIOF14forwmoss.dat"
data = []
data = pylab.loadtxt(filename)
rednumber=4
channels=1024
midpoint= 512//rednumber
data = dataFolder(data,1024)
new = dataReducer(data,rednumber)
newvelo = dataVelocity(new,midpoint,16/channels*rednumber)
#mossplot(new,newvelo)
ind_bg_low = (newvelo > min(newvelo)) & (newvelo < -3)
ind_bg_high = (newvelo > +3) & (newvelo < max(newvelo))
x_bg = np.concatenate((newvelo[ind_bg_low],newvelo[ind_bg_high]))
y_bg = np.concatenate((new[ind_bg_low],new[ind_bg_high]))
ind_bg_mid=(newvelo > -8) & (newvelo < 8)
m, c = np.polyfit(x_bg, y_bg, 1)
background = m*newvelo + c
y_bg_corr = new - background
p = [0.34,-0.2,-800,0.34,-0.15,-750,0.34,0.11,-500,0.34,-0.25,-888,0.34,-0.2,-800,0.34,-0.15,-750,0.34,0.11,-500,0.34,-0.25,-888,0.34,-0.2,-800,0.34,-0.15,-750,0.34,0.11,-500,0.34,-0.25,-888]  # [hwhm, peak center, intensity] #
pbest = scipy.optimize.leastsq(multipleResiduals,p,args=(newvelo[ind_bg_mid],y_bg_corr[ind_bg_mid]),full_output=1)
fitsum=np.zeros(len(newvelo))
for i in range(0,len(pbest[0][:]),3):
    fit = lorentzian(newvelo,pbest[0][i],pbest[0][i+1],pbest[0][i+2],background)
    fitsum= np.add(fitsum,lorentzian(newvelo,pbest[0][i],pbest[0][i+1],pbest[0][i+2]))
    pylab.plot(newvelo,fit,'r-',lw=2, label=i)
pylab.plot(newvelo,new,'b-')
fitsum=np.add(fitsum,background)
pylab.plot(newvelo,fitsum,'g-',lw=5)
pylab.legend()

pylab.show()
input("Press Enter to continue...")

   