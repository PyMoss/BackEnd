import numpy as np
import os,random
import matplotlib.pyplot as plt
import pylab as pylab
import scipy.optimize
#SECTION: Functions
#This function takes the data and sums multiple points together and returns the new matrix
def dataReducer(data, reductionlevel):
    #Checks to see if amount of data is divisible by the reducion level. If not then this method will not work
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
#This function takes the data and performs folding to remove a parabolic baseline.
def dataFolder(data, foldingpoint):
    newfold = []
    for i in range(0,len(data)//2):
        newfold.append(data[i]+data[len(data)-i-1])
    return newfold
#This function coverts the X axis to velocity and returns the new values.
def dataVelocity(data,zeroPoint,rate):
    pos = []
    neg = []
    for i in range(zeroPoint,len(data)):
        pos.append((i-zeroPoint)*rate)
    for i in range(-zeroPoint,0):
        neg.append((i)*rate)
    newVelocitycityMatrix = np.append(neg,pos)
    return newVelocitycityMatrix
#This function is wrapper for plotting.  
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
#This function calculates Doppler effect from velocity as an input and resonance energy
def dopplerEffect (velocity,resonanceEnergy):
    c=2.997924588E8
    eDoppler= (1+velocity)*resonanceEnergy/c
    return eDoppler
#This function calculates effective mossbauer thickness -- not to be implemented for a long time
def effectiveThickness (LMfactor,cross,abundance,density,thickness,molarmass):
    N=6.022E23
    result= LMfactor*N*abundance*density*thickness/molarmass
#The parameters for Lorentzian fits.
def lorentzian(x,hwhm,cent,intense,back=0):
    numerator =  (hwhm**2 )
    denominator = ( x - (cent) )**2 + hwhm**2
    y = intense*(numerator/denominator)+back
    return y
#Residual function for fitting.
def residuals(p,y,x):
    err = y - lorentzian(x,p)
    return err
#Residual function for fitting multiple parameters.
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
#END SECTION: Functions
#Input parameters ideally not hard coded, but for testing and dev they will be.
filename = "H:\Google\Research\Software\wmoss\TPNIOF14forwmoss.dat" #Test Data Path - Go Through GUI
data = []
data = pylab.loadtxt(filename)
redNumber=4
channels=1024
midpoint= 512//redNumber
#End input parameters
data = dataFolder(data,channels)
new = dataReducer(data,redNumber)
newVelocity = dataVelocity(new,midpoint,16/channels*redNumber)
ind_bg_low = (newVelocity > min(newVelocity)) & (newVelocity < -3)
ind_bg_high = (newVelocity > +3) & (newVelocity < max(newVelocity))
x_bg = np.concatenate((newVelocity[ind_bg_low],newVelocity[ind_bg_high]))
y_bg = np.concatenate((new[ind_bg_low],new[ind_bg_high]))
ind_bg_mid=(newVelocity > -8) & (newVelocity < 8)
m, c = np.polyfit(x_bg, y_bg, 1)
background = m*newVelocity + c
y_bg_corr = new - background
#These are the test values for parameters to be fitted.
#It should be able to accept an unlimited number.
#It is set for 12 sets of parameters. It needs to be able to do 1,2,5,6,12
p = [0.34,-0.2,-800,0.34,-0.15,-750,0.34,0.11,-500,0.34,-0.25,-888,0.34,-0.2,-800,0.34,-0.15,-750,0.34,0.11,-500,0.34,-0.25,-888,0.34,-0.2,-800,0.34,-0.15,-750,0.34,0.11,-500,0.34,-0.25,-888]  # [hwhm, peak center, intensity] #
pbest = scipy.optimize.leastsq(multipleResiduals,p,args=(newVelocity[ind_bg_mid],y_bg_corr[ind_bg_mid]),full_output=1)
fitsum=np.zeros(len(newVelocity))
for i in range(0,len(pbest[0][:]),3):
    fit = lorentzian(newVelocity,pbest[0][i],pbest[0][i+1],pbest[0][i+2],background)
    fitsum= np.add(fitsum,lorentzian(newVelocity,pbest[0][i],pbest[0][i+1],pbest[0][i+2]))
    pylab.plot(newVelocity,fit,'r-',lw=2, label=i)
pylab.plot(newVelocity,new,'b-')
fitsum=np.add(fitsum,background)
pylab.plot(newVelocity,fitsum,'g-',lw=5)
pylab.legend()
pylab.show()

   