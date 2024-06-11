import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Haring&Rix empirical relation data points:
x = [9.48904602563293, 12.649993686470106, 9.880974177662731, 10.38739819433045, 11.05372813940274, 11.542205947345161, 11.937117242250142, 12.238367321169266, 12.413851884588674, 10.654649914767347, 10.101916156323]
y = [6.49868832769683, 10.053146433199533, 6.946663927431334, 7.508802427383925, 8.25204652485856, 8.807073548468662, 9.243623376212902, 9.585969847340309, 9.783574069976927, 7.818625746705016, 7.1893075002370495]

if "Bardati" == "Bardati":
    #Jaeden's relation data points
    x = [3637478143.737464, 4508696720.419435, 8399841594.732584, 29931977349.46067, 57255167589.73556, 174107955231.8671, 337259876154.1427, 466562089950.7191]
    y = [10767571.246586034, 13745564.288824985, 28369849.03122879, 125120661.94171678, 265115454.53887838, 966454716.7444683, 2082718682.3809564, 3029719793.5819206]

    x, y = np.log10(x), np.log10(y)

#We fit a line to it (the slope is written textually in the paper, so we fix it):
def line(x,b):
    return 1.12*x+b

res = curve_fit(line, x, y)[0]

def log10_stellarMass_mBH(log10_mBH):
    '''
    Calculate the stellar Mass log (in solar masses) from the black hole mass log (in solar masses)
    '''
    return (log10_mBH-res[0])/1.12

def stellarMass_mBH(mBH):
    '''
    Calculate the stellar Mass (in solar masses) from the black hole mass (in solar masses)
    '''
    if False:
        if mBH == None:
            return None
        if type(mBH) == int or type(mBH) == float:
            return 10**(log10_stellarMass_mBH(np.log10(mBH)))
        else:
            #Propagate the uncertainties
            step1 = np.log10(mBH[0])
            lo = mBH[1]/(mBH[0]*np.log(10))
            hi = mBH[2]/(mBH[0]*np.log(10))

            step2 = (step1-res[0])/1.12 #Let's consider no uncertainty on res[0]
            lo = np.abs(step2)*np.sqrt(()**2+()**2)
        
        return (10**(log10_stellarMass_mBH(np.log10(mBH[0]))), )
    return 10**(log10_stellarMass_mBH(np.log10(mBH)))

if __name__ == "__main__":
    print(res)
    print(log10_stellarMass_mBH(8))
    print(np.log10(stellarMass_mBH(10**8)))
    x_fit = np.linspace(5,13,100)
    y_fit = line(x_fit, res[0])
    plt.plot(x, y, "o")
    plt.plot(x_fit, y_fit, "--")
    plt.xlabel(r"$log(M_\star) \quad [M_\odot]$", fontsize=16)
    plt.ylabel(r"$log(M_\mathrm{BH}) \quad [M_\odot]$", fontsize=16)
    plt.show()

    plt.plot(y, x, "o")
    plt.plot(y_fit, log10_stellarMass_mBH(y_fit), "--")
    plt.ylabel(r"$log(M_\star) \quad [M_\odot]$", fontsize=16)
    plt.xlabel(r"$log(M_\mathrm{BH}) \quad [M_\odot]$", fontsize=16)
    plt.show()
    