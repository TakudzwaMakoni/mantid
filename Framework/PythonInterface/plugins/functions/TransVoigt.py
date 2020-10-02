from mantid.api import IFunction1D, FunctionFactory
import math
import numpy as np
from numpy import vectorize

#Approximation to Voigt function 
#written by C J Ridley for Transfit v2 in Mantid
#Last updated: 28th July 2020
#
#Define approximation to Voigt function consistent with transfit v1 see Igor TN026 for technical details and references


def VoigtFunction(X,Y):
    Y=abs(Y)
    S=abs(X)+Y
    T= Y - (X * 1j)
        
    #Determine values based on value of S
    #REGION 1
    if (S >= 15):
        W= T*0.5641896/(0.5+T*T)
        
    #REGION 2
    else:
        if (S >= 5.5):
            U = T*T
            W= T*(1.410474+U*0.5641896)/(0.75+U*(3+U))
    #REGION 3
        else:
            if (Y>=(0.195*np.abs(X) - 0.176)):
                W= (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+T*0.5642236))))
                W /= (16.4955+T*(38.82363+T*(39.27121+T*(21.69274+T*(6.699398+T)))))
    #REGION 4
            else:
                U= T*T
                W= T*(36183.31-U*(3321.9905-U*(1540.787-U*(219.0313-U*(35.76683-U*(1.320522-U*0.56419))))))
                W /= (32066.6-U*(24322.84-U*(9022.228-U*(2186.181-U*(364.2191-U*(61.57037-U*(1.841439-U)))))))
                W= (np.exp(np.real(U))*np.cos(np.imag(U)) +0j)- W
    return np.real(W)


class TransVoigt(IFunction1D):
    
    def init(self):
        # Starting parameters as fitted from run PRL111643
        self.declareParameter("Position", 1096.3)
        self.declareParameter("Lorentzian FWHM",45.8)
        self.declareParameter("Gaussian FWHM",25.227)
        self.declareParameter("Amplitude",2.08)
        #Legacy parameters from Transfit v1:
        self.declareParameter("b0",25)
        self.declareParameter("b1",0.015)
        self.declareParameter("b2",0)
        #self.declareParameter("b3",0)
        #self.declareParameter("b4",0)
        #self.declareParameter("offset",1)
        #self.declareParameter("shift",0)
       
    def function1D(self, xvals):
        #Vectorise the function to allow if/else statements to work quickly
        vVoigtFunction = vectorize(VoigtFunction)
        #Declare fitting parameters
        pos = self.getParameterValue("Position")
        amp = self.getParameterValue("Amplitude")
        Lor = self.getParameterValue("Lorentzian FWHM")
        Gauss = self.getParameterValue("Gaussian FWHM")
        #Legacy parameters from Transfit v1:
        b0 = self.getParameterValue("b0")
        b1 = self.getParameterValue("b1")
        b2 = self.getParameterValue("b2")
        #b3 = self.getParameterValue("b3")
        #b4 = self.getParameterValue("b4")
        #offset = self.getParameterValue("offset")
        #shift = self.getParameterValue("shift")
        #
        #Legacy background function included from Transfit v1
        #Define background function
        bkgrnd = b0 + b1*xvals + b2*xvals*xvals
        #+ b3*xvals*xvals*xvals + b4*xvals*xvals*xvals*xvals
        #
        #Correct using Beer's law to fit measured absorption, not Xsection
        #np.sqrt(np.log(2)) replaced with 1 as legacy
        width = 1 / Gauss
        shape =  Lor * width
        out=np.exp(-1*np.abs(amp*vVoigtFunction(width*(xvals-pos), shape)))
        out=bkgrnd*out
        #out=bkgrnd*((1-offset)+offset*out)+shift
        return out
    
FunctionFactory.subscribe(TransVoigt)