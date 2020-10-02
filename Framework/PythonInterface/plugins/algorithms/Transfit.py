# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
# -----------------------------------------------------------------------
# Transfit v2 - translated from Fortran77 to Python/Mantid
# Original Author: Christopher J Ridley
# Last updated: 2nd October 2020
#
# Program fits a Voigt function to PEARL transmission data, using the doppler broadening
# of the gaussian component to determine the sample temperature
#
#
# References:
# 'Temperature measurement in a Paris-Edinburgh cell by neutron resonance spectroscopy' - Journal Of Applied Physics 98,
# 064905 (2005)
# 'Remote determination of sample temperature by neutron resonance spectroscopy' - Nuclear Instruments and Methods in
# Physics Research A 547 (2005) 601-615
# -----------------------------------------------------------------------
from mantid.kernel import Direction
from mantid.api import PythonAlgorithm, MultipleFileProperty, AlgorithmFactory
from mantid.simpleapi import LoadRaw, DeleteWorkspace
import numpy as np
import re


class Transfit(PythonAlgorithm):

    def version(self):
        return 1

    def name(self):
        return 'Transfit'

    def category(self):
        return 'PEARL tools'

    def summary(self):
        return TODO

    def PyInit(self):
        # TODO File directory can probably be removed, replace with ADS functionality
        self.declareProperty(MultipleFileProperty(name='File(s)',
                                                  extensions=[".raw", ".s0x"]),
                             direction=Direction.Input,
                             doc='Files of calibration runs (numors). Must be detector scans.')
        self.declareProperty(name='Foil type',
                             defaultValue='Hf01',
                             direction=Direction.Input,
                             doc="Hf01, Hf02, Ta10, Irp6, Iro5, Iro9")
        self.declareProperty(name='Calibration',
                             defaultValue=False,
                             direction=Direction.Input,
                             doc="Calibration flag, default is False in which case temperature measured")
        self.declareProperty(name='Reference temperature',
                             defaultValue='290',
                             direction=Direction.Input,
                             doc="Enter reference temperature in K")
        # TODO this should be removed
        # self.declareProperty('File extension', 'raw', doc="Enter file extension required eg. s01, s02, raw")
        self.declareProperty(name='Debug',
                             defaultvalue=False,
                             direction=Direction.Input,
                             doc="True/False - provides more verbose output of procedure for debugging purposes")
        
    def PyExec(self):
        # ----------------------------------------------------------
        # Imports resonance parameters from ResParam dictionary depending on the selected foil type.
        # ----------------------------------------------------------
        files = self.getProperty("File(s)").value
        foilType = self.getProperty("Foil type").value
        runNo = str(self.getProperty("Run number").value)
        isCalib = self.getProperty("Calibration").value
        mass = ResParamsDict[foilType + '_Mass']
        TD = ResParamsDict[foilType + '_TD'] # TODO what is TD?
        # Energy parameters are in eV
        energy = ResParamsDict[foilType + '_En']
        TwogG = ResParamsDict[foilType+'_TwogG']
        Gg = ResParamsDict[foilType+'_Gg']
        startE = ResParamsDict[foilType+'_startE']

        #Ediv = ResParamsDict[foilType+'_Ediv']
        Ediv = 0.0025  # TODO which is the correct Ediv? Hardcode or from ref

        endE = ResParamsDict[foilType+'_endE']
        refTemp = float(self.getProperty("Reference temperature").value) # TODO is this used?
        isDebug = self.getProperty("Debug").value # TODO is this used?

        #mon_pos= 13.43637495
        #tdiv = (0.5* mn * mon_pos**2)/(Ediv**2) # TODO are these used?

        monitorTransfit(runNo)
        
        # Define the gaussian width at the reference temperature
        # Factor of 1e3 converting to meV for Mantid

        width_300 = 1000.0 * np.sqrt(4*En*k*T_ref*mass/(e * (1+mass)**2)) #TODO what is 300? En?

        # ----------------------------------------------------------
        # Perform fits based on whether isCalib is flagged (calibration or measurement)
        # ----------------------------------------------------------
        if isCalib:
            lorentzFWHM = 1000.0 * (0.5 * TwogG + Gg)
            #
            # Old Voigt function as built into Mantid
            # Fit(Function='name=Voigt,LorentzAmp=-40.3305,LorentzPos=1096.81,LorentzFWHM=142.947,GaussianFWHM='+str(width_300)+',ties=(GaussianFWHM='+str(width_300)+');name=Polynomial,n=2,A0=110.767,A1=0.0270698,A2=-6.65217e-06', InputWorkspace=runNo+'_monitor', MaxIterations=1000, Output='S_fit', OutputCompositeMembers=True, StartX=600, EndX=1700)
            #
            # For guessing initial background values, read in y-data and use starting value as value for b0
            wsBgGuess = mtd[runNos + '_monitor']
            bg0guess = 0.86 * wsBgGuess.readY(0)[0]
            # Take peak position starting guess from tabulated value
            peakPosGuess = energy * 1000
            # New Voigt function as from Igor pro function
            # Fit(Function='name=TransVoigt,Position='+str(peakPosguess)+',Lorentzian FWHM='+str(Lorentz_FWHM)+',Gaussian FWHM='+str(width_300)+',Amplitude=2.07943,b0='+str(b0guess)+',b1=0.0063252,b2=0,constraints=(1<Lorentzian FWHM),ties=(Gaussian FWHM='+str(width_300)+', b2=0)', InputWorkspace=runNo+'_monitor', MaxIterations=20, Output='S_fit')
            # try fitting with Gaussian resolution component?
            Fit(Function='name=TransVoigt,Position=' + str(peakPosGuess) + ',Lorentzian FWHM=' + str(
                lorentzFWHM) + ',Gaussian FWHM=' + str(width_300) + ',Amplitude=1.6,b0=' + str(
                bg0guess) + ',b1=0.0063252,b2=0,constraints=(1<Lorentzian FWHM,' + str(width_300) + '<Gaussian FWHM)',
                InputWorkspace=runNos + '_monitor', MaxIterations=200, Output='S_fit')
            #
            DeleteWorkspace('S_fit_NormalisedCovarianceMatrix')
            DeleteWorkspace(runNos + '_monitor')
            # DeleteWorkspace('S_fit_Workspace')
            
        else:
            #self.log().information("Measuring temperature using S_paramTable for standard values.")
            S_fit= mtd['S_fit_Parameters']
            lorentzFWHM = S_fit.column(1)[1]
            gaussianFWHM = S_fit.column(1)[2]
            
            #calculates the gaussian resolution contribution, sometimes constraints are broken and the value can drop below that from width_300 alone, in which case the instrument contribution is set to zero.
            if gaussianFWHM > width_300:
                gaussianFWHM_inst = np.sqrt(gaussianFWHM**2 - width_300**2)
            else:
                gaussianFWHM_inst = 0
            #
            #Old Voigt function as built into Mantid
            #Fit(Function='name=Voigt,LorentzAmp='+str(S_fit.column(1)[0])+',LorentzPos='+str(S_fit.column(1)[1])+',LorentzFWHM='+str(Lorentz_FWHM)+',GaussianFWHM='+str(Gaussian_FWHM)+',ties=(LorentzFWHM='+str(Lorentz_FWHM)+');name=Polynomial,n=2,A0='+str(S_fit.column(1)[4])+',A1='+str(S_fit.column(1)[5])+',A2='+str(S_fit.column(1)[6]), InputWorkspace=runNo+'_monitor', MaxIterations=1000, Output='T_fit', OutputCompositeMembers=True, StartX=600, EndX=1700)
            #
            #New Voigt function as from Igor pro function
            Fit(Function='name=TransVoigt,Position='+str(S_fit.column(1)[0])+',Lorentzian FWHM='+str(lorentzFWHM)+',Gaussian FWHM='+str(GaussianFWHM)+',Amplitude='+str(S_fit.column(1)[3])+',b0='+str(S_fit.column(1)[4])+',b1='+str(S_fit.column(1)[5])+',b2='+str(S_fit.column(1)[6])+',constraints=('+str(gaussianFWHM)+'<Gaussian FWHM),ties=(Lorentzian FWHM='+str(lorentzFWHM)+')', InputWorkspace=runNos+'_monitor', MaxIterations=200, Output='T_fit')
            #
            DeleteWorkspace('T_fit_NormalisedCovarianceMatrix')
            DeleteWorkspace(runNos+'_monitor')
            #DeleteWorkspace('T_fit_Workspace')
            T_fit= mtd['T_fit_Parameters']
            Gaussian_FWHM_fitted = T_fit.column(1)[2]
            width_T = np.sqrt(Gaussian_FWHM_fitted**2 - Gaussian_FWHM_inst**2)
            #Factor of 1e-3 converts back from meV to eV
            Teff = (  ((width_T * 1e-3)**2) * e * ((1+Mass)**2)   )   /   (4*(1e-3)*T_fit.column(1)[0]*k*Mass)
            #err_d2En = (((Gaussian_FWHM_fitted * 1e-3)**2) / ((1e-3)*T_fit.column(1)[0]))* np.sqrt(  2*((T_fit.column(2)[2])/(Gaussian_FWHM_fitted))**2   + (T_fit.column(2)[0]/T_fit.column(1)[0])**2   )
            #errTeff = err_d2En * (e * (1+Mass)**2)  /   (4*k*Mass)
            Teff_low=(  (( (width_T - T_fit.column(2)[2]) * 1e-3)**2) * e * ((1+Mass)**2)   )   /   (4*(1e-3)*(T_fit.column(1)[0]+T_fit.column(2)[0])*k*Mass)
            Teff_high=(  (( (width_T + T_fit.column(2)[2]) * 1e-3)**2) * e * ((1+Mass)**2)   )   /   (4*(1e-3)*(T_fit.column(1)[0]-T_fit.column(2)[0])*k*Mass)
            errTeff = 0.5*(Teff_high-Teff_low)
        #----------------------------------------------------------
        #If the temperature is too far below the Debye temperature, then the result is inaccurate. Else the temperature is calcualted assuming free gas formulation
        #----------------------------------------------------------
            if 8*Teff < 3*TD:
                self.log().information("The effective temperature is currently too far below the Debye temperature to give an accurate measure.")
                Tactual = Teff
                Terror = errTeff
            else:
                Tactual = 3*TD / (4*np.log((8*Teff + 3*TD)/(8*Teff - 3*TD)))
                Tactual_high = 3*TD / (4*np.log((8*(Teff+errTeff) + 3*TD)/(8*(Teff+errTeff) - 3*TD)))
                Tactual_low = 3*TD / (4*np.log((8*(Teff-errTeff) + 3*TD)/(8*(Teff-errTeff) - 3*TD)))
                Terror = 0.5*(np.abs( Tactual - Tactual_high  ) + np.abs( Tactual - Tactual_low  ))
            #Introduce a catchment for unphysically small determined errors, setting to defualt value if below a set threshold
            Terror_flag=0
            if Terror<5.0:
                Terror_flag=1
                Terror=10.0
            if debug_flag == 'True':
                self.log().information("-----------------------------")
                self.log().information("Debugging....")
                self.log().information("The Debye temperature is "+str(TD)+" K")
                self.log().information("The effective temperature is: {:.1f}".format(Teff)+"+/- {:.1f}".format(errTeff)+" K")
                self.log().information("Energy bin width set to " + str(1000*Ediv)+" meV" )
                self.log().information("E range is between " + str(startE)+" and "+str(endE))
                self.log().information("Reference temperature is: {:.1f}".format(T_ref)+" K")
                self.log().information("Gaussian width at this reference temperature is: {:.2f}".format(width_300)+" meV")      
                self.log().information("Lorentzian FWHM is fixed: {:.2f}".format(Lorentz_FWHM)+" meV")
                self.log().information("Gaussian FWHM is fitted as: {:.2f}".format(Gaussian_FWHM_fitted)+" meV")
                self.log().information("Instrumental contribution is: {:.2f}".format(Gaussian_FWHM_inst)+" meV")
                self.log().information("Temperature contribution is: {:.2f}".format(width_T)+" meV")
                self.log().information("-----------------------------")
            
            
            self.log().information("Sample temperature is: {:.1f}".format(Tactual)+" +/- {:.1f}".format(Terror)+" K")
            if Terror_flag==1:
                self.log().information("(the default error, as determined error unphysically small)")\


    # ----------------------------------------------------------
    # Define function for importing raw monitor data, summing, normalising, converting units, and cropping
    # ----------------------------------------------------------

    def monitorTransfit(self, runNo):
        if "_" not in runNo:
            if len(str(runNo)) == 6:
                LoadRaw(Filename=directory + '/PEARL00' + runNo + '.' + extension, OutputWorkspace=runNo + '_raw')
            if len(str(runNo)) == 5:
                LoadRaw(Filename=directory + '/PEARL000' + runNo + '.' + extension, OutputWorkspace=runNo + '_raw')
            CropWorkspace(InputWorkspace=runNo + '_raw', OutputWorkspace=runNo + '_raw', XMin=100, XMax=19990)
            NormaliseByCurrent(InputWorkspace=runNo + '_raw', OutputWorkspace=runNo + '_raw')
            ExtractSingleSpectrum(InputWorkspace=runNo + '_raw', OutputWorkspace=runNo + '_3', WorkspaceIndex=3)
            DeleteWorkspace(runNo + '_raw')
            ConvertUnits(InputWorkspace=runNo + '_3', Target='Energy', OutputWorkspace=runNo + '_3')
            # Problem with rebinning?
            # Rebin(InputWorkspace=runNo+'_3', OutputWorkspace=runNo+'_3', Params=str(Ediv))
            TransfitRebin(runNo, runNo + '_monitor')
            # CropWorkspace(InputWorkspace=runNo+'_3', OutputWorkspace=runNo+'_monitor', XMin=1000*startE, XMax=1000*endE)
            DeleteWorkspace(runNo + '_3')
        else:
            start = int(re.split('_+', runNo)[0])
            end = int(re.split('_+', runNo)[1]) + 1
            for i in range(start, end):
                # LoadRaw(Filename='/Users/Chris/Documents/transfit/PEARL00'+str(i)+'.raw', OutputWorkspace=str(i)+'_raw')
                if len(str(i)) == 6:
                    LoadRaw(Filename=directory + '/PEARL00' + str(i) + '.' + extension,
                            OutputWorkspace=str(i) + '_raw')
                if len(str(i)) == 5:
                    LoadRaw(Filename=directory + '/PEARL000' + str(i) + '.' + extension,
                            OutputWorkspace=str(i) + '_raw')
                CropWorkspace(InputWorkspace=str(i) + '_raw', OutputWorkspace=str(i) + '_raw', XMin=100, XMax=19990)
                NormaliseByCurrent(InputWorkspace=str(i) + '_raw', OutputWorkspace=str(i) + '_raw')
                ExtractSingleSpectrum(InputWorkspace=str(i) + '_raw', OutputWorkspace=str(i) + '_3',
                                      WorkspaceIndex=3)
                DeleteWorkspace(str(i) + '_raw')
                ConvertUnits(InputWorkspace=str(i) + '_3', Target='Energy', OutputWorkspace=str(i) + '_3')
                # Problem with rebinning?
                # Rebin(InputWorkspace=str(i)+'_3', OutputWorkspace=str(i)+'_3', Params=str(Ediv))
                TransfitRebin(str(i), str(i) + '_3')
                # CropWorkspace(InputWorkspace=str(i)+'_3', OutputWorkspace=str(i)+'_3', XMin=1000*startE, XMax=1000*endE)
            for i in range(start + 1, end):
                Plus(LHSWorkspace=str(start) + '_3', RHSWorkspace=str(i) + '_3', OutputWorkspace=str(start) + '_3')
                DeleteWorkspace(str(i) + '_3')
            scale_val = (end - start - 1) ** (-1)
            # self.log().information("Scale val is "+str(scale_val))
            CreateSingleValuedWorkspace(OutputWorkspace='scale', DataValue=scale_val)
            Multiply(LHSWorkspace=str(start) + '_3', RHSWorkspace='scale', OutputWorkspace=runNo + '_monitor')
            DeleteWorkspace('scale')
            DeleteWorkspace(str(start) + '_3')


#Physical constants
k = 1.38066e-23
e = 1.60218e-19
Pi = np.pi
m = 1.e0
mn = 1.67493e-27
h = 6.62608e-34
Na = 6.02214e+23
flightpath = 13.805
t0 = 1598e-9
t_stnd = 19.0

# Resonance constants broadly consistent with those reported from
# 'Neutron cross sections Vol 1, Resonance Parameters'
# S.F. Mughabghab and D.I. Garber
# June 1973, BNL report 325
# Mass in atomic units, En in eV, temperatures in K, gamma factors in eV
ResParamsDict = {
    # Hf01
    "Hf01_Mass": 177.0,
    "Hf01_En": 1.098,
    "Hf01_TD": 252.0,
    "Hf01_TwogG": 0.00192,
    "Hf01_Gg": 0.0662,
    "Hf01_startE": 0.6,
    "Hf01_Ediv": 0.01,
    "Hf01_endE": 1.7,
    # Hf02
    "Hf02_Mass": 177.0,
    "Hf02_En":  2.388,
    "Hf02_TD": 252.0,
    "Hf02_TwogG": 0.009,
    "Hf02_Gg": 0.0608,
    "Hf02_startE": 2.0,
    "Hf02_Ediv": 0.01,
    "Hf02_endE": 2.7,
    # Ta10
    "Ta10_Mass": 181.0,
    "Ta10_En": 10.44,
    "Ta10_TD": 240.0,
    "Ta10_TwogG": 0.00335,
    "Ta10_Gg": 0.0069,
    "Ta10_startE": 9.6,
    "Ta10_Ediv": 0.01,
    "Ta10_endE": 11.4,
    # Irp6
    "Irp6_Mass": 191.0,
    "Irp6_En": 0.6528,
    "Irp6_TD": 420.0,
    "Irp6_TwogG": 0.000547,
    "Irp6_Gg": 0.072,
    "Irp6_startE": 0.1,
    "Irp6_Ediv": 0.01,
    "Irp6_endE": 0.9,
    # Iro5
    "Iro5_Mass": 191.0,
    "Iro5_En": 5.36,
    "Iro5_TD": 420.0,
    "Iro5_TwogG": 0.006,
    "Iro5_Gg": 0.082,
    "Iro5_startE": 4.9,
    "Iro5_Ediv": 0.01,
    "Iro5_endE": 6.3,
    # Iro9
    "Iro9_Mass": 191.0,
    "Iro9_En": 9.3,
    "Iro9_TD": 420.0,
    "Iro9_TwogG": 0.0031,
    "Iro9_Gg": 0.082,
    "Iro9_startE": 8.7,
    "Iro9_Ediv": 0.01,
    "Iro9_endE": 9.85
}

AlgorithmFactory.subscribe(Transfit)

