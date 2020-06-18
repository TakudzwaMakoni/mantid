# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
from mantid.simpleapi import Fit, DoublePulseFit, CompareWorkspaces, GausOsc, CreateWorkspace
from mantid.api import AnalysisDataService, IFunction1D, FunctionFactory
import numpy as np


class SingleDomainDoublePulseFitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        delta = 0.33
        x = np.linspace(0.,15.,100)
        x_offset = np.linspace(delta, 15. + delta, 100)

        testFunction = GausOsc()
        y1 = testFunction(x)
        y2 = testFunction(x_offset)
        y = y1+y2
        ws = CreateWorkspace(x,y)

        convolution =  FunctionFactory.createCompositeFunction('Convolution')
        innerFunction = FunctionFactory.createInitialized('name=GausOsc,A=0.2,Sigma=0.2,Frequency=1,Phi=0')
        deltaFunctions = FunctionFactory.createInitialized('(name=DeltaFunction,Height=1,Centre={},ties=(Height=1,Centre={});name=DeltaFunction,Height=1,Centre=0,ties=(Height=1,Centre=0))'.format(-delta, -delta))
        convolution.setAttributeValue('FixResolution', False)
        convolution.add(innerFunction)
        convolution.add(deltaFunctions)

        innerFunctionSingle = FunctionFactory.createInitialized('name=GausOsc,A=0.2,Sigma=0.2,Frequency=1,Phi=0')

        DoublePulseFit(Function=innerFunctionSingle, InputWorkspace=ws, CreateOutput = True, PulseOffset = delta, StartX=0.0, EndX=15.0, Output='DoublePulseFit', MaxIterations=1)
        Fit(Function=convolution, InputWorkspace=ws, CreateOutput = True, StartX=0.0, EndX=15.0, Output='Fit', MaxIterations=1)

    @classmethod
    def tearDownClass(cls):
        AnalysisDataService.clear()

    def test_that_simulated_output_data_is_the_same(self):
        result, message = CompareWorkspaces('Fit_Workspace', 'DoublePulseFit_Workspace')
        self.assertTrue(result)

    def test_that_covariance_matricies_are_the_same(self):
        result, message = CompareWorkspaces('Fit_NormalisedCovarianceMatrix', 'DoublePulseFit_NormalisedCovarianceMatrix')
        self.assertTrue(result)

    def test_that_output_parameters_are_the_same(self):
        result, message = CompareWorkspaces('Fit_Parameters', 'DoublePulseFit_Parameters')
        self.assertTrue(result)


class MultiDomainDoublePulseFitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        delta = 0.33
        x = np.linspace(0.,15.,100)
        x_offset = np.linspace(delta, 15. + delta, 100)

        testFunction = GausOsc()
        y1 = testFunction(x)
        y2 = testFunction(x_offset)
        y = y1+y2
        ws = CreateWorkspace(x,y)

        convolution =  FunctionFactory.createCompositeFunction('Convolution')
        innerFunction = FunctionFactory.createInitialized('name=GausOsc,A=0.2,Sigma=0.2,Frequency=1,Phi=0')
        deltaFunctions = FunctionFactory.createInitialized('(name=DeltaFunction,Height=1,Centre={},ties=(Height=1,Centre={});name=DeltaFunction,Height=1,Centre=0,ties=(Height=1,Centre=0))'.format(-delta, -delta))
        convolution.setAttributeValue('FixResolution', False)
        convolution.add(innerFunction)
        convolution.add(deltaFunctions)

        MultiDomainSingleFunction = FunctionFactory.createInitializedMultiDomainFunction('name=GausOsc,A=0.2,Sigma=0.2,Frequency=1,Phi=0', 2)
        MultiDomainConvolutionFunction = FunctionFactory.createInitializedMultiDomainFunction(str(convolution), 2)

        DoublePulseFit(Function=MultiDomainSingleFunction, InputWorkspace=ws, InputWorkspace_1=ws, CreateOutput = True, PulseOffset = delta, StartX=0.0, EndX=15.0, Output='DoublePulseFit', MaxIterations=1)
        Fit(Function=MultiDomainConvolutionFunction, InputWorkspace=ws, InputWorkspace_1=ws, CreateOutput = True, StartX=0.0, EndX=15.0, Output='Fit', MaxIterations=1)

    @classmethod
    def tearDownClass(cls):
        AnalysisDataService.clear()

    def test_that_simulated_output_data_is_the_same(self):
        result0, message0 = CompareWorkspaces('Fit_Workspace_0', 'DoublePulseFit_Workspace_0')
        result1, message1 = CompareWorkspaces('Fit_Workspace_1', 'DoublePulseFit_Workspace_1')

        self.assertTrue(result0)
        self.assertTrue(result1)

    def test_that_covariance_matricies_are_the_same(self):
        result, message = CompareWorkspaces('Fit_NormalisedCovarianceMatrix', 'DoublePulseFit_NormalisedCovarianceMatrix')
        self.assertTrue(result)

    def test_that_output_parameters_are_the_same(self):
        result, message = CompareWorkspaces('Fit_Parameters', 'DoublePulseFit_Parameters')
        self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()
