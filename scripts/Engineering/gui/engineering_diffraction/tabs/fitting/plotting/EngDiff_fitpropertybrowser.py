# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
#

from qtpy.QtCore import Slot
from mantidqt.utils.qt import import_qt
from mantidqt.utils.observer_pattern import GenericObservable
from mantidqt.widgets.fitpropertybrowser import FitPropertyBrowser
from mantid.api import AnalysisDataService

BaseBrowser = import_qt('.._common', 'mantidqt.widgets', 'FitPropertyBrowser')


class EngDiffFitPropertyBrowser(FitPropertyBrowser):
    """
    A wrapper around python FitPropertyBrowser altered for
    engineering diffraction UI
    """

    def __init__(self, canvas, toolbar_manager, parent=None):
        super(EngDiffFitPropertyBrowser, self).__init__(canvas, toolbar_manager, parent)
        self.fit_notifier = GenericObservable()

    def set_output_window_names(self):
        """
         Override this function as no window
        :return: None
        """
        return None

    @Slot(str)
    def fitting_done_slot(self, name):
        """
        This is called after Fit finishes to update the fit curves.
        :param name: The name of Fit's output workspace.
        """

        ws = AnalysisDataService.retrieve(name)
        self.do_plot(ws, plot_diff=self.plotDiff())
        # self.loadFunction("name=BackToBackExponential,I=1000,A=0.0361794,B=0.0215356,X0=38821.8,S=30.5841") # I=1677.09
        self.saveFunction(self.workspaceName())
        results_dict = eval(self.getFitAlgorithmParameters().replace('true', 'True').replace('false', 'False'))
        self.fit_notifier.notify_subscribers(results_dict)
