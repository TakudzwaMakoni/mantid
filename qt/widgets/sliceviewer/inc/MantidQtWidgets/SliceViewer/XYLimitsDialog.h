// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllOption.h"
#include "MantidGeometry/MDGeometry/IMDDimension.h"
#include "ui_XYLimitsDialog.h"
#include <QDialog>

/** Dialog that prompts the user
 * to enter the X and Y limits to the axes
 *
 * For use with the SliceViewer.
 */
class EXPORT_OPT_MANTIDQT_SLICEVIEWER XYLimitsDialog : public QDialog {
  Q_OBJECT

public:
  XYLimitsDialog(QWidget *parent = nullptr);
  ~XYLimitsDialog() override;

  void setXDim(const Mantid::Geometry::IMDDimension_const_sptr &dim);
  void setYDim(const Mantid::Geometry::IMDDimension_const_sptr &dim);
  void setLimits(double x0, double x1, double y0, double y1);
  double getXMin();
  double getXMax();
  double getYMin();
  double getYMax();

private:
  Ui::XYLimitsDialogClass ui;
};
