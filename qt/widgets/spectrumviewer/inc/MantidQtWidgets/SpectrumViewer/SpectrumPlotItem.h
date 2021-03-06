// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <QPainter>
#include <QRect>

#include <qwt_plot_item.h>
#include <qwt_scale_map.h>

#include "MantidQtWidgets/SpectrumViewer/DataArray.h"
#include "MantidQtWidgets/SpectrumViewer/DllOptionSV.h"

/**
    @class SpectrumPlotItem

    This class is responsible for actually drawing the image data onto
    a QwtPlot for the SpectrumView data viewer.

    @author Dennis Mikkelson
    @date   2012-04-03
 */

namespace MantidQt {
namespace SpectrumView {

class EXPORT_OPT_MANTIDQT_SPECTRUMVIEWER SpectrumPlotItem : public QwtPlotItem {

public:
  /// Construct basic plot item with NO data to plot.
  SpectrumPlotItem();

  ~SpectrumPlotItem() override;

  /// Specify the data to be plotted and the color table to use
  void setData(const DataArray_const_sptr &dataArray,
               std::vector<QRgb> *positiveColorTable,
               std::vector<QRgb> *negativeColorTable);

  /// Set a non-linear lookup table to scale data values before mapping to color
  void setIntensityTable(std::vector<double> *intensityTable);

  /// Draw the image (this is called by QWT and must not be called directly.)
  void draw(QPainter *painter, const QwtScaleMap &xMap, const QwtScaleMap &yMap,
            const QRect &canvasRect) const override;

protected:
  int m_bufferID;                    // set to 0 or 1 to select buffer
  DataArray_const_sptr m_dataArray0; // these provide double buffers
  DataArray_const_sptr m_dataArray1; // for the float data.

private:
  /* This class just uses the following */
  /* but they are created and deleted */
  /* in the upper level classes */
  std::vector<QRgb> *m_positiveColorTable;
  std::vector<QRgb> *m_negativeColorTable;
  std::vector<double> *m_intensityTable;
};

} // namespace SpectrumView
} // namespace MantidQt
