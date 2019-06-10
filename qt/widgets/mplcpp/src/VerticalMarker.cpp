// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidQtWidgets/MplCpp/VerticalMarker.h"
#include "MantidPythonInterface/core/CallMethod.h"
#include "MantidQtWidgets/Common/Python/QHashToDict.h"
#include "MantidQtWidgets/Common/Python/Sip.h"

using Mantid::PythonInterface::GlobalInterpreterLock;
// using Mantid::PythonInterface::PythonException;
using Mantid::PythonInterface::callMethodNoCheck;
using namespace MantidQt::Widgets::Common;
using namespace MantidQt::Widgets::MplCpp;

namespace {

Python::Object
newMarker(FigureCanvasQt *canvas, QString const &colour, double x,
          boost::optional<QHash<QString, QVariant>> const &otherKwargs) {
  GlobalInterpreterLock lock;

  Python::Object markersModule{Python::NewRef(
      PyImport_ImportModule("mantidqt.widgets.fitpropertybrowser.markers"))};

  auto const args = Python::NewRef(Py_BuildValue(
      "(Osd)", canvas->pyobj().ptr(), colour.toLatin1().constData(), x));
  Python::Dict kwargs = Python::qHashToDict(otherKwargs.get());

  auto const marker = markersModule.attr("VerticalMarker")(*args, **kwargs);
  return marker;
}

} // namespace

namespace MantidQt {
namespace Widgets {
namespace MplCpp {

/**
 * @brief Create a VerticalMarker instance
 * @param canvas The canvas to draw the vertical marker on to
 * @param colour The colour of the vertical marker
 * @param x The x coordinate of the marker
 */
VerticalMarker::VerticalMarker(FigureCanvasQt *canvas, QString const &colour,
                               double x,
                               QHash<QString, QVariant> const &otherKwargs)
    : InstanceHolder(newMarker(canvas, colour, x, otherKwargs)) {}

void VerticalMarker::redraw() { callMethodNoCheck<void>(pyobj(), "redraw"); }

void VerticalMarker::remove() { callMethodNoCheck<void>(pyobj(), "remove"); }

void VerticalMarker::setColor(QString const &color) {
  callMethodNoCheck<void>(pyobj(), "set_color", color.toLatin1().constData());
}

bool VerticalMarker::setXPosition(double x) {
  GlobalInterpreterLock lock;

  auto const movedPy = Python::Object(pyobj().attr("set_x_position")(x));
  return PyLong_AsLong(movedPy.ptr()) > 0;
}

void VerticalMarker::setXMinimum(double x) {
  callMethodNoCheck<void>(pyobj(), "set_x_minimum", x);
}

void VerticalMarker::setXMaximum(double x) {
  callMethodNoCheck<void>(pyobj(), "set_x_maximum", x);
}

void VerticalMarker::mouseMoveStart(double x, double y) {
  callMethodNoCheck<void>(pyobj(), "mouse_move_start", x, y);
}

void VerticalMarker::mouseMoveStop() {
  callMethodNoCheck<void>(pyobj(), "mouse_move_stop");
}

bool VerticalMarker::mouseMove(double x) {
  GlobalInterpreterLock lock;

  auto const movedPy = Python::Object(pyobj().attr("mouse_move")(x));
  return PyLong_AsLong(movedPy.ptr()) > 0;
}

bool VerticalMarker::isMoving() {
  GlobalInterpreterLock lock;

  auto const isMovingPy = Python::Object(pyobj().attr("is_marker_moving")());
  return PyLong_AsLong(isMovingPy.ptr()) > 0;
}

std::tuple<double, double>
VerticalMarker::transformPixelsToCoords(int xPixels, int yPixels) {
  GlobalInterpreterLock lock;

  auto const toDouble = [](Python::Object const &value) {
    return PyFloat_AsDouble(value.ptr());
  };

  auto const coords =
      pyobj().attr("transform_pixels_to_coords")(xPixels, yPixels);
  return std::make_tuple<double, double>(toDouble(coords[0]),
                                         toDouble(coords[1]));
}

// boost::optional<QScopedPointer<QCursor>>
// VerticalMarker::overrideCursor(double x, double y) {
//  GlobalInterpreterLock lock;
//
//  auto cursor = Python::Object(pyobj().attr("override_cursor")(x, y));
//
//  if (cursor.is_none())
//    return boost::optional<QScopedPointer<QCursor>>();
//  return QScopedPointer<QCursor>(Python::extract<QCursor>(cursor));
//}

} // namespace MplCpp
} // namespace Widgets
} // namespace MantidQt
