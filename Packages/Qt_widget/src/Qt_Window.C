// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : src/Qt_Window.C
// package       : QT_window
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_Window.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/IO/Qt_Window_tool.h>

namespace CGAL {

Qt_widget::Qt_widget(QWidget *parent, const char *name) :
  QWidget(parent, name), initialized(false), Locked(), _pointSize(4),
  _pointStyle(DISC), _has_tool(false), current_tool(0)
{
  setCaption("CGAL::Qt_widget");
  initialize();
  paint.begin(&pixmap);
  setBackgroundColor(Qt::white);
  paint.setPen(QPen(Qt::black,3));
  clear();
}

void Qt_widget::initialize()
{
  xmin=0;
  xmax=width()-1;
  ymin=0;
  ymax=height()-1;
  xscal=1;
  yscal=1;
  pixmap.resize(size());
}

void Qt_widget::setScales()
{
  xscal=(width()-1)/(xmax-xmin);
  yscal=(height()-1)/(ymax-ymin);
}

void Qt_widget::resizeEvent(QResizeEvent *e)
{
  // save paint state
  QFont f=paint.font();
  QBrush b=paint.brush();
  QPen p=paint.pen();
  QColor bc=paint.backgroundColor();

  paint.end();  // end painting on pixmap

  /*
    the only difference between an initialized Qt_widget and a
    non-initialized one is here:
    if the widget has been initialized, a resizeEvent modifies the
    scalings where as it modifies z_min(), z_max() dimensions if not.
  */
  if (!isInitialized())
    initialize();
  else
    {
      pixmap.resize(size());
      setScales();
    }
  paint.begin(&pixmap); // begin again painting on pixmap

  // restore paint state
  paint.setFont(f);
  paint.setBrush(b);
  paint.setPen(p);
  paint.setBackgroundColor(bc);

  clear();
  setScales();
  emit(resized());
}

void Qt_widget::paintEvent(QPaintEvent *e)
{
  // save paint state
  QFont f=paint.font();
  QBrush b=paint.brush();
  QPen p=paint.pen();
  QColor bc=paint.backgroundColor();

  paint.end();  // end painting on pixmap
  bitBlt(this, 0, 0, &pixmap); // copy pixmap to the Qt_widget
  paint.begin(&pixmap); // begin again painting on pixmap

  // restore paint state
  paint.setFont(f);
  paint.setBrush(b);
  paint.setPen(p);
  paint.setBackgroundColor(bc);
}

void Qt_widget::mousePressEvent(QMouseEvent *e)
{
  emit(mousePressed(e));
  if (has_tool())
    current_tool->mousePressEvent(e);
}

void Qt_widget::mouseReleaseEvent(QMouseEvent *e)
{
  emit(mouseReleased(e));
  if (has_tool())
    current_tool->mouseReleaseEvent(e);
}

void Qt_widget::mouseMoveEvent(QMouseEvent *e)
{
  emit(mouseMoved(e));
  if (has_tool())
    current_tool->mouseMoveEvent(e);
}

void Qt_widget::wheelEvent(QMouseEvent *e)
{
  if (has_tool())
    current_tool->wheelEvent(e);
}

void Qt_widget::mouseDoubleClickEvent(QMouseEvent *e)
{
  if (has_tool())
    current_tool->mouseDoubleClickEvent(e);
}

void Qt_widget::keyPressEvent(QKeyEvent *e)
{
  if (has_tool())
    current_tool->keyPressEvent(e);
}

void Qt_widget::keyReleaseEvent(QKeyEvent *e)
{
  if (has_tool())
    current_tool->keyReleaseEvent(e);
}

void Qt_widget::enterEvent(QEvent *e)
{
  if (has_tool())
    current_tool->enterEvent(e);
}

void Qt_widget::leaveEvent(QEvent *e)
{
  if (has_tool())
    current_tool->leaveEvent(e);
}

void Qt_widget::init(double x_min, double x_max, double y_min)
{
  double y_max=y_min+(height()-1)*(x_max-x_min)/(width()-1);
  init(x_min, x_max, y_min, y_max);
}

void Qt_widget::init(double  x_min, double x_max, double y_min, double y_max)
{
  xmin=x_min;
  xmax=x_max;
  ymin=y_min;
  ymax=y_max;
  setScales();
  initialized=true;
}

double Qt_widget::x_real(int x) const
{
  return(xmin+x/xscal);
}

double Qt_widget::y_real(int y) const
{
  return(ymax-y/yscal);
}

double Qt_widget::x_real_dist(double d) const
{
  return(d/xscal);
}

double Qt_widget::y_real_dist(double d) const
{
  return(d/yscal);
}

int Qt_widget::x_pixel(double x) const
{
  return( static_cast<int>((x-xmin)*xscal+0.5) );
}

int Qt_widget::y_pixel(double y) const
{
  return( - static_cast<int>((y-ymax)*yscal+0.5) );
}

int Qt_widget::x_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*xscal+0.5) );
  else
    return( static_cast<int>(d*xscal-0.5) );
}

int Qt_widget::y_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*yscal+0.5) );
  else
    return( static_cast<int>(d*yscal-0.5) );
}

Qt_widget& Qt_widget::operator<<(Qt_widget_tool* tool)
{
  if (has_tool())
    current_tool->detach();
  current_tool=tool;
  _has_tool=true;
  tool->attach(this);
  return *this;
}

Qt_widget& Qt_widget::operator<<(const Color& c)
{
  setColor(CGAL2Qt_Color(c));
  return *this;
}

Qt_widget& Qt_widget::operator<<(const PointStyle& ps)
{
  setPointStyle(ps);
  return *this;
}

void Qt_widget::clear() {
  painter().eraseRect(rect());
}

Qt_widget& operator<<(Qt_widget& w, const Bbox_2& r)
{
  int
    xmin = w.x_pixel(r.xmin()),
    ymin = w.y_pixel(r.ymin()),
    xmax = w.x_pixel(r.xmax()),
    ymax = w.y_pixel(r.ymax());

  w.painter().drawWinFocusRect(xmin, ymin, xmax-xmin, ymax-ymin);
  w.doPaint();
  return w;
}

} // namespace CGAL

#include "Qt_Window.moc"

#endif
