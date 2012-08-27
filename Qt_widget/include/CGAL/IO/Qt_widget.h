// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_WIDGET_H
#define CGAL_QT_WIDGET_H

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
//temporary, should remove next line!!
#include <CGAL/Triangle_2_Iso_rectangle_2_intersection.h>
#include <CGAL/IO/Color.h>

#include <vector>
#include <list>
#include <map>
#include <cmath>

#include <qwidget.h>
#include <qpainter.h>
#include <qcolor.h>
#include <qpixmap.h>
#include <qprinter.h>

#include <CGAL/auto_link/CGALQt3.h>

namespace CGAL {

class Qt_widget_layer;
enum PointStyle { PIXEL, CROSS, PLUS, CIRCLE, DISC, RECT, BOX };

class Qt_widget : public QWidget {
  Q_OBJECT
public:
  // constructor
  Qt_widget(QWidget *parent = 0, const char *name = 0);
  // destructor
  ~Qt_widget() {};

  // initialization of coordinates system
  void set_window(const double x_min,
		  const double x_max, 
		  const double y_min, 
		  const double y_max,
		  bool const_ranges = false);
  void zoom(double ratio);
  void zoom(double ratio, double xc, double yc);
  void set_x_scale(const double xscale){ xscal = xscale; }
  void set_y_scale(const double yscale){ yscal = yscale; }

  void move_center(const double distx, const double disty);
  void set_center(const double x, const double y);

  // painting system
  inline QPainter& get_painter() { return (*painter); };
  inline QPixmap& get_pixmap() { return (*pixmap); };
  inline QWMatrix& get_matrix() { return (*matrix); };
  void lock() { ++Locked; };
  void unlock() { if (Locked>0) --Locked; do_paint(); };
  void do_paint() { if (Locked==0) repaint(FALSE); };

  virtual QSize sizeHint() const {return QSize(geometry().width(),
					geometry().height());} 
  
  // properties
  // ~~~~~~~~~~
  // color
  QColor color() const;
  void setColor(const QColor c);
  // backGroundColor
  QColor backgroundColor() const;
  void setBackgroundColor(const QColor& c);
  // fillColor
  QColor fillColor() const;
  void setFillColor(const QColor c);
  // isFilled
  bool isFilled() const;
  void setFilled(const bool f);
  // lineWidth
  uint lineWidth() const;
  void setLineWidth(const uint i);
  // pointSize
  uint pointSize() const;
  void setPointSize(const uint i);
  // pointStyle
  typedef CGAL::PointStyle PointStyle;
  PointStyle pointStyle() const;
  void setPointStyle(const PointStyle s);
  // rasterOp
  RasterOp rasterOp() {return painter->rasterOp();}
  void setRasterOp(const RasterOp r) {painter->setRasterOp(r);}

  // CGAL version of setFooColor
  // used by the manipulators system
  // DO NOT USE THESE THREE UNDOCUMENTED FUNCTIONS !!
  inline void setColor(const Color c)
    { setColor(CGAL2Qt_Color(c)); };
  inline void setBackgroundColor(const Color c)
    { setBackgroundColor(CGAL2Qt_Color(c)); };
  inline void setFillColor(const Color c)
    { setFillColor(CGAL2Qt_Color(c)); };

  // set pen() color to c, cf. manipulators below for setting
  // backgroundColor and fillColor
  Qt_widget& operator<<(const Color& c);
  // set point style
  Qt_widget& operator<<(const PointStyle& ps);
  // clear the Widget, fill it with backgroundColor()
  void clear();

  // coordinates system
  // ~~~~~~~~~~~~~~~~~~
  // real world coordinates
  double x_real(int x) const;
  double y_real(int y) const;
  template <class FT>
  void x_real(int, FT&) const;
  template <class FT>
  void y_real(int y, FT&) const;


  double x_real_dist(double d) const;
  double y_real_dist(double d) const;


  // pixel coordinates
  int x_pixel(double x) const;
  int y_pixel(double y) const;
  int x_pixel_dist(double d) const;
  int y_pixel_dist(double d) const;

  inline double x_min() const { return xmin; };
  inline double y_min() const { return ymin; };
  inline double x_max() const { return xmax; };
  inline double y_max() const { return ymax; };

  inline double x_scal() { return xscal; }
  inline double y_scal() { return yscal; }

  void new_object(CGAL::Object obj) { emit(new_cgal_object(obj)); };
  
  //layers
  
  void attach(Qt_widget_layer *layer);
  
  
  // remove a layer from the list of displayable scenes
  void detach(Qt_widget_layer* s);

signals:
  void s_mousePressEvent(QMouseEvent *e);
  void s_mouseReleaseEvent(QMouseEvent *e);
  void s_mouseMoveEvent(QMouseEvent *e);
  void s_paintEvent(QPaintEvent *e);
  void s_resizeEvent(QResizeEvent *e);
  void s_wheelEvent(QWheelEvent *e);
  void s_mouseDoubleClickEvent(QMouseEvent *e);
  void s_keyPressEvent(QKeyEvent *e);
  void s_keyReleaseEvent(QKeyEvent *e);
  void s_enterEvent(QEvent *e);
  void s_leaveEvent(QEvent *e);
  void s_event(QEvent *e);
  
  void custom_redraw(); //deprecated:  if user want to draw something
                        //after layers replaced by redraw_on_front
  void redraw_on_front(); //called by redraw at the end
  void redraw_on_back();  //called by redraw at the beginning


  void new_cgal_object(CGAL::Object);	//this signal is emited every time an
					//attached tool constructed an object

  void rangesChanged(); 
  // triggered when ranges (xmin, xmax, ymin,...) are changed

public slots:
  void print_to_ps();
  virtual void redraw();

// backward-compatibility with CGAL-2.4, back() and forth() are
// deprecated, as well as add_to_history() or clear_history().
signals:
  void internal_back();
  void internal_forth();
  void internal_add_to_history();
  void internal_clear_history();
public slots:
  bool back() { emit(internal_back()); return true; }
  bool forth() { emit(internal_forth()); return true; }
public:
  void add_to_history() { emit(internal_add_to_history()); }
  void clear_history() { emit(internal_clear_history()); }

protected:
  void paintEvent(QPaintEvent *e);
  void resizeEvent(QResizeEvent *e);
  void showEvent(QShowEvent *e);
  void mousePressEvent(QMouseEvent *e);
  void mouseReleaseEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void wheelEvent(QWheelEvent *e);
  void mouseDoubleClickEvent(QMouseEvent *e);
  void keyPressEvent(QKeyEvent *e);
  void keyReleaseEvent(QKeyEvent *e);
  void enterEvent(QEvent *e);
  void leaveEvent(QEvent *e);
  bool event(QEvent *e);


private:
  // private functions
  // ~~~~~~~~~~~~~~~~~

  void resize_pixmap();
  // resize properly the pixmap size, saving then restoring the
  // painter properties

  void	  set_scales(); 
  // set xscal and yscal. Update ranges if const_ranges is false.

  // color types convertors
  static QColor CGAL2Qt_Color(Color c);
  static Color Qt2CGAL_color(QColor c);

  void attach_standard(Qt_widget_layer *layer);
  bool is_standard_active();
  bool does_standard_eat_events();
  friend class Qt_widget_standard_toolbar;


  // private member datas
  // ~~~~~~~~~~~~~~~~~~~~
  bool    set_scales_to_be_done;
  // this flag is set when the widget is not visible and should
  // postpone the set_scales() call.

  unsigned int Locked;
  // point style and size
  uint	      _pointSize;
  PointStyle  _pointStyle;

  QPixmap     *pixmap;	// the pixmap on which paints the painter
  QPainter    *painter;	// the painter
  QPrinter    *printer;	// the printer
  QWMatrix    *matrix;  // the world matrix

  QBrush      savedBrush; // saved brush, to be able to restore it on
  // setFilled(true)

  double    xmin, xmax, ymin, ymax; // real dimensions
  double    xmin_old, xmax_old, ymin_old, ymax_old;
            //backup ranges for resize
  double    xscal, yscal; // scales int/double
  bool      constranges; // tell if the ranges should be const

  //for layers
  std::list<Qt_widget_layer*>	qt_layers;
  std::list<Qt_widget_layer*> qt_standard_layers;
};//end Qt_widget class

// manipulators
// ~~~~~~~~~~~~
// single manipulators
inline
Qt_widget& operator<<(Qt_widget& w, Qt_widget& (*m)(Qt_widget&))
{
  return m(w);
}

// w << noFill << ... stop the filling of geometrical object
inline
Qt_widget& noFill(Qt_widget& w)
{
  w.setFilled(false);
  return w;
}

// manipulators with one argument
template <class Param>
struct Qt_widgetManip {
  Qt_widget& (*f)(Qt_widget&, Param);
  Param p;
  Qt_widgetManip(Qt_widget& (*ff)(Qt_widget&, Param),
		 Param pp) : f(ff), p(pp) {}
};

// usage: w << manip(Param) f ...
template <class Param>
Qt_widget& operator<<(Qt_widget& w, Qt_widgetManip<Param> m)
{
  return m.f(w, m.p);
}

#define CGAL_QTWIDGET_MANIP(param,function) \
inline \
Qt_widget& __Qt_widgetManip##function##Aux (Qt_widget& w, param p) \
{ w.set##function(p); return w; } \
inline \
Qt_widgetManip<param> function(param p) \
{ return Qt_widgetManip<param>( __Qt_widgetManip##function##Aux, p); }

// w << BackgroundColor(c) << ... sets the background color
CGAL_QTWIDGET_MANIP( Color, BackgroundColor )

// w << FillColor(c) << ... sets the fill color
CGAL_QTWIDGET_MANIP( Color, FillColor )

// w << LineWidth(i) << ... sets lines width
CGAL_QTWIDGET_MANIP( unsigned int, LineWidth )

// w << PointSize(i) << ... sets points size
CGAL_QTWIDGET_MANIP( unsigned int, PointSize )

// color types convertors
// ~~~~~~~~~~~~~~~~~~~~~~
inline
QColor Qt_widget::CGAL2Qt_Color(Color c)
{
  return QColor(c.red(), c.green(), c.blue());
}

inline
Color Qt_widget::Qt2CGAL_color(QColor c)
{
  return Color(c.red(),c.green(),c.blue());
}

// properties
// ~~~~~~~~~~
inline
QColor Qt_widget::color() const
{
  return painter->pen().color();
}


inline
void Qt_widget::setColor(const QColor c)
{
  QPen p=get_painter().pen();
  p.setColor(c);
  get_painter().setPen(p);
}

inline
QColor Qt_widget::backgroundColor() const
{
  return painter->backgroundColor();
}

inline
void Qt_widget::setBackgroundColor(const QColor& c)
{
  QWidget::setPaletteBackgroundColor(c);
  get_painter().setBackgroundColor(c);
  clear();
}

inline
QColor Qt_widget::fillColor() const
{
  return painter->brush().color();
}

inline
void Qt_widget::setFillColor(const QColor c)
{
  setFilled(true);
  get_painter().setBrush(c);
}

inline
bool Qt_widget::isFilled() const
{
  return( painter->brush().style()==Qt::NoBrush );
}

inline
void Qt_widget::setFilled(const bool f)
{
  if (f)
    painter->setBrush(savedBrush);
  else
    {
      savedBrush=painter->brush();
      painter->setBrush(QBrush());
    };
}

inline
uint Qt_widget::lineWidth() const
{
  return( painter->pen().width());
}

inline
void Qt_widget::setLineWidth(const unsigned int i)
{
  QPen p=get_painter().pen();
  p.setWidth(i);
  get_painter().setPen(p);
}

inline
uint Qt_widget::pointSize() const
{
  return _pointSize;
}

inline
void Qt_widget::setPointSize(const unsigned int i)
{
  _pointSize=i;
}

inline
PointStyle Qt_widget::pointStyle() const
{
  return _pointStyle;
}

inline
void Qt_widget::setPointStyle(const PointStyle ps)
{
  _pointStyle=ps;
}

// drawing methods
// ~~~~~~~~~~~~~~~

template <class R>
Qt_widget& operator<<(Qt_widget& w, const Point_2<R>& p)
{
  int x = w.x_pixel(CGAL::to_double(p.x()));
  int y = w.y_pixel(CGAL::to_double(p.y()));

  uint size=w.pointSize();
  PointStyle ps=w.pointStyle();

  switch (ps)
  {
    case PIXEL:
    {
       w.get_painter().drawPoint(x,y);
       break;
    }
    case CROSS:
    {
       w.get_painter().drawLine(x-size/2, y-size/2, x+size/2, y+size/2);
       w.get_painter().drawLine(x-size/2, y+size/2, x+size/2, y-size/2);
       break;
    }
    case PLUS:
    {
       w.get_painter().drawLine(x, y-size/2, x, y+size/2);
       w.get_painter().drawLine(x-size/2, y, x+size/2, y);
       break;
    }
    case CIRCLE:
    {
       QBrush old_brush=w.get_painter().brush();
       w.get_painter().setBrush(QBrush());
       w.get_painter().drawEllipse(x-size/2, y-size/2, size, size);
       w.get_painter().setBrush(old_brush);
       break;
    }
    case DISC:
    {
       QBrush old_brush=w.get_painter().brush();
       w.get_painter().setBrush(w.get_painter().pen().color());
       w.get_painter().drawEllipse(x-size/2, y-size/2, size, size);
       w.get_painter().setBrush(old_brush);
       break;
    }
    case RECT:
    {
      QBrush old_brush=w.get_painter().brush();
      w.get_painter().setBrush(QBrush());
      w.get_painter().drawRect(x-size/2, y-size/2, size, size);
      w.get_painter().setBrush(old_brush);
      break;
    }
    case BOX:
    {
      QBrush old_brush=w.get_painter().brush();
      w.get_painter().setBrush(w.get_painter().pen().color());
      w.get_painter().drawRect(x-size/2, y-size/2, size, size);
      w.get_painter().setBrush(old_brush);
      break;
    }
  };
  w.do_paint();
  return w;
}

#ifdef CGAL_SEGMENT_2_H
template <class R>
Qt_widget& operator<<(Qt_widget& w, const Segment_2<R>& s)
{
  typedef Simple_cartesian<double> RT;

  double xr1, yr1, xr2, yr2;
  double scs_x, scs_y, sct_x, sct_y;
  scs_x = CGAL::to_double(s.source().x());
  scs_y = CGAL::to_double(s.source().y());
  sct_x = CGAL::to_double(s.target().x());
  sct_y = CGAL::to_double(s.target().y());

  xr1 = w.x_real(0); xr2 = w.x_real(w.geometry().width());
  //next condition true if is outside on the X axes
  if((scs_x < xr1 && sct_x < xr1) ||
     (scs_x > xr2 && sct_x > xr2))
    return w;
  else{
    yr2 = w.y_real(0); yr1 = w.y_real(w.geometry().height());
    //next condition true if is outside on the Y axes
    if((scs_y < yr1 && sct_y < yr1) ||
       (scs_y > yr2 && sct_y > yr2))
      return w;
  }
  
  //if is here, the segment intersect the screen boundaries or is inside
  int x1, y1, x2, y2;
  Segment_2<RT>  sr;
  sr = Segment_2<RT>(Point_2<RT>(scs_x, scs_y), Point_2<RT>(sct_x, sct_y));
  //next condition true if the segment is inside
  if(!(scs_x >= xr1 && scs_x <= xr2 &&
     sct_x >= xr1 && sct_x <= xr2 && 
     scs_y >= yr1 && scs_y <= yr2 &&
     sct_y >= yr1 && sct_y <= yr2))
    {
    Iso_rectangle_2<RT> r = Iso_rectangle_2<RT>(Point_2<RT>(xr1, yr1),
                                              Point_2<RT>(xr2, yr2));
    CGAL::Object obj = CGAL::intersection(r, sr);  
    if (const Point_2<RT> *p = object_cast<Point_2<RT> >(&obj)){
      return w << *p;
    }
    else if (const Segment_2<RT> *s = object_cast<Segment_2<RT> >(&obj)) {
      sr = *s;
    }
    else {
      CGAL_assertion(obj.is_empty());
      return w;
    }
  }
  x1 = w.x_pixel(CGAL::to_double(sr.source().x()));
  x2 = w.x_pixel(CGAL::to_double(sr.target().x()));
  y1 = w.y_pixel(CGAL::to_double(sr.source().y()));
  y2 = w.y_pixel(CGAL::to_double(sr.target().y()));
  w.get_painter().drawLine(x1, y1, x2, y2);
  w.do_paint();
  return w;
}
#endif // CGAL_SEGMENT_2_H

#ifdef CGAL_LINE_2_H

template <class R>
Qt_widget& operator<<(Qt_widget& w, const Line_2<R>& l)
{
  typedef Simple_cartesian<double> Rep;
  typedef Point_2<Rep> Point;

  const Point_2<R>
    p1=l.point(),
    p2=p1+l.direction().vector();

  const Point
    p1d=Point(CGAL::to_double(p1.x()),CGAL::to_double(p1.y())),
    p2d=Point(CGAL::to_double(p2.x()),CGAL::to_double(p2.y()));

  double
    x1=w.x_min(),
    y1=w.y_min(),
    x2=w.x_max(),
    y2=w.y_max();

  const double
    dx=p1d.x()-p2d.x(),
    dy=p1d.y()-p2d.y();

  if (dx==0 && dy==0) return w;

  if (CGAL::abs(dx)>CGAL::abs(dy))
    {
      y1=p1d.y()+(x1-p1d.x())*dy/dx;
      y2=p1d.y()+(x2-p1d.x())*dy/dx;
    }
  else
    {
      x1=p1d.x()+(y1-p1d.y())*dx/dy;
      x2=p1d.x()+(y2-p1d.y())*dx/dy;
    }

  w.get_painter().drawLine(w.x_pixel(x1),w.y_pixel(y1),
		       w.x_pixel(x2),w.y_pixel(y2));
  return w;
}

#endif // CGAL_LINE_2_H

#ifdef CGAL_RAY_2_H
template <class R>
Qt_widget& operator<<(Qt_widget& w, const Ray_2<R>& r)
{
  typedef Simple_cartesian<double> Rep;
  typedef Point_2<Rep> Point;

  const Point_2<R>
    p1=r.point(0),
    p2=r.point(1);

  const Point
    p1d=Point(CGAL::to_double(p1.x()),CGAL::to_double(p1.y())),
    p2d=Point(CGAL::to_double(p2.x()),CGAL::to_double(p2.y()));


  const double
    dx=p1d.x()-p2d.x(),
    dy=p1d.y()-p2d.y();

  if (dx==0 && dy==0) return w;

  double x,y;

  if (CGAL::abs(dx)>CGAL::abs(dy))
    {
      if (p1d.x()<p2d.x())
	x = w.x_max();
      else
	x = w.x_min();
      y=p1d.y()+(x-p1d.x())*dy/dx;
    }
  else
    {
      if (p1d.y()<p2d.y())
	y = w.y_max();
      else
	y = w.y_min();
      x=p1d.x()+(y-p1d.y())*dx/dy;
    }
  w.get_painter().drawLine(w.x_pixel(p1d.x()),w.y_pixel(p1d.y()),
		       w.x_pixel(x),w.y_pixel(y));
  return w;

}
#endif //CGAL_RAY_2_H

#ifdef CGAL_TRIANGLE_2_H
template< class R >
Qt_widget&
operator<<(Qt_widget& w, const Triangle_2<R>& t)
{
  CGAL::Iso_rectangle_2<R> r( Point_2<R>(w.x_real(0), w.y_real(0)), 
                              Point_2<R>(w.x_real(w.geometry().width()), 
                              w.y_real(w.geometry().height())));
  CGAL::Object obj = CGAL::intersection(t, r);
  Point_2<R> pi;
  Segment_2<R> si;
  Triangle_2<R> ti;
  typedef Point_2<R> Point;
  std::vector<Point> vi;
  if(CGAL::assign(pi, obj))
    w << pi;
  if(CGAL::assign(si, obj))
    w << si;
  if(CGAL::assign(ti, obj))
  {
    QPointArray array(3);
    array[0] = QPoint(w.x_pixel(CGAL::to_double(t.vertex(0).x())), 
                                w.y_pixel(CGAL::to_double(t.vertex(0).y())));
    array[1] = QPoint(w.x_pixel(CGAL::to_double(t.vertex(1).x())), 
                                w.y_pixel(CGAL::to_double(t.vertex(1).y())));
    array[2] = QPoint(w.x_pixel(CGAL::to_double(t.vertex(2).x())), 
                                w.y_pixel(CGAL::to_double(t.vertex(2).y())));
    w.get_painter().drawPolygon(array);
  }   
  if(CGAL::assign(vi, obj)){
    QPointArray array(int(vi.size()));
    typename std::vector<Point>::const_iterator it = vi.begin();
    int pos = 0;
    while(it != vi.end()){
      array[pos] = QPoint(w.x_pixel(CGAL::to_double((*it).x())), 
                          w.y_pixel(CGAL::to_double((*it).y())));
      pos++;
      it++;
    }  
  w.get_painter().drawPolygon(array);
  }
  w.do_paint();

  return w;}
#endif

#ifdef CGAL_CIRCLE_2_H
template < class R>
Qt_widget& operator<<(Qt_widget& w, const Circle_2<R>& c)
{
  int 
    cx=w.x_pixel(CGAL::to_double(c.center().x())),
    cy=w.y_pixel(CGAL::to_double(c.center().y())),
    rx=w.x_pixel_dist((std::sqrt(CGAL::to_double(c.squared_radius())))),
    ry=w.y_pixel_dist((std::sqrt(CGAL::to_double(c.squared_radius()))));

  w.get_painter().drawEllipse(cx-rx,cy-ry,2*rx,2*ry);
  w.do_paint();
  return w;
}
#endif // CGAL_CIRCLE_2_H

#ifdef CGAL_ISO_RECTANGLE_2_H
template< class R >
Qt_widget&
operator<<(Qt_widget& w, const Iso_rectangle_2<R>& r)
{
  int xmin = w.x_pixel(CGAL::to_double(r.xmin()));
  int ymin = w.y_pixel(CGAL::to_double(r.ymin()));
  int xmax = w.x_pixel(CGAL::to_double(r.xmax()));
  int ymax = w.y_pixel(CGAL::to_double(r.ymax()));
  w.get_painter().drawRect(xmin,ymin,xmax-xmin,ymax-ymin);
  w.do_paint();
  return w;
}
#endif // CGAL_ISO_RECTANGLE_2_H

#ifdef CGAL_BBOX_2_H
Qt_widget& operator<<(Qt_widget& w, const Bbox_2& r);
// see Qt_widget for the implementation of this non-template function
#endif // CGAL_BBOX_2_H

// templated x_real and y_real

template <class FT>
void Qt_widget::x_real(int x, FT& return_t) const
{
  if(xscal<1)
    return_t = static_cast<FT>(xmin+(int)(x/xscal));
  else{
    return_t = static_cast<FT>(xmin+x/xscal);
  }
}

template <class FT>
void Qt_widget::y_real(int y, FT& return_t) const
{
    if(yscal<1)
      return_t = static_cast<FT>(ymax-(int)(y/yscal));
    else{
      return_t = static_cast<FT>(ymax-y/yscal);
  }  
}

} // namespace CGAL

#endif // CGAL_QT_WIDGET_H
