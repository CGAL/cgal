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
// file          : include/CGAL/IO/Qt_widget.h
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_H
#define CGAL_QT_WIDGET_H

#include <qwidget.h>
#include <qpainter.h>
#include <qcolor.h>
#include <qpixmap.h>
#include <qmessagebox.h>
#include <qprinter.h>


#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <vector>
#include <list>
#include <map>

#ifndef CGAL_QT_WIDGET_HISTORY_H
  #include <CGAL/IO/Qt_widget_history.h>
#endif
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
  void set_window(double x_min, double x_max, double y_min, double y_max,
		  bool const_ranges = false);
  void zoom_in(double ratio);
  void zoom_out(double ratio);
  void zoom_in(double ratio, double xc, double yc);
  void zoom_out(double ratio, double xc, double yc);
  void set_x_scale(double xscale){ xscal = xscale; }
  void set_y_scale(double yscale){ yscal = yscale; }

  void add_to_history(){history.add_to_history(xmin, xmax, ymin, ymax);}
  void clear_history(){history.clear();};

  inline void move_center(double distx, double disty) 
  {    
    xcentre += distx;
    ycentre += disty;
    set_scale_center(xcentre, ycentre);
    add_to_history(); //add the current viewport to history
    configure_history_buttons();
  }
  inline void set_center(double x, double y) 
  {
    xcentre = x; ycentre = y;
    set_scale_center(xcentre, ycentre);
  };

  // painting system
  inline QPainter& get_painter() { return (*painter); };
  inline QPixmap& get_pixmap() { return (*pixmap); };
  inline QWMatrix& get_matrix() { return (*matrix); };
  void lock() { ++Locked; };
  void unlock() { if (Locked>0) --Locked; do_paint(); };
  void do_paint() { if (Locked==0) repaint( FALSE ); };
  
  
  // properties
  // ~~~~~~~~~~
  // color
  QColor color() const;
  void setColor(QColor c);
  // backGroundColor
  QColor backgroundColor() const;
  void setBackgroundColor(QColor c);
  // fillColor
  QColor fillColor() const;
  void setFillColor(QColor c);
  // isFilled
  bool isFilled() const;
  void setFilled(bool f);
  // lineWidth
  uint lineWidth() const;
  void setLineWidth(uint i);
  // pointSize
  uint pointSize() const;
  void setPointSize(uint i);
  // pointStyle
  typedef CGAL::PointStyle PointStyle;
  PointStyle pointStyle() const;
  void setPointStyle(PointStyle s);
  // rasterOp
  RasterOp rasterOp() {return painter->rasterOp();}
  void setRasterOp(RasterOp r) {painter->setRasterOp(r);}

  // CGAL version of setFooColor
  // used by the manipulators system
  // DO NOT USE THESE THREE UNDOCUMENTED FUNCTIONS !!
  inline void setColor(Color c)
    { setColor(CGAL2Qt_Color(c)); };
  inline void setBackgroundColor(Color c)
    { setBackgroundColor(CGAL2Qt_Color(c)); };
  inline void setFillColor(Color c)
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
  // real coordinates
  double x_real(int x) const;
  double y_real(int y) const;
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
  void s_wheelEvent(QMouseEvent *e);
  void s_mouseDoubleClickEvent(QMouseEvent *e);
  void s_keyPressEvent(QKeyEvent *e);
  void s_keyReleaseEvent(QKeyEvent *e);
  void s_enterEvent(QEvent *e);
  void s_leaveEvent(QEvent *e);
  void s_event(QEvent *e);
  
  void custom_redraw(); // if user want to draw something after layers
  void new_cgal_object(CGAL::Object);	//this signal is emited every time an
					//attached tool constructed an object
//private signals (not documented)
  void set_back_enabled(bool i);    //used by the standard toolbar
  void set_forward_enabled(bool i); //used by the standard toolbar

public slots:
  void print_to_ps();
  virtual void redraw();
  void back(){
    if(history.back()){
      set_window_p(history.get_atom()->x1(), history.get_atom()->x2(),
                history.get_atom()->y1(), history.get_atom()->y2());
      configure_history_buttons();
    }
  }
  void forth(){
    if(history.forward()) {
      set_window_p(history.get_atom()->x1(), history.get_atom()->x2(),
                history.get_atom()->y1(), history.get_atom()->y2());
      configure_history_buttons();
    }
  }
protected:
  void paintEvent(QPaintEvent *e);
  void resizeEvent(QResizeEvent *e);
  void mousePressEvent(QMouseEvent *e);
  void mouseReleaseEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void wheelEvent(QMouseEvent *e);
  void mouseDoubleClickEvent(QMouseEvent *e);
  void keyPressEvent(QKeyEvent *e);
  void keyReleaseEvent(QKeyEvent *e);
  void enterEvent(QEvent *e);
  void leaveEvent(QEvent *e);
  bool event(QEvent *e);


private:
  void	  set_scales(); // set xscal and yscal
  void	  set_scale_center(double xc, double yc);
  double  xcentre, ycentre; //the center of the axex
  
  Qt_widget_history history;
  void set_window_p(double x_min, double x_max, double y_min, double y_max);
  void configure_history_buttons(); //change the enabled state of
                                    //the history buttons

  // color types convertors
  static QColor CGAL2Qt_Color(Color c);
  static Color Qt2CGAL_color(QColor c);

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

  double xmin, xmax, ymin, ymax; // real dimensions
  double xscal, yscal; // scalings int/double
  bool constranges; // tell if the ranges should be const
  bool is_the_first_time;

  //for layers
  std::list<Qt_widget_layer*>	qt_layers;
  std::list<Qt_widget_layer*> qt_standard_layers;
  void attach_standard(Qt_widget_layer *layer);
  bool is_standard_active();
  friend class Qt_widget_standard_toolbar;
};//end Qt_widget class

// manipulators
// ~~~~~~~~~~~~
// single manipulators
inline
Qt_widget& operator<<(Qt_widget& w, Qt_widget& (*m)(Qt_widget&))
{
  return m(w);
};

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

// usage: w << manip(Param) << ...
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
};


inline
void Qt_widget::setColor(QColor c)
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
void Qt_widget::setBackgroundColor(QColor c)
{
  QWidget::setBackgroundColor(c);
  get_painter().setBackgroundColor(c);
  clear();
}

inline
QColor Qt_widget::fillColor() const
{
  return painter->brush().color();
}

inline
void Qt_widget::setFillColor(QColor c)
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
void Qt_widget::setFilled(bool f)
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
void Qt_widget::setLineWidth(unsigned int i)
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
void Qt_widget::setPointSize(unsigned int i)
{
  _pointSize=i;
}

inline
PointStyle Qt_widget::pointStyle() const
{
  return _pointStyle;
}

inline
void Qt_widget::setPointStyle(PointStyle ps)
{
  _pointStyle=ps;
}

// drawing methods
// ~~~~~~~~~~~~~~~

template <class R>
Qt_widget& operator<<(Qt_widget& w, const Point_2<R>& p)
{
  int x = w.x_pixel(to_double(p.x()));
  int y = w.y_pixel(to_double(p.y()));

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
  const int
    x1=w.x_pixel(to_double(s.source().x())),
    y1=w.y_pixel(to_double(s.source().y())),
    x2=w.x_pixel(to_double(s.target().x())),
    y2=w.y_pixel(to_double(s.target().y()));
  w.get_painter().drawLine(x1,y1,x2,y2);
  w.do_paint();
  return w;
}
#endif // CGAL_SEGMENT_2_H

#ifdef CGAL_LINE_2_H

template <class R>
Qt_widget& operator<<(Qt_widget& w, const Line_2<R>& l)
{
  typedef Cartesian<double> Rep;
  typedef Point_2<Rep> Point;

  const Point_2<R>
    p1=l.point(),
    p2=p1+l.direction().vector();

  const Point
    p1d=Point(to_double(p1.x()),to_double(p1.y())),
    p2d=Point(to_double(p2.x()),to_double(p2.y()));

  double
    x1=w.x_min(),
    y1=w.y_min(),
    x2=w.x_max(),
    y2=w.y_max();

  const double
    dx=p1d.x()-p2d.x(),
    dy=p1d.y()-p2d.y();

  if (dx==0 && dy==0) return w;

  if (fabs(dx)>fabs(dy))
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
  typedef Cartesian<double> Rep;
  typedef Point_2<Rep> Point;

  const Point_2<R>
    p1=r.point(0),
    p2=r.point(1);

  const Point
    p1d=Point(to_double(p1.x()),to_double(p1.y())),
    p2d=Point(to_double(p2.x()),to_double(p2.y()));


  const double
    dx=p1d.x()-p2d.x(),
    dy=p1d.y()-p2d.y();

  if (dx==0 && dy==0) return w;

  double x,y;

  if (fabs(dx)>fabs(dy))
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
  const int
    ax = w.x_pixel(to_double(t.vertex(0).x())),
    ay = w.y_pixel(to_double(t.vertex(0).y())),
    bx = w.x_pixel(to_double(t.vertex(1).x())),
    by = w.y_pixel(to_double(t.vertex(1).y())),
    cx = w.x_pixel(to_double(t.vertex(2).x())),
    cy = w.y_pixel(to_double(t.vertex(2).y()));

  QPointArray array;

  array.setPoints(3,ax,ay,bx,by,cx,cy);
  w.get_painter().drawPolygon(array);
  w.do_paint();
  return w;
}
#endif

#ifdef CGAL_CIRCLE_2_H
template < class R>
Qt_widget& operator<<(Qt_widget& w, const Circle_2<R>& c)
{
  int 
    cx=w.x_pixel(to_double(c.center().x())),
    cy=w.y_pixel(to_double(c.center().y())),
    rx=w.x_pixel_dist((std::sqrt(to_double(c.squared_radius())))),
    ry=w.y_pixel_dist((std::sqrt(to_double(c.squared_radius()))));

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
  int
    xmin = w.x_pixel(to_double(r.min().x())),
    ymin = w.y_pixel(to_double(r.min().y())),
    xmax = w.x_pixel(to_double(r.max().x())),
    ymax = w.y_pixel(to_double(r.max().y()));
  w.get_painter().drawRect(xmin,ymin,xmax-xmin,ymax-ymin);
  w.do_paint();
  return w;
}
#endif // CGAL_ISO_RECTANGLE_2_H

#ifdef CGAL_BBOX_2_H
Qt_widget& operator<<(Qt_widget& w, const Bbox_2& r);
// see Qt_widget for the implementation of this non-template function
#endif // CGAL_BBOX_2_H




} // namespace CGAL

#endif // CGAL_QT_WIDGET_H
