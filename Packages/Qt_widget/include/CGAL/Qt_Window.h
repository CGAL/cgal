#ifndef CGAL_QT_WINDOW_H
#define CGAL_QT_WINDOW_H

#include <qwidget.h>
#include <qpainter.h>
#include <qcolor.h>
#include <qpixmap.h>

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>

CGAL_BEGIN_NAMESPACE
enum PointStyle { PIXEL, CROSS, PLUS, CIRCLE, DISC, RECT, BOX };
CGAL_END_NAMESPACE

//CGAL_BEGIN_NAMESPACE

class QCGALWidget : public QWidget {
  Q_OBJECT
    Q_PROPERTY( QColor color READ color WRITE setColor )
    Q_PROPERTY( QColor backGroundColor READ backgroundColor WRITE
 		setBackgroundColor )
    Q_PROPERTY( QColor fillColor READ fillColor WRITE
 		setFillColor )
    Q_PROPERTY( bool isFilled READ isFilled WRITE setFilled )
    Q_PROPERTY( uint lineWidth READ lineWidth WRITE setLineWidth )
    Q_PROPERTY( uint pointSize READ pointSize WRITE setPointSize )
/*     Q_PROPERTY( PointStyle pointStyle READ pointStyle WRITE */
/* 		setPointStyle ) */
/*     Q_ENUMS( PointStyle ) */
public:
  // constructor
  QCGALWidget(QWidget *parent = 0, const char *name = 0);
  // destructor
  ~QCGALWidget() {};

  // initialization of coordinates system
  void init(double x_min, double x_max, double y_min, double y_max);
  void init(double x_min, double x_max, double y_min);
  bool isInitialized() const; // tell if init has been called

  // color types convertors
  static QColor CGAL2Qt_Color(CGAL::Color c);
  static CGAL::Color Qt2CGAL_color(QColor c);

  // painting system
  inline QPainter& painter() { return paint; };
  void lock() { ++Locked; };
  void unlock() { if (Locked>0) --Locked; doPaint(); };
  void doPaint() { if (Locked==0) repaint(FALSE); };

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

  // CGAL version of setFooColor
  // used by the manipulators system
  inline void setColor(CGAL::Color c)
    { setColor(CGAL2Qt_Color(c)); };
  inline void setBackgroundColor(CGAL::Color c)
    { setBackgroundColor(CGAL2Qt_Color(c)); };
  inline void setFillColor(CGAL::Color c)
    { setFillColor(CGAL2Qt_Color(c)); };

  // set pen() color to c, cf. manipulators below for setting
  // backgroundColor and fillColor
  QCGALWidget& operator<<(const CGAL::Color& c);
  // set point style
  QCGALWidget& operator<<(const CGAL::PointStyle& ps);
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

signals:
  void mousePressed(QMouseEvent *e);
  void mouseReleased(QMouseEvent *e);
  void mouseMoved(QMouseEvent *e);
  void redraw();

protected:
  void paintEvent(QPaintEvent *e);
  void resizeEvent(QResizeEvent *e);
  void mousePressEvent(QMouseEvent *e);
  void mouseReleaseEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);

private:
  void initialize(); // initialize initiale dimensions
  bool initialized;
  void setScales(); // set xscal and yscal

  unsigned int Locked;
  // point style and size
  uint _pointSize;
  PointStyle _pointStyle;

  QPixmap pixmap; // the pixmap on which paints the painter
  QPainter paint; // the painter
  QBrush savedBrush; // saved brush, to be able to restore it on
  // setFilled(true)

  double xmin, xmax, ymin, ymax; // real dimensions
  double xscal, yscal; // scalings int/double
};

typedef QCGALWidget        Window_stream;

inline
bool QCGALWidget:: isInitialized() const
{
  return initialized;
}

// manipulators
// ~~~~~~~~~~~~
// single manipulators
inline
QCGALWidget& operator<<(QCGALWidget& w, QCGALWidget& (*m)(QCGALWidget&))
{
  return m(w);
};

// w << noFill << ... stop the filling of geometrical object
inline
QCGALWidget& noFill(QCGALWidget& w)
{
  w.setFilled(false);
  return w;
}

// manipulators with one argument
template <class Param>
struct QCGALWidgetManip {
  QCGALWidget& (*f)(QCGALWidget&, Param);
  Param p;
  QCGALWidgetManip(QCGALWidget& (*ff)(QCGALWidget&, Param),
		 Param pp) : f(ff), p(pp) {}
};

// usage: w << manip(Param) << ...
template <class Param>
QCGALWidget& operator<<(QCGALWidget& w, QCGALWidgetManip<Param> m)
{
  return m.f(w, m.p);
}

#define __DEFINE_QCGALWIDGET_MANIP(param,function) \
inline \
QCGALWidget& __QCGALWidgetManip##function##Aux (QCGALWidget& w, param p) \
{ w.set##function(p); return w; } \
inline \
QCGALWidgetManip<param> function(param p) \
{ return QCGALWidgetManip<param>( __QCGALWidgetManip##function##Aux, p); }

// w << BackgroundColor(c) << ... sets the background color
__DEFINE_QCGALWIDGET_MANIP( CGAL::Color, BackgroundColor )

// w << FillColor(c) << ... sets the fill color
__DEFINE_QCGALWIDGET_MANIP( CGAL::Color, FillColor )

// w << LineWidth(i) << ... sets lines width
__DEFINE_QCGALWIDGET_MANIP( uint, LineWidth )

// w << PointSize(i) << ... sets points size
__DEFINE_QCGALWIDGET_MANIP( uint, PointSize )

// color types convertors
// ~~~~~~~~~~~~~~~~~~~~~~
inline
QColor QCGALWidget::CGAL2Qt_Color(CGAL::Color c)
{
  return QColor(c.red(), c.green(), c.blue());
}

inline
CGAL::Color QCGALWidget::Qt2CGAL_color(QColor c)
{
  return CGAL::Color(c.red(),c.green(),c.blue());
}

// properties
// ~~~~~~~~~~
inline
QColor QCGALWidget::color() const
{
  return paint.pen().color();
};


inline
void QCGALWidget::setColor(QColor c)
{
  QPen p=painter().pen();
  p.setColor(c);
  painter().setPen(p);
}

inline
QColor QCGALWidget::backgroundColor() const
{
  return paint.backgroundColor();
}

inline
void QCGALWidget::setBackgroundColor(QColor c)
{
  QWidget::setBackgroundColor(c);
  painter().setBackgroundColor(c);
  clear();
}

inline
QColor QCGALWidget::fillColor() const
{
  return paint.brush().color();
}

inline
void QCGALWidget::setFillColor(QColor c)
{
  setFilled(true);
  painter().setBrush(c);
}

inline
bool QCGALWidget::isFilled() const
{
  return( paint.brush().style()==Qt::NoBrush );
}

inline
void QCGALWidget::setFilled(bool f)
{
  if (f)
    paint.setBrush(savedBrush);
  else
    {
      savedBrush=paint.brush();
      paint.setBrush(QBrush());
    };
}

inline
uint QCGALWidget::lineWidth() const
{
  return( paint.pen().width());
}

inline
void QCGALWidget::setLineWidth(uint i)
{
  QPen p=painter().pen();
  p.setWidth(i);
  painter().setPen(p);
}

inline
uint QCGALWidget::pointSize() const
{
  return _pointSize;
}

inline
void QCGALWidget::setPointSize(uint i)
{
  _pointSize=i;
}

inline
QCGALWidget::PointStyle QCGALWidget::pointStyle() const
{
  return _pointStyle;
}

inline
void QCGALWidget::setPointStyle(PointStyle ps)
{
  _pointStyle=ps;
}

// drawing methods
// ~~~~~~~~~~~~~~~

template <class R>
QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Point_2<R>& p)
{
  int
    x=w.x_pixel(CGAL::to_double(p.x())),
    y=w.y_pixel(CGAL::to_double(p.y()));

  uint size=w.pointSize();
  CGAL::PointStyle ps=w.pointStyle();

  switch (ps)
    {
    case CGAL::PIXEL:
      {
	w.painter().drawPoint(x,y);
	break;
      }
    case CGAL::CROSS:
      {
	w.painter().drawLine(x-size/2, y-size/2, x+size/2, y+size/2);
	w.painter().drawLine(x-size/2, y+size/2, x+size/2, y-size/2);
	break;
      }
    case CGAL::PLUS:
      {
	w.painter().drawLine(x, y-size/2, x, y+size/2);
	w.painter().drawLine(x-size/2, y, x+size/2, y);
	break;
      }
    case CGAL::CIRCLE:
      {
	QBrush old_brush=w.painter().brush();
	w.painter().setBrush(QBrush());
	w.painter().drawEllipse(x-size/2, y-size/2, size, size);
	w.painter().setBrush(old_brush);
	break;
      }
    case CGAL::DISC:
      {
	QBrush old_brush=w.painter().brush();
	w.painter().setBrush(w.painter().pen().color());
	w.painter().drawEllipse(x-size/2, y-size/2, size, size);
	w.painter().setBrush(old_brush);
	break;
      }
    case CGAL::RECT:
      {
	QBrush old_brush=w.painter().brush();
	w.painter().setBrush(QBrush());
	w.painter().drawRect(x-size/2, y-size/2, size, size);
	w.painter().setBrush(old_brush);
	break;
      }
    case CGAL::BOX:
      {
	QBrush old_brush=w.painter().brush();
	w.painter().setBrush(w.painter().pen().color());
	w.painter().drawRect(x-size/2, y-size/2, size, size);
	w.painter().setBrush(old_brush);
	break;
      }
    };
  w.doPaint();
  return w;
}

#ifdef CGAL_SEGMENT_2_H
template <class R>
QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Segment_2<R>& s)
{
  const int
    x1=w.x_pixel(CGAL::to_double(s.source().x())),
    y1=w.y_pixel(CGAL::to_double(s.source().y())),
    x2=w.x_pixel(CGAL::to_double(s.target().x())),
    y2=w.y_pixel(CGAL::to_double(s.target().y()));
  w.painter().drawLine(x1,y1,x2,y2);
  w.doPaint();
  return w;
}
#endif // CGAL_SEGMENT_2_H

#ifdef CGAL_LINE_2_H

template <class R>
QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Line_2<R>& l)
{
  typedef CGAL::Cartesian<double> Rep;
  typedef CGAL::Point_2<Rep> Point;

  const CGAL::Point_2<R>
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

  w.painter().drawLine(w.x_pixel(x1),w.y_pixel(y1),
		       w.x_pixel(x2),w.y_pixel(y2));
  return w;
}

#endif // CGAL_LINE_2_H

#ifdef CGAL_RAY_2_H
template <class R>
QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Ray_2<R>& r)
{
  typedef CGAL::Cartesian<double> Rep;
  typedef CGAL::Point_2<Rep> Point;

  const CGAL::Point_2<R>
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
  w.painter().drawLine(w.x_pixel(p1d.x()),w.y_pixel(p1d.y()),
		       w.x_pixel(x),w.y_pixel(y));
  return w;

}
#endif //CGAL_RAY_2_H

#ifdef CGAL_TRIANGLE_2_H
template< class R >
QCGALWidget&
operator<<(QCGALWidget& w, const CGAL::Triangle_2<R>& t)
{
  const int
    ax = w.x_pixel(CGAL::to_double(t.vertex(0).x())),
    ay = w.y_pixel(CGAL::to_double(t.vertex(0).y())),
    bx = w.x_pixel(CGAL::to_double(t.vertex(1).x())),
    by = w.y_pixel(CGAL::to_double(t.vertex(1).y())),
    cx = w.x_pixel(CGAL::to_double(t.vertex(2).x())),
    cy = w.y_pixel(CGAL::to_double(t.vertex(2).y()));

  QPointArray array;

  array.setPoints(3,ax,ay,bx,by,cx,cy);
  w.painter().drawPolygon(array);
  w.doPaint();
  return w;
}
#endif

#ifdef CGAL_CIRCLE_2_H
template < class R>
QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Circle_2<R>& c)
{
  int 
    cx=w.x_pixel(CGAL::to_double(c.center().x())),
    cy=w.y_pixel(CGAL::to_double(c.center().y())),
    rx=w.x_pixel_dist((std::sqrt(CGAL::to_double(c.squared_radius())))),
    ry=w.y_pixel_dist((std::sqrt(CGAL::to_double(c.squared_radius()))));

  w.painter().drawEllipse(cx-rx,cy-ry,2*rx,2*ry);
  w.doPaint();
  return w;
}
#endif // CGAL_CIRCLE_2_H

#ifdef CGAL_ISO_RECTANGLE_2_H
template< class R >
QCGALWidget&
operator<<(QCGALWidget& w, const CGAL::Iso_rectangle_2<R>& r)
{
  int
    xmin = w.x_pixel(CGAL::to_double(r.min().x())),
    ymin = w.y_pixel(CGAL::to_double(r.min().y())),
    xmax = w.x_pixel(CGAL::to_double(r.max().x())),
    ymax = w.y_pixel(CGAL::to_double(r.max().y()));

  w.painter().drawRect(xmin,ymin,xmax-xmin,ymax-ymin);
  w.doPaint();
  return w;
}
#endif // CGAL_ISO_RECTANGLE_2_H

#ifdef CGAL_BBOX_2_H
QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Bbox_2& r);
// see Qt_Window for the implementation of this non-template function
#endif // CGAL_BBOX_2_H

#ifdef CGAL_POLYGON_2_H
template <class Tr,class Co>
QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Polygon_2<Tr,Co>& pol)
{
  typedef CGAL::Polygon_2<Tr,Co>::Vertex_const_iterator VI;
  QPointArray array;

  array.resize(pol.size());

  unsigned int n=0;
  for(VI i=pol.vertices_begin();i!=pol.vertices_end();i++)
    {
      array.setPoint(n++,w.x_pixel(CGAL::to_double(i->x())),
		     w.y_pixel(CGAL::to_double(i->y())));
    }
  w.painter().drawPolygon(array);
  w.doPaint();
  return w;
}
#endif // CGAL_POLYGON_2_H

#ifdef CGAL_TRIANGULATION_2_H
template < class Gt, class Tds>
QCGALWidget&
operator<<(QCGALWidget& w,  const CGAL::Triangulation_2<Gt, Tds> &t)
{
  w.lock();
  t.draw_triangulation(w);
  w.unlock();
  return w;
}
#endif // CGAL_TRIANGULATION_2_H

#ifdef CGAL_DELAUNAY_TRIANGULATION_2_H
template < class Gt, class Tds >
QCGALWidget&
operator<<(QCGALWidget& w,  const CGAL::Delaunay_triangulation_2<Gt,Tds> &dt)
{
  w.lock();
  dt.draw_triangulation(w);
  w.unlock();
  w.doPaint();
  return w;
}
#endif // CGAL_DELAUNAY_TRIANGULATION_2_H

#ifdef CGAL_CONSTRAINED_TRIANGULATION_2_H
template < class Gt, class Tds>
QCGALWidget&
operator<<(QCGALWidget& w,  const CGAL::Constrained_triangulation_2<Gt,Tds> &t)
{
  w.lock();
  t.draw_triangulation(w);
  w.unlock();
  return w;
}
#endif // CGAL_CONSTRAINED_TRIANGULATION_2_H

#ifdef CGAL_REGULAR_TRIANGULATION_2_H
template < class Gt, class Tds >
QCGALWidget&
operator<<(QCGALWidget& w, CGAL::Regular_triangulation_2<Gt,Tds> &t)
{
  w.lock();
  t.draw_triangulation(w);
  w.unlock();
  return w;
}
#endif // CGAL_REGULAR_TRIANGULATION_2_H


//CGAL_END_NAMESPACE

#endif // QT_WINDOW_H
