#include "Qt_Window.h"
#include <CGAL/Bbox_2.h>
#include <CGAL/IO/esprit_logo.xpm>

//CGAL_BEGIN_NAMESPACE

QCGALWidget::QCGALWidget(QWidget *parent, const char *name) :
  QWidget(parent, name), initialized(false), Locked(), _pointSize(4),
  _pointStyle(CGAL::DISC)
{
  setCaption("Qt CGAL::Window_stream");
  setIcon(QPixmap(CGAL::esprit_logo));
  initialize();
  paint.begin(&pixmap);
  setBackgroundColor(Qt::white);
  paint.setPen(QPen(Qt::black,3));
  clear();
}

void QCGALWidget::initialize()
{
  xmin=0;
  xmax=width()-1;
  ymin=0;
  ymax=height()-1;
  xscal=1;
  yscal=1;
  pixmap.resize(size());
}

void QCGALWidget::setScales()
{
  xscal=(width()-1)/(xmax-xmin);
  yscal=(height()-1)/(ymax-ymin);
}

void QCGALWidget::resizeEvent(QResizeEvent *e)
{
  // save paint state
  QFont f=paint.font();
  QBrush b=paint.brush();
  QPen p=paint.pen();
  QColor bc=paint.backgroundColor();

  paint.end();  // end painting on pixmap

  /*
    the only difference between an initialized QCGALWidget and a
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
  emit(redraw());
}

void QCGALWidget::paintEvent(QPaintEvent *e)
{
  // save paint state
  QFont f=paint.font();
  QBrush b=paint.brush();
  QPen p=paint.pen();
  QColor bc=paint.backgroundColor();

  paint.end();  // end painting on pixmap
  bitBlt(this, 0, 0, &pixmap); // copy pixmap to the QCGALWidget
  paint.begin(&pixmap); // begin again painting on pixmap

  // restore paint state
  paint.setFont(f);
  paint.setBrush(b);
  paint.setPen(p);
  paint.setBackgroundColor(bc);
}

void QCGALWidget::mousePressEvent(QMouseEvent *e)
{
  emit(mousePressed(e));
}

void QCGALWidget::mouseReleaseEvent(QMouseEvent *e)
{
  emit(mouseReleased(e));
}

void QCGALWidget::mouseMoveEvent(QMouseEvent *e)
{
  emit(mouseMoved(e));
}

void QCGALWidget::init(double x_min, double x_max, double y_min)
{
  double y_max=y_min+(height()-1)*(x_max-x_min)/(width()-1);
  init(x_min, x_max, y_min, y_max);
}

void QCGALWidget::init(double  x_min, double x_max, double y_min, double y_max)
{
  xmin=x_min;
  xmax=x_max;
  ymin=y_min;
  ymax=y_max;
  setScales();
  initialized=true;
}

double QCGALWidget::x_real(int x) const
{
  return(xmin+x/xscal);
}

double QCGALWidget::y_real(int y) const
{
  return(ymax-y/yscal);
}

double QCGALWidget::x_real_dist(double d) const
{
  return(d/xscal);
}

double QCGALWidget::y_real_dist(double d) const
{
  return(d/yscal);
}

int QCGALWidget::x_pixel(double x) const
{
  return( static_cast<int>((x-xmin)*xscal+0.5) );
}

int QCGALWidget::y_pixel(double y) const
{
  return( - static_cast<int>((y-ymax)*yscal+0.5) );
}

int QCGALWidget::x_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*xscal+0.5) );
  else
    return( static_cast<int>(d*xscal-0.5) );
}

int QCGALWidget::y_pixel_dist(double d) const
{
  if (d>0)
    return( static_cast<int>(d*yscal+0.5) );
  else
    return( static_cast<int>(d*yscal-0.5) );
}

QCGALWidget& QCGALWidget::operator<<(const CGAL::Color& c)
{
  setColor(CGAL2Qt_Color(c));
  return *this;
}

QCGALWidget& QCGALWidget::operator<<(const CGAL::PointStyle& ps)
{
  setPointStyle(ps);
  return *this;
}

void QCGALWidget::clear() {
  painter().eraseRect(rect());
}

QCGALWidget& operator<<(QCGALWidget& w, const CGAL::Bbox_2& r)
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

#include "Qt_Window.moc"

//CGAL_END_NAMESPACE
