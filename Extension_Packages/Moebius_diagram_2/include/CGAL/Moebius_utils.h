#ifndef CHR_MOEBIUS_UTILS_H
#define CHR_MOEBIUS_UTILS_H

CGAL_BEGIN_NAMESPACE

#ifdef CHR_NOTRACE
#  define TRACE(output) do {} while (0)
#else
#  define TRACE(output) do (std::cout << output); while (0)
#endif

const int M_NOINTER = -1;
const int M_ABOVE = -2;
const int M_LEFT = -3;
const int M_RIGHT = -4;

CGAL_END_NAMESPACE

#ifdef CGAL_USE_QT


#include <CGAL/basic.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Point_2.h>


#include <CGAL/IO/Qt_widget.h>

#include <math.h>

CGAL_BEGIN_NAMESPACE

template <class K>
static double angle (double x, double y)
{
  static const double raddeg = 180 / 3.14159265358979;
#if 1
  double a = raddeg * atan2 (y, x);
  //if (a < 0) a = a + 360;
#else
  static const double threshold = std::sqrt (2) / 2;
  double norm = 1 / CGAL_NTS sqrt (x*x + y*y);
  double x_ = x * norm;
  double y_ = y * norm;
  
  double a = raddeg;

  if (x_ > threshold || x_ < - threshold) {
    a *= std::acos (x_);
    if (y_ < 0) a = 360 - a;
  } else {
    a *= std::asin (y_);
    if (x_ < 0) a = 360 - a;
  }
#endif
  TRACE ("      angle ("<<x<<", "<<y<<") = "<<a<<"\n");
  //  CGAL_postcondition (0 <= a && a <= 360);
  return a;
}

template <class K>
int my_round (double x)
{
  int i = (int) x;
  if (((double) i) - x > 0.5) return i - 1;
  if (((double) i) - x < -0.5) return i + 1;
  return i;
}

template <class K>
void qt_draw_arc (Qt_widget& w, const Circle_2<K>& c,
		  const Point_2<K> &source,
		  const Point_2<K> &target)
{
  static const double raddeg = 180 / 3.14159265358979;
  double cx = CGAL::to_double (c.center().x());
  double cy = CGAL::to_double (c.center().y());
  int rx = w.x_pixel_dist (std::sqrt (CGAL::to_double (c.squared_radius())));
  int ry = w.y_pixel_dist (std::sqrt (CGAL::to_double (c.squared_radius())));

  int icx = w.x_pixel (cx);
  int icy = w.y_pixel (cy);

  double x1 = CGAL::to_double (source.x()) - cx;
  double y1 = CGAL::to_double (source.y()) - cy;
  double x2 = CGAL::to_double (target.x()) - cx;
  double y2 = CGAL::to_double (target.y()) - cy;

  double da1 = raddeg * atan2 (y1, x1);
  double da2 = raddeg * atan2 (y2, x2);

  int a1 = my_round<K> (16.0 * da1);
  int a2 = my_round<K> (16.0 * da2);

  int da = a2 - a1;
  while (da < 0) da += 360 * 16;
  while (da > 360 * 16) da -= 360 * 16;
  CGAL_assertion (da >= 0);
  CGAL_assertion (da <= 360 * 16);
  int start = a1;

  if (c.orientation () == CGAL::CLOCKWISE) { 
    da = 360 * 16 - da;
    start = a2;
  }

  w.get_painter().drawArc (icx-rx, icy-ry, 2*rx, 2*ry, start, da); 
}

CGAL_END_NAMESPACE

#endif




#endif//CHR_MOEBIUS_UTILS_H
