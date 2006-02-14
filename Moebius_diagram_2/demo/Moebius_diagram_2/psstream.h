#ifndef CHR_PSSTREAM_H
#define CHR_PSSTREAM_H

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>

#include <CGAL/Point_2.h>

class psstream : public std::ofstream
{
public:
  typedef enum { CROSS } POINT_STYLE;

  psstream (const char *filename,
	    double xmin, double ymin,
	    double xmax, double ymax) :
    _style (CROSS),
    _xmin (xmin), _ymin (ymin),
    _xmax (xmax), _ymax (ymax),
    _scale (),
    _os (filename)
    {
      assert (xmin < xmax);
      assert (ymin < ymax);

      _os << "%!PS-Adobe-2.0 EPSF-2.0\n"
	  << "%%BoundingBox: 0 0 300 200\n"
	  << "%%EndComments\n";

      double sx = 300 / (xmax-xmin);
      double sy = 200 / (ymax-ymin);
      _scale = (sx < sy) ? sx : sy;
      _os << _scale << " " << _scale << " scale\n";

      _os << -xmin << " " << -xmax << " translate\n";
    }

  double xmin () const { return _xmin; }
  double ymin () const { return _ymin; }
  double xmax () const { return _xmax; }
  double ymax () const { return _ymax; }
  double scale () const { return _scale; }
  POINT_STYLE point_style () const { return _style; }
  //  std::ostream &stream const () { return *this; }

  psstream &operator<< (psstream::POINT_STYLE style)
    {
      _style = style;
      return *this;
    }

  template <class K>
  psstream &operator<< (CGAL::Point_2<K> p) {
    double x = CGAL::to_double (p.x());
    double y = CGAL::to_double (p.y());
    if (x < xmin() || x > xmax() || y < ymin() || y > ymax()) return *this;
    
    switch (_style) {
    case CROSS: {
      double s = ((xmax()-xmin()) + (ymax()-ymin())) / 1000.0;
      _os << "newpath\n  "
	  << x - s << " " << y - s << " moveto "
	  << 2 * s << " " << 2 * s << " rlineto "
	  << x - s << " " << y + s << " moveto "
	  << 2 * s << " " << - 2 * s << " rlineto\nclosepath stroke\n";
    } break;
    }
    return *this;
  }

  template <class K>
  void draw_arc (const CGAL::Circle_2<K> &c,
		 const CGAL::Point_2<K> &source,
		 const CGAL::Point_2<K> &target) {
    static const double raddeg = 180 / 3.14159265358979;
    double cx = CGAL::to_double (c.center().x());
    double cy = CGAL::to_double (c.center().y());
    double r2 = CGAL::to_double (c.squared_radius ());
    double r = std::sqrt (r2);

    if (cx + r < xmin () ||
	cx - r > xmax () ||
	cy + r < ymin () ||
	cy - r > ymax ())
      return;
    
    double xmin_c = xmin () - cx;
    double ymin_c = ymin () - cy;
    double xmax_c = xmax () - cx;
    double ymax_c = ymax () - cy;

    if (xmin_c > 0) {
      if (ymin_c > 0)
	if (xmin_c*xmin_c + ymin_c*ymin_c > r2) return;
      else if (ymax_c < 0)
	if (xmin_c*xmin_c + ymax_c*ymax_c > r2) return;
    } else if (xmax_c < 0) {
      if (ymin_c > 0)
	if (xmax_c*xmax_c + ymin_c*ymin_c > r2) return;
      else if (ymax_c < 0)
	if (xmax_c*xmax_c + ymax_c*ymax_c > r2) return;
    }

    double x1 = CGAL::to_double (source.x()) - cx;
    double y1 = CGAL::to_double (source.y()) - cy;
    double x2 = CGAL::to_double (target.x()) - cx;
    double y2 = CGAL::to_double (target.y()) - cy;

    double a1, a2;
    if (c.orientation() == CGAL::COUNTERCLOCKWISE) {
      a1 = raddeg * atan2 (y1, x1);
      a2 = raddeg * atan2 (y2, x2);
    } else {
      a2 = raddeg * atan2 (y1, x1);
      a1 = raddeg * atan2 (y2, x2);
    }

    _os << cx << " " << cy << " " << r << "  "
	<< a1 << " " << a2 << " arc stroke\n";
  }

  ~psstream () {
    _os << "showpage\n";
  }

private:
  POINT_STYLE _style;
  double _xmin, _ymin, _xmax, _ymax;
  double _scale;
  std::ofstream _os;
};

#endif//CHR_PSSTREAM_H
