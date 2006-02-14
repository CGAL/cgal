#include <CGAL/Cartesian.h>

#include <iostream>
#include <fstream>

#include <CGAL/IO/PS_Stream.h>

typedef CGAL::Cartesian<double>  K;
typedef K::Point_2               Point;
typedef K::Segment_2             Segment;
typedef K::Ray_2                 Ray;
typedef K::Triangle_2            Triangle;
typedef K::Iso_rectangle_2       Rect;
typedef K::Circle_2              Circle;
typedef CGAL::Bbox_2             BBox;

int main()
{
  CGAL::PS_Stream::PS_BBox bb(-2,-2,2,2);
  CGAL::PS_Stream ps(bb,300,"toto.ps",CGAL::PS_Stream::QUIET_EPS);

  CGAL::PS_Stream::Grid g(0.5,0.5);
  CGAL::PS_Stream::Axis a(1.0,1.0);
  Point p(-1,1), q(1,1), r(1,-1), s(-1,-1),som(0,2), b1(-0.5,1.5),
        b2(0.5,1.5), o(0,0);

  Triangle tr(b1,b2,som);
  Circle ci(o,1.0);
  CGAL::PS_Stream::Label l1("Centre");

  ps << p;
  ps << CGAL::line_width(2);    
  ps << CGAL::point_style(CGAL::PS_Stream::FDOT) << q;
  ps << CGAL::point_style(CGAL::PS_Stream::FBOX) << r;
  ps << CGAL::point_style(CGAL::PS_Stream::ICROSS) << s;
  ps << CGAL::fill(true) << CGAL::fill_color(CGAL::Color(0,255,0))
     << CGAL::border_color(CGAL::Color(0,0,255)) << ci;
  ps << CGAL::move_to(Point(0.15,0.1)) << CGAL::font("Helvetica-Oblique") << l1;
  ps << CGAL::fill(false) << CGAL::border_color(CGAL::Color(255,0,255)) << tr;
  ps << CGAL::current_context(CGAL::CTXT_DEFAULT);
  ps << CGAL::move_to(Point(0.6,1.5)) << CGAL::ps_label("triangle");
  ps << a;
  ps << CGAL::show_grid(g);
  ps << CGAL::border(3);

  return 0;
}
