#include <iostream>
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/IO/PS_Stream.h>

typedef CGAL::Point_2< CGAL::Cartesian<double> >     Point;
typedef CGAL::Segment_2< CGAL::Cartesian<double> >   Segment;
typedef CGAL::Ray_2< CGAL::Cartesian<double> >       Ray;
typedef CGAL::Triangle_2< CGAL::Cartesian<double> >  Triangle;
typedef CGAL::Iso_rectangle_2< CGAL::Cartesian<double> > Rect;
typedef CGAL::Circle_2< CGAL::Cartesian<double> >    Circle;
typedef CGAL::Bbox_2 BBox;

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
  ps << point_style(CGAL::PS_Stream::FDOT) << q;
  ps << point_style(CGAL::PS_Stream::FBOX) << r;
  ps << point_style(CGAL::PS_Stream::ICROSS) << s;
  ps << CGAL::fill(true) << fill_color(CGAL::Color(0,255,0))
     << border_color(CGAL::Color(0,0,255)) << ci;
  ps << move_to(Point(0.15,0.1)) << CGAL::font("Helvetica-Oblique")<< l1;
  ps << CGAL::fill(false)<<border_color(CGAL::Color(255,0,255))<< tr;
  ps << current_context(CGAL::CTXT_DEFAULT);
  ps << move_to(Point(0.6,1.5)) <<CGAL::ps_label("triangle");
  ps <<a;
  ps <<show_grid(g);
  ps <<CGAL::border(3);

  return 0;
}
