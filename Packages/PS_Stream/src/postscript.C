#include <CGAL/Cartesian.h>
#include <iostream.h>
#include <fstream.h>

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/IO/Postscript_stream.h>


typedef CGAL::Point_2< CGAL::Cartesian<double> >     Point;
typedef CGAL::Segment_2< CGAL::Cartesian<double> >   Segment;
typedef CGAL::Ray_2< CGAL::Cartesian<double> >       Ray;
typedef CGAL::Triangle_2< CGAL::Cartesian<double> >  Triangle;
typedef CGAL::Iso_rectangle_2< CGAL::Cartesian<double> > Rect;
typedef CGAL::Circle_2< CGAL::Cartesian<double> >    Circle;
typedef CGAL::Bbox_2 BBox;


int  main()
{
  BBox bb(-2,-2,2,2);

  CGAL::PS_Stream ps(bb,15*CGAL::PS_Stream::CM,cout,CGAL::PS_Stream::QUIET_EPS);
  // CGAL::PS_Stream ps(bb,"toto.ps");
  CGAL::PS_Stream::Grid g(1,1);
  CGAL::PS_Stream::Axis a(0.5,0.5);
  Point p(-0.2,-0.7), q(0.6,0.3), r(0,0);
  Ray ra(q,p);
  Segment seg(q,r);
  Triangle tr(p,q,Point(-1,1));
  Rect re(p,q);
  Circle ci(r,1.0);
  CGAL::PS_Stream::Label l1("This is Standard Label (postscript)"), 
                        l2("Point p Label with another size");
  CGAL::PS_Stream::Latex_Label l3("This is a Latex Label {$\\alpha^{y}_{i}$}");
  CGAL::PS_Stream::Latex_Label l4("{Another Latex Label : $\\sum_{i=1}^{n} x_{i} = \\int_{0}^{1}f$}");
  CGAL::PS_Stream::Context c;
      
  ps << set_point_style(CGAL::PS_Stream::NONE) << p;
  ps << set_point_size(2) << set_line_width(0) << ra << ci;
  ps << set_fill(true) << set_fill_color(CGAL::Color(0,255,0));
  ps << set_point_style(CGAL::PS_Stream::EDOT) << r;
  ps << re <<  set_fill_color(CGAL::Color(0,0,255));
  ps << set_border_color(CGAL::Color(0,0,0)) << p << q;
  ps << move_to(r);
  c=ps.context();
  c.set_fill_color(CGAL::Color(255,255,0));
  ps << set_font("Helvetica-Oblique") << l1;
  ps << move_to(p);
  ps << set_font_size(8) << l2;
  ps << move_to(Point(-1.9,0)) << l3;
  ps << move_to(Point(-1.9,-1.9)) << l4;
  ps << set_current_context(CGAL::CTXT_DEFAULT);
  ps << set_line_style(CGAL::PS_Stream::DASH3) << tr;
  ps << Point(-1,-1) << Point(1,1) << Point(1,-1) << Point(-1,1)<<g<<a;
  ps << CGAL::PS_Stream::Border(1)<<seg;
  cout << flush;
  return 1;
}

