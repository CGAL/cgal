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


typedef CGAL_Point_2< CGAL_Cartesian<double> >     Point;
typedef CGAL_Segment_2< CGAL_Cartesian<double> >   Segment;
typedef CGAL_Ray_2< CGAL_Cartesian<double> >       Ray;
typedef CGAL_Triangle_2< CGAL_Cartesian<double> >  Triangle;
typedef CGAL_Iso_rectangle_2< CGAL_Cartesian<double> > Rect;
typedef CGAL_Circle_2< CGAL_Cartesian<double> >    Circle;
typedef CGAL_Bbox_2 BBox;


int  main()
{
  BBox bb(-2,-2,2,2);

  CGAL_PS_Stream ps(bb,15*CGAL_PS_Stream::CM,cout,CGAL_PS_Stream::QUIET_EPS);
  // CGAL_PS_Stream ps(bb,"toto.ps");
  CGAL_PS_Stream::Grid g(1,1);
  CGAL_PS_Stream::Axis a(0.5,0.5);
  Point p(-0.2,-0.7), q(0.6,0.3), r(0,0);
  Ray ra(q,p);
  Segment seg(q,r);
  Triangle tr(p,q,Point(-1,1));
  Rect re(p,q);
  Circle ci(r,1.0);
  CGAL_PS_Stream::Label l1("This is Standard Label (postscript)"), 
                        l2("Point p Label with another size");
  CGAL_PS_Stream::Latex_Label l3("This is a Latex Label {$\\alpha^{y}_{i}$}");
  CGAL_PS_Stream::Latex_Label l4("{Another Latex Label : $\\sum_{i=1}^{n} x_{i} = \\int_{0}^{1}f$}");
  CGAL_PS_Stream::Context c;
      
  ps << set_point_style(CGAL_PS_Stream::NONE) << p;
  ps << set_point_size(2) << set_line_width(0) << ra << ci;
  ps << set_fill(true) << set_fill_color(CGAL_Color(0,255,0));
  ps << set_point_style(CGAL_PS_Stream::EDOT) << r;
  ps << re <<  set_fill_color(CGAL_Color(0,0,255));
  ps << set_border_color(CGAL_Color(0,0,0)) << p << q;
  ps << move_to(r);
  c=ps.context();
  c.set_fill_color(CGAL_Color(255,255,0));
  ps << set_font("Helvetica-Oblique") << l1;
  ps << move_to(p);
  ps << set_font_size(8) << l2;
  ps << move_to(Point(-1.9,0)) << l3;
  ps << move_to(Point(-1.9,-1.9)) << l4;
  ps << set_current_context(CGAL_CTXT_DEFAULT);
  ps << set_line_style(CGAL_PS_Stream::DASH3) << tr;
  ps << Point(-1,-1) << Point(1,1) << Point(1,-1) << Point(-1,1)<<g<<a;
  ps << CGAL_PS_Stream::Border(1)<<seg;
  cout << flush;
  return 1;
}

