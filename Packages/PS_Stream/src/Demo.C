#include <iostream.h>
#include <fstream.h>

#include <CGAL/cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/IO/Postscript_stream.h>

typedef CGAL_Point_2< C<double> >     Point;
typedef CGAL_Segment_2< C<double> >   Segment;
typedef CGAL_Ray_2< C<double> >       Ray;
typedef CGAL_Triangle_2< C<double> >  Triangle;
typedef CGAL_Iso_rectangle_2< C<double> > Rect;
typedef CGAL_Circle_2< C<double> >    Circle;
typedef CGAL_Bbox_2 BBox;


int  main(int argc, char* argv[])
{
  BBox bb(-0.5,0,10,9);

  CGAL_PS_stream ps(bb,8*CGAL_PS_stream::CM,argv[1],CGAL_PS_stream::QUIET_EPS);
  Point p[7];
  Point q[7];
  Point r[7];
  Segment s[7];
  for (int i=1;i<8;i++)
  {
    p[i-1]=Point(0,i+0.1);
    q[i-1]=Point(1,i);
    r[i-1]=Point(8,i);
    s[i-1]=Segment(Point(5,i), Point(7,i));
  }
  
  ps << set_point_style(CGAL_PS_Context::NONE) << p[0];
  ps << set_point_style(CGAL_PS_Context::XCROSS) << p[1];
  ps << set_point_style(CGAL_PS_Context::ICROSS) << p[2];
  ps << set_point_style(CGAL_PS_Context::EDOT) << p[3];
  ps << set_point_style(CGAL_PS_Context::FDOT) << p[4];
  ps << set_point_style(CGAL_PS_Context::EBOX) << p[5];
  ps << set_point_style(CGAL_PS_Context::FBOX) << p[6];
  ps << move_to(q[0]) << put_ps_label("NONE");
  ps << move_to(q[1]) << put_ps_label("XCROSS");
  ps << move_to(q[2]) << put_ps_label("ICROSS");
  ps << move_to(q[3]) << put_ps_label("EDOT");
  ps << move_to(q[4]) << put_ps_label("FDOT");
  ps << move_to(q[5]) << put_ps_label("EBOX");
  ps << move_to(q[6]) << put_ps_label("FBOX");
  ps << set_line_style(CGAL_PS_Context::SOLID) << s[0];
  ps << set_line_style(CGAL_PS_Context::DASH1) << s[1];
  ps << set_line_style(CGAL_PS_Context::DASH2) << s[2];
  ps << set_line_style(CGAL_PS_Context::DASH3) << s[3];
  ps << set_line_style(CGAL_PS_Context::DASH4) << s[4];
  ps << set_line_style(CGAL_PS_Context::DASH5) << s[5];
  ps << set_line_style(CGAL_PS_Context::DASH6) << s[6];
  ps << move_to(r[0]) << put_ps_label("SOLID");
  ps << move_to(r[1]) << put_ps_label("DASH1");
  ps << move_to(r[2]) << put_ps_label("DASH2");
  ps << move_to(r[3]) << put_ps_label("DASH3");
  ps << move_to(r[4]) << put_ps_label("DASH4");
  ps << move_to(r[5]) << put_ps_label("DASH5");
  ps << move_to(r[6]) << put_ps_label("DASH6");
}
