#include <CGAL/Cartesian.h>

#include <iostream>
#include <fstream>

#include <CGAL/IO/PS_Stream.h>

typedef CGAL::Cartesian<double>   K;
typedef K::Point_2                Point;
typedef K::Segment_2              Segment;
typedef K::Ray_2                  Ray;
typedef K::Triangle_2             Triangle;
typedef K::Iso_rectangle_2        Rect;
typedef K::Circle_2               Circle;
typedef CGAL::Bbox_2              BBox;

int main()
{
 CGAL::PS_Stream::PS_BBox bb(-2,-2,2,2);

 //  CGAL::PS_Stream ps(bb,400,500,"titi.ps",CGAL::PS_Stream::QUIET_EPS);
 // CGAL::PS_Stream ps(bb,10,20,cout,CGAL::PS_Stream::QUIET);
  CGAL::PS_Stream ps(bb,300,"toto.ps",CGAL::PS_Stream::QUIET_EPS);
 // CGAL::PS_Stream ps(bb,10,cout,CGAL::PS_Stream::QUIET);
 // CGAL::PS_Stream ps(bb,"titi.ps",CGAL::PS_Stream::QUIET_EPS);
 // CGAL::PS_Stream ps(bb,cout,CGAL::PS_Stream::QUIET);


 CGAL::PS_Stream::Grid g(0.5,0.5);
 CGAL::PS_Stream::Axis a(1.0,1.0,2);
 Point p(-1,1), q(1,1), r(1,-1), s(-1,-1),som(0,2), b1(-0.5,1.5),
   b2(0.5,1.5), o(0,0);

 Triangle tr(b1,b2,som);
 // Rect re(p,q);
 Circle ci(o,1.0);
 CGAL::PS_Stream::Label l1("Centre");
 // CGAL::PS_Stream::Label l2("Point p Label with another size");
 // CGAL::PS_Stream::Latex_Label l3("This is a Latex Label {$\\alpha^{y}_{i}$}");
//  CGAL::PS_Stream::Latex_Label l4("{Another Latex Label :
// 		// 	 $\\sum_{i=1}^{n} x_{i} =
// 			/\\int_{0}^{1}f$}");	


 ps << p;
 ps << CGAL::line_width(2);    
 ps << point_style(CGAL::PS_Stream::FDOT) << q;
 ps << point_style(CGAL::PS_Stream::FBOX) << r;
 ps << point_style(CGAL::PS_Stream::ICROSS) << s;
 ps << point_style(CGAL::PS_Stream::FDOT) << CGAL::point_size(2) << o;
 ps << CGAL::fill(true) <<  fill_color(CGAL::Color(0,255,0))<<border_color(CGAL::Color(0,0,255))<<ci;
//  ps << CGAL::show_direction(true)<<seg;
 // ps << set_point_style(CGAL::PS_Stream::FDOT) << CGAL:: set_point_size(3)<<set_border_color(CGAL::Color(0,0,255))<< r;
 // ps <<set_fill_color(CGAL::Color(0,0,255))<<re;
 // ps << set_border_color(CGAL::Color(255,0,0)) << p << q;
 ps << move_to(Point(0.15,0.1)) << CGAL:: font("Helvetica-Oblique") << l1;
 // ps <<CGAL:: set_font("Helvetica-Oblique") << l1;
 // ps << CGAL::put_latex_label("latex_label");
 ps << CGAL::fill(true)<<border_color(CGAL::Color(255,0,255))<< tr;

// ps << set_font_size(8) << l2; 
 // ps <<  CGAL:: set_font_size(8)<< l3 ;
 ps << current_context(CGAL::CTXT_DEFAULT);
 //ps.set_default_context();
 ps << move_to(Point(0.6,1.5)) <<CGAL::latex_label("triangle");
//  ps << q;
 // ps << Point(-1,-1) << Point(1,1) << Point(1,-1) << Point(-1,1);
//  ps << CGAL::PS_Stream::Border(1)<<seg;
//  ps << CGAL::put_border(3);
//  ps << CGAL::set_font("Helvetica-Oblique");
//  ps << CGAL::set_font_size(8);
 ps << a;
 ps <<show_grid(g);
 ps <<CGAL::border(3);
// ps <<g;
 // ps <<set_border_color(CGAL::Color(255,170,0))<<l3;
 
 return 0;
}
