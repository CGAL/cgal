#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <cassert>
#include <fstream>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel               Gt;
typedef CGAL::Tag_true                                                    Alpha_cmp_tag;
//We use CGAL::Default to skip the complete declaration of base classes
typedef CGAL::Alpha_shape_vertex_base_3<Gt,CGAL::Default,Alpha_cmp_tag>   Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt,CGAL::Default,Alpha_cmp_tag>     Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>                       Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>                            Triangulation_3;
//Alpha shape with ExactAlphaComparisonTag set to true (note that the tag is also
//set to true for Vb and Fb)
typedef CGAL::Alpha_shape_3<Triangulation_3,Alpha_cmp_tag>                Alpha_shape_3;
typedef Gt::Point_3                                                       Point;

int main()
{
  //Set of points for which the alpha shapes cannot be computed with
  //a floating point alpha value (on certain platforms)
  std::list<Point> lp;
  lp.push_back(Point(49.2559,29.1767,6.7723));
  lp.push_back(Point(49.3696,31.4775,5.33777));
  lp.push_back(Point(49.4264,32.6279,4.6205));
  lp.push_back(Point(49.3127,30.3271,6.05503));

   // compute alpha shape
  Alpha_shape_3 as(lp.begin(),lp.end(),0,Alpha_shape_3::GENERAL);

  return 0;
}
