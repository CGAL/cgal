// This file tests the deprecated headers
// Triangulation_euclidean_traits_.._3 on purpose.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iostream>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Projection_traits_xy_3<Kernel>   Gt_xy;
typedef CGAL::Projection_traits_yz_3<Kernel>   Gt_yz;
typedef CGAL::Projection_traits_xz_3<Kernel>   Gt_xz;
typedef CGAL::Delaunay_triangulation_2<Gt_xy>  DT_xy;
typedef CGAL::Delaunay_triangulation_2<Gt_xz>  DT_xz;
typedef CGAL::Delaunay_triangulation_2<Gt_yz>  DT_yz;
typedef Gt_yz::Point_2                         Point;

typedef CGAL::Triangulation_data_structure_2 <
                       CGAL::Triangulation_vertex_base_2<Gt_xy>,
		       CGAL::Constrained_triangulation_face_base_2<Gt_xy> > CDT_TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt_xy,CDT_TDS,CGAL::Exact_predicates_tag> CDT;

int main(){
  DT_xy t_xy;
  DT_yz t_yz;
  DT_xz t_xz;
  
  Point pts[3]={Point(2,3,5),Point(5,1,4),Point(4,2,12)};
  
  t_xy.insert (pts,pts+3); t_xy.insert(Point(4,3,5)); t_xy.insert(Point(4,7,12)); 
  t_yz.insert (pts,pts+3); t_yz.insert(Point(4,3,5)); 
  t_xz.insert (pts,pts+3); t_xz.insert(Point(4,7,12)); 

  
  assert( t_xy.number_of_vertices()==5 );
  assert( t_yz.number_of_vertices()==3 );
  assert( t_xz.number_of_vertices()==3 );
  
  CDT cdt;

  cdt.insert(pts,pts+3);
  cdt.insert_constraint(Point(1,1,4),Point(0,0,4));
  cdt.insert_constraint(Point(0,1,2),Point(1,0,2));
  
  assert (cdt.number_of_vertices() == 8 );
  return 0;
}
