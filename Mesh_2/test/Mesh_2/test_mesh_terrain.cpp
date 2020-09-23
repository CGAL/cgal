#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<Gt> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, TDS,
                                                   CGAL::Exact_predicates_tag> Delaunay;
typedef CGAL::Delaunay_mesh_size_criteria_2<Delaunay> Criteria;

typedef K::Point_3   Point;

int main()
{
  Delaunay dt;
  typedef Delaunay::Vertex_handle Vertex_handle;
  Vertex_handle va = dt.insert(Point(-4,0, 0));
  Vertex_handle vb = dt.insert(Point(0,-1, 0));
  Vertex_handle vc = dt.insert(Point(4,0, 0));
  Vertex_handle vd = dt.insert(Point(0,1,0));
  dt.insert(Point(2, 0.6, 0));

  dt.insert_constraint(va, vb);
  dt.insert_constraint(vb, vc);
  dt.insert_constraint(vc, vd);
  dt.insert_constraint(vd, va);

  CGAL::refine_Delaunay_mesh_2(dt, Criteria(0.125, 0.5));
  // dt.insert(begin, end);
  std::cout << dt.number_of_vertices() << std::endl;
  return 0;
}
