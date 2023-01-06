#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Projection_traits_3.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K1;
typedef CGAL::Projection_traits_3<K1> K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

int main()
{
  K gt{ { 0, 0, 1} };
  CDT cdt(gt);

  Vertex_handle va = cdt.insert(Point(-4,0,0));
  Vertex_handle vb = cdt.insert(Point(0,-1,0));
  Vertex_handle vc = cdt.insert(Point(4,0,0));
  Vertex_handle vd = cdt.insert(Point(0,1,0));
  cdt.insert(Point(2, 0.6, 0));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::cout << "Meshing the triangulation..." << std::endl;
  Criteria criteria(0.125, 0.5, gt);
  CGAL::refine_Delaunay_mesh_2(cdt, CGAL::parameters::criteria(criteria));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
}
