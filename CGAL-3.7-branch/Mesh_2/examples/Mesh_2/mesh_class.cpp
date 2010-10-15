#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

int main()
{
  CDT cdt;

  Vertex_handle va = cdt.insert(Point(-4,0));
  Vertex_handle vb = cdt.insert(Point(0,-1));
  Vertex_handle vc = cdt.insert(Point(4,0));
  Vertex_handle vd = cdt.insert(Point(0,1));
  cdt.insert(Point(2, 0.6));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::cout << "Meshing the triangulation with default criterias..."
            << std::endl;

  Mesher mesher(cdt);
  mesher.refine_mesh();

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::cout << "Meshing with new criterias..." << std::endl;
  // 0.125 is the default shape bound. It corresponds to abound 20.6 degree.
  // 0.5 is the upper bound on the length of the longuest edge.
  // See reference manual for Delaunay_mesh_size_traits_2<K>.
  mesher.set_criteria(Criteria(0.125, 0.5));
  mesher.refine_mesh();

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
}
