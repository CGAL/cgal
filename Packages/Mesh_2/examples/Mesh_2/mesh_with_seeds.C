// file examples/Mesh_2/mesh_with_seeds.C
#include <CGAL/basic.h>
#include <iostream>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_traits_2.h>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_mesh_size_traits_2<K> Meshtraits;
typedef CGAL::Constrained_Delaunay_triangulation_2<Meshtraits, Tds,
  CGAL::Exact_predicates_tag> CDT;
typedef CGAL::Delaunay_mesh_2<CDT> Mesh;

typedef CDT::Vertex_handle Vertex_handle;
typedef K::Point_2 Point;

int main(int, char**)
{
  Mesh mesh;

  Vertex_handle va = mesh.insert(Point(2,0));
  Vertex_handle vb = mesh.insert(Point(0,2));
  Vertex_handle vc = mesh.insert(Point(-2,0));
  Vertex_handle vd = mesh.insert(Point(0,-2));

  mesh.insert_constraint(va, vb);
  mesh.insert_constraint(vb, vc);
  mesh.insert_constraint(vc, vd);
  mesh.insert_constraint(vd, va);

  va = mesh.insert(Point(3,3));
  vb = mesh.insert(Point(-3,3));
  vc = mesh.insert(Point(-3,-3));
  vd = mesh.insert(Point(3,0-3));

  mesh.insert_constraint(va, vb);
  mesh.insert_constraint(vb, vc);
  mesh.insert_constraint(vc, vd);
  mesh.insert_constraint(vd, va);

  std::list<Point> list_of_seeds;

  list_of_seeds.push_back(Point(0, 0));

  mesh.set_seeds(list_of_seeds.begin(), list_of_seeds.end());

  std::cout << "Number of vertices: " << mesh.number_of_vertices() 
	    << std::endl;

  std::cout << "Meshing the domain..."
	    << std::endl;
  mesh.refine_mesh();

  std::cout << "Number of vertices: " << mesh.number_of_vertices() 
	    << std::endl;
}
