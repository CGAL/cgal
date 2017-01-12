#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/internal/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Triangle_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Triangle_mesh::Property_map<Triangle_mesh::Edge_index,bool> Constrained_edge_map;


int main()
{
  {
  // test open clipping with Surface_mesh
  Triangle_mesh tm1, tm2;
  std::ifstream input("data-coref/elephant.off");
  input >> tm1;
  input.close();
  input.open("data-coref/sphere.off");
  input >> tm2;
  input.close();

  Constrained_edge_map ecm1 =
    tm1.add_property_map<Triangle_mesh::Edge_index,bool>("e:cst", false).first;

  PMP::clip(tm1, tm2, false, params::edge_is_constrained_map(ecm1));
  std::ofstream output("clipped_opened.off");
  output << tm1;

  // test open clipping with Polyhedron
  Polyhedron P, Q;
  input.open("data-coref/elephant.off");
  input >> P;
  input.close();
  input.open("data-coref/sphere.off");
  input >> Q;

  PMP::clip(P, Q, false,
            params::face_index_map(get(CGAL::face_external_index, P)).
                    vertex_index_map(get(CGAL::vertex_external_index, P)),
            params::face_index_map(get(CGAL::face_external_index, Q)));
  assert(P.size_of_vertices() == tm1.number_of_vertices());
  }
  {
  Triangle_mesh tm1, tm2;
  std::ifstream input("data-coref/elephant.off");
  input >> tm1;
  input.close();
  input.open("data-coref/sphere.off");
  input >> tm2;

  PMP::clip(tm1, tm2, true);
  std::ofstream output("clipped_closed.off");
  output << tm1;
  }

  return 0;
}
