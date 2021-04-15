#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;


int main()
{
  CGAL::Surface_mesh<Kernel::Point_3> mesh;
  std::ifstream input("data/test_mesh_6.off");
  input >> mesh;
  input.close();

  for (Triangle_mesh::Vertex_index v1 : vertices(mesh))
    for (Triangle_mesh::Vertex_index v2 : vertices(mesh))
    {
      Surface_mesh_shortest_path shortest_paths(mesh);
      shortest_paths.add_source_point(v1);
      double dist = shortest_paths.shortest_distance_to_source_points(v2).first;
      assert (dist==0 || v1!=v2);
      CGAL_USE(dist);
    }

  return 0;
}
