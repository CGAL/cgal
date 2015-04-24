#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>

#include <CGAL/Random.h>
#include <CGAL/Default.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh_shortest_path.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
// default property maps
typedef boost::property_map<Polyhedron_3,
                            boost::vertex_external_index_t>::type  Vertex_index_map;
typedef boost::property_map<Polyhedron_3,
                            CGAL::halfedge_external_index_t>::type Halfedge_index_map;
typedef boost::property_map<Polyhedron_3,
                            CGAL::face_external_index_t>::type     Face_index_map;
typedef CGAL::Surface_mesh_shortest_path<Traits,
                                         Vertex_index_map,
                                         Halfedge_index_map,
                                         Face_index_map>  Surface_mesh_shortest_path;
typedef boost::graph_traits<Polyhedron_3> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::halfedge_iterator halfedge_iterator;
typedef Graph_traits::face_iterator face_iterator;


int main(int argc, char** argv)
{
  // read input polyhedron
  Polyhedron_3 polyhedron;
  std::ifstream input((argc>1)?argv[1]:"data/elephant.off");
  input >> polyhedron;
  input.close();

  // pick up a random face
  const size_t randSeed = argc > 2 ? std::atoi(argv[2]) : 7915421;
  CGAL::Random rand(randSeed);
  const int target_face_index = rand.get_int(0, num_faces(polyhedron));
  face_iterator face_it = faces(polyhedron).first;
  std::advance(face_it,target_face_index);
  // ... and define a barycentric coordinate inside the face
  Traits::Barycentric_coordinate face_location = {{0.25, 0.5, 0.25}};

  // construct a shortest path query object and add a source point
  // Note that the external index property map are automatically initialized
  Surface_mesh_shortest_path shortest_paths(polyhedron,
                                            get(boost::vertex_external_index, polyhedron),
                                            get(CGAL::halfedge_external_index, polyhedron),
                                            get(CGAL::face_external_index, polyhedron),
                                            get(CGAL::vertex_point, polyhedron));
  shortest_paths.add_source_point(*face_it, face_location);

  // For all vertices in the polyhedron, compute the points of
  // the shortest path to the source point and write them
  // into a file readable using the CGAL Polyhedron demo
  std::ofstream output("shortest_paths_no_id.cgal");
  vertex_iterator vit, vit_end;
  for ( boost::tie(vit, vit_end) = vertices(polyhedron);
        vit != vit_end; ++vit)
  {
    std::vector<Traits::Point_3> points;
    shortest_paths.shortest_path_points_to_source_points(*vit, std::back_inserter(points));

    // print the points
    output << points.size() << " ";
    for (std::size_t i = 0; i < points.size(); ++i)
      output << " " << points[i];
    output << std::endl;
  }

  return 0;
}
