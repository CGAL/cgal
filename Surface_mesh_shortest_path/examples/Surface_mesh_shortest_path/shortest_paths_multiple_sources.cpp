#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_shortest_path.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef Surface_mesh_shortest_path::Face_location Face_location;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;
typedef Graph_traits::face_descriptor face_descriptor;

int main(int argc, char** argv)
{
  // read input tmesh
  Triangle_mesh tmesh;
  std::ifstream input((argc>1)?argv[1]:"data/elephant.off");
  input >> tmesh;
  input.close();

  // pick up some source points inside faces,
  const unsigned int randSeed = argc > 2 ? boost::lexical_cast<unsigned int>(argv[2]) : 7915421;
  CGAL::Random rand(randSeed);
  // by copying the faces in a vector to get a direct access to faces
  std::size_t nb_faces=num_faces(tmesh);
  face_iterator fit, fit_end;
  boost::tie(fit, fit_end) = faces(tmesh);
  std::vector<face_descriptor> face_vector(fit, fit_end);
  // and creating a vector of Face_location objects
  const std::size_t nb_source_points = 30;
  Traits::Barycentric_coordinates face_location = {{0.25, 0.5, 0.25}};
  std::vector<Face_location> faceLocations(nb_source_points, Face_location(face_descriptor(), face_location));
  for (std::size_t i = 0; i < nb_source_points; ++i)
  {
    faceLocations[i].first=face_vector[rand.get_int(0, static_cast<int>(nb_faces))];
  }

  // construct a shortest path query object and add a range of source points
  Surface_mesh_shortest_path shortest_paths(tmesh);
  shortest_paths.add_source_points(faceLocations.begin(), faceLocations.end());

  // For all vertices in the tmesh, compute the points of
  // the shortest path to the source point and write them
  // into a file readable using the CGAL Tmesh demo
  std::ofstream output("shortest_paths_multiple_sources.cgal");
  vertex_iterator vit, vit_end;
  for ( boost::tie(vit, vit_end) = vertices(tmesh);
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
