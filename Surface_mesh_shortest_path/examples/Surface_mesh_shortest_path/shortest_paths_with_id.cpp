#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>


#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

int main(int argc, char** argv)
{
  Traits::Construct_barycentric_coordinate construct_barycentric_coordinate;
  
  Polyhedron_3 polyhedron;
  
  std::ifstream inStream(argv[1]);
  
  inStream >> polyhedron;
  
  inStream.close();

  CGAL::set_halfedgeds_items_id(polyhedron);
  
  const size_t randSeed = argc > 2 ? std::atoi(argv[2]) : 7915421;
  CGAL::Random rand(randSeed);
  
  const size_t targetFaceIndex = rand.get_int(0, num_faces(polyhedron));
  
  face_iterator facesCurrent, facesEnd;
  boost::tie(facesCurrent, facesEnd) = faces(polyhedron);
  
  size_t currentFaceIndex = 0;
  
  while (currentFaceIndex < targetFaceIndex)
  {
    ++facesCurrent;
    ++currentFaceIndex;
  }

  face_descriptor targetFace = *facesCurrent;
  
  Traits::Barycentric_coordinate faceLocation = construct_barycentric_coordinate(Traits::FT(0.25), Traits::FT(0.5), Traits::FT(0.25));
  
  Traits traits;
  Surface_mesh_shortest_path shortestPaths(polyhedron, traits);
 
  shortestPaths.construct_sequence_tree(targetFace, faceLocation);
  
  vertex_iterator verticesCurrent, verticesEnd;
  
  for (boost::tie(verticesCurrent, verticesEnd) = boost::vertices(polyhedron); verticesCurrent != verticesEnd; ++verticesCurrent)
  {
    std::vector<Traits::Point_3> points;
    
    shortestPaths.shortest_path_points_to_source_points(*verticesCurrent, std::back_inserter(points));
    
    std::cout << points.size();
    
    for (size_t i = 0; i < points.size(); ++i)
    {
      std::cout << " " << points[i];
    }
    
    std::cout << std::endl;
  }
  
  return 0;
}
