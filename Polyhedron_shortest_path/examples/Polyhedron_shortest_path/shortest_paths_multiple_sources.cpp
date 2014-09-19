#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polyhedron_shortest_path.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef CGAL::Polyhedron_shortest_path_traits<Kernel, Polyhedron_3> Traits;
typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

Traits::Barycentric_coordinate random_coordinate(CGAL::Random& rand)
{
  typename Traits::Construct_barycentric_coordinate construct_barycentric_coordinate;
  Traits::FT u = rand.uniform_real(Traits::FT(0.0), Traits::FT(1.0));
  Traits::FT v = rand.uniform_real(Traits::FT(0.0), Traits::FT(Traits::FT(1.0) - u));
  return construct_barycentric_coordinate(u, v, Traits::FT(Traits::FT(1.0) - u - v));
}

int main(int argc, char** argv)
{
  Polyhedron_3 polyhedron;
  
  std::ifstream inStream(argv[1]);
  
  inStream >> polyhedron;
  
  inStream.close();
  
  CGAL::set_halfedgeds_items_id(polyhedron);
  
  const size_t randSeed = argc > 2 ? std::atoi(argv[2]) : 6065626;
  CGAL::Random rand(randSeed);
  
  face_iterator facesStart, facesEnd;
  boost::tie(facesStart, facesEnd) = faces(polyhedron);
  
  std::vector<face_descriptor> faceList;
  
  for (face_iterator facesCurrent = facesStart; facesCurrent != facesEnd; ++facesCurrent)
  {
    faceList.push_back(*facesCurrent);
  }

  const size_t numSamplePoints = 30;
  
  std::vector<Polyhedron_shortest_path::Face_location> faceLocations;
  
  for (size_t i = 0; i < numSamplePoints; ++i)
  {
    faceLocations.push_back(Polyhedron_shortest_path::Face_location(faceList[rand.get_int(0, CGAL::num_faces(polyhedron))], random_coordinate(rand)));
  }
  
  Traits traits;
  Polyhedron_shortest_path shortestPaths(polyhedron, traits);

  shortestPaths.construct_sequence_tree(faceLocations.begin(), faceLocations.end());
  
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
