// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef GraphTraits::halfedge_iterator halfedge_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

template <class Traits>
struct Point_sequence_collector
{
  typedef typename Traits::Point_3 Point_3;
  
  std::vector<Point_3> m_points;
  
  void point(const Point_3& point)
  {
    m_points.push_back(point);
  }
};

int main(int argc, char** argv)
{
  Polyhedron_3 polyhedron;
  
  std::ifstream inStream("data/elephant.off");
  
  inStream >> polyhedron;
  
  inStream.close();

  const size_t targetFaceIndex = 432;
  
  face_iterator facesCurrent, facesEnd;
  boost::tie(facesCurrent, facesEnd) = CGAL::faces(polyhedron);
  
  size_t currentFaceIndex = 0;
  
  while (currentFaceIndex < targetFaceIndex)
  {
    ++facesCurrent;
    ++currentFaceIndex;
  }
  
  face_descriptor targetFace = *facesCurrent;
  
  Traits::Barycentric_coordinate faceLocation(Traits::FT(0.25), Traits::FT(0.5), Traits::FT(0.25));
  
  Traits traits;
  Polyhedron_shortest_path shortestPaths(traits, polyhedron);

  shortestPaths.compute_shortest_paths(targetFace, faceLocation);
  
  vertex_iterator verticesCurrent, verticesEnd;
  
  for (boost::tie(verticesCurrent, verticesEnd) = boost::vertices(polyhedron); verticesCurrent != verticesEnd; ++verticesCurrent)
  {
    Point_sequence_collector<Traits> collector;
    
    shortestPaths.shortest_path_points(*verticesCurrent, collector);
    
    std::cout << collector.m_points.size();
    
    for (size_t i = 0; i < collector.m_points.size(); ++i)
    {
      std::cout << " " << collector.m_points[i];
    }
    
    std::cout << std::endl;
  }
  
  return 0;
}