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

#include <CGAL/Default.h>

#include <fstream>
#include <iterator>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
typedef boost::property_map<Polyhedron_3, boost::vertex_external_index_t>::type VertexIndexMap;
typedef boost::property_map<Polyhedron_3, CGAL::halfedge_external_index_t>::type HalfedgeIndexMap;
typedef boost::property_map<Polyhedron_3, CGAL::face_external_index_t>::type FaceIndexMap;
typedef CGAL::Polyhedron_shortest_path<Traits, VertexIndexMap, HalfedgeIndexMap, FaceIndexMap, CGAL::Default> Polyhedron_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

int main()
{
  Traits::Construct_barycentric_coordinate construct_barycentric_coordinate;
  
  Polyhedron_3 polyhedron;
  
  std::ifstream inStream("data/elephant.off");
  
  inStream >> polyhedron;
  
  inStream.close();

  const size_t targetFaceIndex = 432;
  
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
  Polyhedron_shortest_path shortestPaths(polyhedron, 
    CGAL::get(boost::vertex_external_index, polyhedron), 
    CGAL::get(CGAL::halfedge_external_index, polyhedron),
    CGAL::get(CGAL::face_external_index, polyhedron),
    CGAL::get(CGAL::vertex_point, polyhedron), 
    traits);

  shortestPaths.construct_sequence_tree(targetFace, faceLocation);
  
  vertex_iterator verticesCurrent, verticesEnd;
  
  std::ofstream outPaths("polylines.cgal");
  
  for (boost::tie(verticesCurrent, verticesEnd) = boost::vertices(polyhedron); verticesCurrent != verticesEnd; ++verticesCurrent)
  {
    std::vector<Traits::Point_3> points;
    
    shortestPaths.shortest_path_points_to_source_points(*verticesCurrent, std::back_inserter(points));
    
    outPaths << points.size();
    
    for (size_t i = 0; i < points.size(); ++i)
    {
      outPaths << " " << points[i];
    }
    
    outPaths << std::endl;
  }
  
  outPaths.close();
  
  return 0;
}
