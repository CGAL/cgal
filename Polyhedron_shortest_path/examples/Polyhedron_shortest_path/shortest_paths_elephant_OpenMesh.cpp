// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>

#include <CGAL/boost/graph/iterator.h>

#include <CGAL/property_map.h>

#include <fstream>
#include <iterator>

#define UNUSED(X) (void)sizeof(X)

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::PolyMesh_ArrayKernelT< /* MyTraits*/> Mesh;

typedef boost::graph_traits<Mesh> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef GraphTraits::halfedge_iterator halfedge_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

typedef std::map<vertex_descriptor, std::size_t> vertex_int_map;
typedef std::map<halfedge_descriptor, std::size_t> halfedge_int_map;
typedef std::map<face_descriptor, std::size_t> face_int_map;

typedef boost::associative_property_map<vertex_int_map> VertexIndexMap;
typedef boost::associative_property_map<halfedge_int_map> HalfedgeIndexMap;
typedef boost::associative_property_map<face_int_map> FaceIndexMap;

typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Mesh> Traits;
typedef CGAL::Polyhedron_shortest_path<Traits, 
  VertexIndexMap, 
  HalfedgeIndexMap, 
  FaceIndexMap, 
  boost::property_map<Mesh, CGAL::vertex_point_t>::type> Polyhedron_shortest_path;


int main(int argc, char** argv)
{
  UNUSED(argc);
  UNUSED(argv);

  Mesh polyhedron;

  OpenMesh::IO::read_mesh(polyhedron, "data/elephant.off");

  const size_t targetFaceIndex = 432;
  
  face_iterator facesCurrent, facesEnd;
  boost::tie(facesCurrent, facesEnd) = faces(polyhedron);
  
  size_t currentFaceIndex = 0;
  
  face_int_map underlyingFaceMap;
  FaceIndexMap faceIndexMap(underlyingFaceMap);
  
  face_descriptor targetFace;
  
  while (currentFaceIndex < num_faces(polyhedron))
  {
    if (currentFaceIndex == targetFaceIndex)
    {
      targetFace = *facesCurrent;
    }
    
    faceIndexMap[*facesCurrent] = currentFaceIndex;
    ++facesCurrent;
    
    ++currentFaceIndex;
  }

  vertex_int_map underlyingVertexMap;
  VertexIndexMap vertexIndexMap(underlyingVertexMap);
  size_t currentVertexIndex = 0;
  vertex_iterator verticesCurrent, verticesEnd;
  boost::tie(verticesCurrent, verticesEnd) = vertices(polyhedron);
  
  while (currentVertexIndex < num_vertices(polyhedron))
  {
    vertexIndexMap[*verticesCurrent] = currentVertexIndex;
    ++currentVertexIndex;
    ++verticesCurrent;
  }
  
  halfedge_int_map underlyingHalfedgeMap;
  HalfedgeIndexMap halfedgeIndexMap(underlyingHalfedgeMap);
  size_t currentHalfedgeIndex = 0;
  halfedge_iterator halfedgesCurrent, halfedgesEnd;
  boost::tie(halfedgesCurrent, halfedgesEnd) = halfedges(polyhedron);
  
  while (currentHalfedgeIndex < num_halfedges(polyhedron))
  {
    halfedgeIndexMap[*halfedgesCurrent] = currentHalfedgeIndex;
    ++currentHalfedgeIndex;
    ++halfedgesCurrent;
  }
  
  Traits::Barycentric_coordinate faceLocation(Traits::FT(0.25), Traits::FT(0.5), Traits::FT(0.25));
  
  Traits traits;
  Polyhedron_shortest_path shortestPaths(polyhedron, vertexIndexMap, halfedgeIndexMap, faceIndexMap, boost::get(boost::vertex_point, polyhedron), traits);

  shortestPaths.m_debugOutput = true;
  
  shortestPaths.construct_sequence_tree(targetFace, faceLocation);
  
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
