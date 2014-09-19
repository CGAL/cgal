#include <cstdlib>
#include <iostream>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;

typedef boost::graph_traits<Mesh> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef GraphTraits::halfedge_iterator halfedge_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Mesh> Traits;
typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;

int main(int argc, char** argv)
{
  Traits::Construct_barycentric_coordinate construct_barycentric_coordinate;
  
  Mesh polyhedron;

  OpenMesh::IO::read_mesh(polyhedron, argv[1]);

  const size_t randSeed = argc > 2 ? std::atoi(argv[2]) : 9601263;
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
  Polyhedron_shortest_path shortestPaths(polyhedron, traits);

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
