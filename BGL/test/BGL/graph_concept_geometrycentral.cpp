
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/surface_mesh_factories.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <CGAL/boost/graph/graph_traits_geometrycentral.h>

#include <map>
#include <unordered_map>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_concepts.hpp>

using SurfaceMesh = geometrycentral::surface::ManifoldSurfaceMesh;
using VertexPositionGeometry = geometrycentral::surface::VertexPositionGeometry;

using vertex_descriptor = boost::graph_traits<SurfaceMesh>::vertex_descriptor;

using vertex_iterator = boost::graph_traits<SurfaceMesh>::vertex_iterator;

void test_concepts()
{
  boost::function_requires< boost::GraphConcept<SurfaceMesh> >();
  boost::function_requires< boost::AdjacencyGraphConcept<SurfaceMesh> >();
  boost::function_requires< boost::VertexListGraphConcept<SurfaceMesh> >();
  boost::function_requires< boost::EdgeListGraphConcept<SurfaceMesh> >();
  boost::function_requires< boost::IncidenceGraphConcept<SurfaceMesh> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<SurfaceMesh> >();
  boost::function_requires< boost::BidirectionalGraphConcept<SurfaceMesh> >();
}


int main()
{
  test_concepts();

  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = geometrycentral::surface::readManifoldSurfaceMesh("./data/tetrahedron.off");
  num_vertices(*mesh);

  for( auto v : mesh->vertices()){
    assert(v == v.halfedge().twin().tipVertex());
  }

  int index = 0;
  std::map<vertex_descriptor,int> vim;
  std::unordered_map<vertex_descriptor,int> vium;
  std::pair<vertex_iterator, vertex_iterator> vr = vertices(*mesh);
  for (vertex_descriptor vd : vertices(*mesh)){
    vim[vd]= index;
    vium[vd]= index;
    ++index;
  }

  // We can record the distance in a map, but also in VertexData
  // std::map<vertex_descriptor,double> distance;
  geometrycentral::surface::VertexData<double> distance(*mesh);

  auto source = *(vertices(*mesh).first);
  boost::breadth_first_search(*mesh,
                              source,
                              boost::visitor(boost::make_bfs_visitor(boost::record_distances(boost::make_assoc_property_map(distance), boost::on_tree_edge())))
                              .vertex_index_map(boost::associative_property_map<geometrycentral::surface::VertexData<std::size_t>>(mesh->getVertexIndices())));

  for(auto v : mesh->vertices()){
    std::cout << v << " at distance " << distance[v] << std::endl;
  }

  return 0;
}