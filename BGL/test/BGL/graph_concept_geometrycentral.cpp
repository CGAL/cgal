
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/surface_mesh_factories.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <CGAL/boost/graph/graph_traits_geometrycentral.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <map>
#include <unordered_map>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_concepts.hpp>

using SurfaceMesh = geometrycentral::surface::ManifoldSurfaceMesh;
using VertexPositionGeometry = geometrycentral::surface::VertexPositionGeometry;

using vertex_descriptor = boost::graph_traits<SurfaceMesh>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<SurfaceMesh>::edge_descriptor;
using halfedge_descriptor = boost::graph_traits<SurfaceMesh>::halfedge_descriptor;
using face_descriptor = boost::graph_traits<SurfaceMesh>::face_descriptor;


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

void test_maps(const SurfaceMesh& sm)
{
  std::map<vertex_descriptor,int> vim;
  std::unordered_map<vertex_descriptor,int> vium;
  std::map<halfedge_descriptor,int> him;
  std::unordered_map<halfedge_descriptor,int> hium;
  std::map<edge_descriptor,int> eim;
  std::unordered_map<edge_descriptor,int> eium;
  std::map<face_descriptor,int> fim;
  std::unordered_map<face_descriptor,int> fium;

  int i = 0, j = 0;
  for(auto v : vertices(sm)){
    vim[v] = i++;
    vium[v] = j++;
  }
  i = 0; j = 0;
  for(auto v : vertices(sm)){
    assert(vim[v] == i++);
    assert(vium[v] == j++);
  }

  i = 0; j = 0;
  for(auto v : halfedges(sm)){
    him[v] = i++;
    hium[v] = j++;
  }
  i = 0; j = 0;
  for(auto v : halfedges(sm)){
    assert(him[v] == i++);
    assert(hium[v] == j++);
  }
  i = 0; j = 0;
  for(auto v : edges(sm)){
    eim[v] = i++;
    eium[v] = j++;
  }
  i = 0; j = 0;
  for(auto v : edges(sm)){
    assert(eim[v] == i++);
    assert(eium[v] == j++);
  }

  i = 0; j = 0;
  for(auto v : faces(sm)){
    fim[v] = i++;
    fium[v] = j++;
  }
  i = 0; j = 0;
  for(auto v : faces(sm)){
    assert(fim[v] == i++);
    assert(fium[v] == j++);
  }
}

int main()
{
  test_concepts();

  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = geometrycentral::surface::readManifoldSurfaceMesh("./data/tetrahedron.off");

  test_maps(*mesh);

  std::cout  << num_vertices(*mesh) << " vertices" << std::endl;

  for(auto h : halfedges(*mesh)){
    assert(halfedge(edge(h,*mesh),*mesh) == h);
    if(CGAL::is_border(h,*mesh)){ std::cout  << "b" << std::endl;}else{ std::cout  << "nb" << std::endl;}
  }
  num_vertices(*mesh);

  for( auto v : mesh->vertices()){
    assert(v == v.halfedge().twin().tipVertex());
    for(auto e : out_edges(v,*mesh)){
      assert(v == source(e,*mesh));
    }
    for(auto e : in_edges(v,*mesh)){
      assert(v == target(e,*mesh));
    }
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
                         //     .vertex_index_map(boost::associative_property_map<geometrycentral::surface::VertexData<std::size_t>>(mesh->getVertexIndices()))
                              );

  for(auto v : mesh->vertices()){
    std::cout << v << " at distance " << distance[v] << std::endl;
  }

  boost::property_map<SurfaceMesh, boost::vertex_point_t>::const_type vpm(*geometry);
  for(auto he : halfedges(*mesh)){
    std::cout << "edge length: " << CGAL::Polygon_mesh_processing::edge_length(he,*mesh, CGAL::parameters::vertex_point_map(vpm)) << std::endl;
  }

  geometry->requireEdgeLengths();

  geometrycentral::surface::VertexData<vertex_descriptor> predecessor(*mesh);
  boost::prim_minimum_spanning_tree(*mesh,
                                    boost::make_assoc_property_map(predecessor),
                                    boost::root_vertex(source).weight_map(boost::make_assoc_property_map(geometry->edgeLengths)));

  return 0;
}
