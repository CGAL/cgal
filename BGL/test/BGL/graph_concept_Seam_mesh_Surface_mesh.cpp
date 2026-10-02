#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>

#include <cstdlib>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor         SM_vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor       SM_halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor           SM_edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor           SM_face_descriptor;

typedef Mesh::Property_map<SM_edge_descriptor, bool>            Seam_edge_pmap;
typedef Mesh::Property_map<SM_vertex_descriptor, bool>          Seam_vertex_pmap;
typedef CGAL::Seam_mesh<Mesh, Seam_edge_pmap>                   Owned_seam_mesh;
typedef CGAL::Seam_mesh<Mesh, Seam_edge_pmap, Seam_vertex_pmap> Seam_mesh;

typedef boost::graph_traits< Seam_mesh > Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::halfedge_descriptor halfedge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;
typedef Traits::face_descriptor face_descriptor;

void concept_check_seam_surface_mesh()
{
  boost::function_requires< boost::GraphConcept<Seam_mesh> >();
  boost::function_requires< boost::VertexListGraphConcept<Seam_mesh> >();
  boost::function_requires< boost::EdgeListGraphConcept<Seam_mesh> >();
  boost::function_requires< boost::IncidenceGraphConcept<Seam_mesh> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Seam_mesh> >();
  boost::function_requires< boost::BidirectionalGraphConcept<Seam_mesh> >();

  boost::function_requires< CGAL::HalfedgeGraphConcept<Seam_mesh> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Seam_mesh> >();
  boost::function_requires< CGAL::FaceGraphConcept<Seam_mesh> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Seam_mesh> >();

  // null
  boost::graph_traits<Seam_mesh>::null_vertex();
  boost::graph_traits<Seam_mesh>::null_halfedge();
  boost::graph_traits<Seam_mesh>::null_face();
}

bool test_computed_seam_vertices()
{
  Mesh mesh;
  const SM_vertex_descriptor v0 = mesh.add_vertex(K::Point_3(0, 0, 0));
  const SM_vertex_descriptor v1 = mesh.add_vertex(K::Point_3(1, 0, 0));
  const SM_vertex_descriptor v2 = mesh.add_vertex(K::Point_3(1, 1, 0));
  const SM_vertex_descriptor v3 = mesh.add_vertex(K::Point_3(0, 1, 0));
  if(mesh.add_face(v0, v1, v2) == Mesh::null_face() ||
     mesh.add_face(v0, v2, v3) == Mesh::null_face())
    return false;

  Seam_edge_pmap seam_edges =
    mesh.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  const auto diagonal_result = CGAL::edge(v0, v2, mesh);
  if(!diagonal_result.second)
    return false;
  const SM_edge_descriptor diagonal = diagonal_result.first;
  put(seam_edges, diagonal, true);

  Owned_seam_mesh seam_mesh(mesh, seam_edges);
  if(!seam_mesh.has_on_seam(v0) || !seam_mesh.has_on_seam(v2) ||
     seam_mesh.has_on_seam(v1) || seam_mesh.has_on_seam(v3) ||
     seam_mesh.number_of_seam_edges() != 1 ||
     num_vertices(seam_mesh) != 6 ||
     num_edges(seam_mesh) != num_edges(mesh) + 1 ||
     num_halfedges(seam_mesh) != num_halfedges(mesh) + 2)
    return false;

  const Owned_seam_mesh copied(seam_mesh);
  if(!copied.has_on_seam(v0) || !copied.has_on_seam(v2))
    return false;

  CGAL::Seam_mesh deduced(mesh, seam_edges);
  if(!deduced.has_on_seam(v0) || !deduced.has_on_seam(v2))
    return false;

  Seam_vertex_pmap seam_vertices =
    mesh.add_property_map<SM_vertex_descriptor, bool>("v:on_seam", false).first;
  put(seam_vertices, v0, true);
  put(seam_vertices, v2, true);
  Seam_mesh with_precomputed_maps(mesh, seam_edges, seam_vertices);
  if(with_precomputed_maps.number_of_seam_edges() != 1 ||
     num_edges(with_precomputed_maps) != num_edges(mesh) + 1)
    return false;

  put(seam_edges, diagonal, false);
  Owned_seam_mesh initially_empty(mesh, seam_edges);
  if(!initially_empty.add_seam(v0, v2) ||
     !initially_empty.has_on_seam(v0) || !initially_empty.has_on_seam(v2))
    return false;

  put(seam_edges, diagonal, false);
  put(seam_vertices, v0, false);
  put(seam_vertices, v2, false);
  Seam_mesh with_external_map(mesh, seam_edges, seam_vertices);
  return with_external_map.add_seam(v0, v2) &&
         get(seam_vertices, v0) && get(seam_vertices, v2);
}

bool test_add_to_precomputed_seams()
{
  Mesh mesh;
  const SM_vertex_descriptor v0 = mesh.add_vertex(K::Point_3(0, 0, 0));
  const SM_vertex_descriptor v1 = mesh.add_vertex(K::Point_3(1, 0, 0));
  const SM_vertex_descriptor v2 = mesh.add_vertex(K::Point_3(1, 1, 0));
  const SM_vertex_descriptor v3 = mesh.add_vertex(K::Point_3(0, 2, 0));
  const SM_vertex_descriptor v4 = mesh.add_vertex(K::Point_3(-1, 1, 0));
  if(mesh.add_face(v0, v1, v2) == Mesh::null_face() ||
     mesh.add_face(v0, v2, v3) == Mesh::null_face() ||
     mesh.add_face(v0, v3, v4) == Mesh::null_face())
    return false;

  Seam_edge_pmap seam_edges =
    mesh.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  const auto first_edge = CGAL::edge(v0, v2, mesh);
  if(!first_edge.second)
    return false;
  put(seam_edges, first_edge.first, true);

  {
    Owned_seam_mesh seam_mesh(mesh, seam_edges);
    if(seam_mesh.number_of_seam_edges() != 1 ||
       !seam_mesh.add_seam(v0, v3) ||
       seam_mesh.number_of_seam_edges() != 2 ||
       !seam_mesh.has_on_seam(v3) || seam_mesh.has_on_seam(v4))
      return false;
  }

  const auto second_edge = CGAL::edge(v0, v3, mesh);
  if(!second_edge.second)
    return false;
  put(seam_edges, second_edge.first, false);
  Seam_vertex_pmap seam_vertices =
    mesh.add_property_map<SM_vertex_descriptor, bool>("v:on_seam", false).first;
  put(seam_vertices, v0, true);
  put(seam_vertices, v2, true);

  Seam_mesh seam_mesh(mesh, seam_edges, seam_vertices);
  return seam_mesh.number_of_seam_edges() == 1 &&
         seam_mesh.add_seam(v0, v3) &&
         seam_mesh.number_of_seam_edges() == 2 &&
         get(seam_vertices, v3) && !get(seam_vertices, v4);
}

int main()
{
  concept_check_seam_surface_mesh();
  return test_computed_seam_vertices() && test_add_to_precomputed_seams()
       ? EXIT_SUCCESS : EXIT_FAILURE;
}
