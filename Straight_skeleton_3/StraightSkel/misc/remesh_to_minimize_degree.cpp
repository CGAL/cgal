

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_traits_3.h>

#include <fstream>
#include <iostream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;
using Plane_3 = K::Plane_3;
using Vector_3 = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
using degree_size_type = boost::graph_traits<Mesh>::degree_size_type;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  Mesh sm;
  CGAL::IO::read_polygon_mesh(argv[1], sm);

  // declare vectors to store mesh properties
  std::vector<std::size_t> region_ids(num_faces(sm));
  std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices

  std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions
  boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected

  // detect planar regions in the mesh
  std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(sm,
                                           CGAL::make_random_access_property_map(region_ids),
                                           CGAL::parameters::cosine_of_maximum_angle(0.98)
                                                            .region_primitive_map(plane_map)
                                                            .maximum_distance(0.001));

  // detect corner vertices on the boundary of planar regions
  std::size_t nb_corners =
    PMP::detect_corners_of_regions(sm,
                                   CGAL::make_random_access_property_map(region_ids),
                                   nb_regions,
                                   CGAL::make_random_access_property_map(corner_id_map),
                                   CGAL::parameters::cosine_of_maximum_angle(0.98).
                                                     maximum_distance(0.001).
                                                     edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

#if 1
  // run the remeshing algorithm using filled properties
  std::vector<Vector_3> normal_map(num_faces(sm));
  for(face_descriptor f : faces(sm))
    normal_map[f] = plane_map[f].orthogonal_vector();

  // mark all vertices as corners
  nb_corners = 0;
  for(vertex_descriptor v : vertices(sm))
    corner_id_map[v] = nb_corners++;

  Mesh out;
  PMP::remesh_almost_planar_patches(sm,
                                    out,
                                    nb_regions, nb_corners,
                                    CGAL::make_random_access_property_map(region_ids),
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::make_random_access_property_map(ecm),
                                    CGAL::parameters::patch_normal_map(CGAL::make_random_access_property_map(normal_map)),
                                    CGAL::parameters::do_not_triangulate_faces(false));
    sm = out;
#else
  std::cout << "Now, flip stuff" << std::endl;

  // Let's flip edges such that the degree of the vertices is minimized.
  // This is a simple greedy algorithm that flips edges until no more
  // improvement can be made.
  // It is not guaranteed to find the optimal solution, but it is fast
  // and works well in practice.
  // The criterion is the following:
  // an edge is to be flipped if the sum of the degrees of its
  // incident vertices is greater than the sum of the degrees of
  // the vertices that would be incident to the edge after the flip.
  // Also, flipping should only happens if:
  // - all four vertices should be in convex position.
  // - the edge does not already exist.

  // - the edge is not a boundary edge.
  // - the edge is not a crease edge.

  // Create a property map to store the degree of each vertex.
  auto degree_map = sm.add_property_map<vertex_descriptor, degree_size_type>("v:degree", 0).first;

  // Initialize the degree of each vertex.
  for (vertex_descriptor v : vertices(sm)) {
    degree_map[v] = degree(v, sm);
  }

  // Flip edges until no more improvement can be made.
  int residual_threshold = 4;
  while (residual_threshold > 2) {
    std::set<edge_descriptor> edges_to_test(edges(sm).begin(), edges(sm).end());
    while (!edges_to_test.empty()) {
      edge_descriptor e = *edges_to_test.begin();
      edges_to_test.erase(edges_to_test.begin());
      if (is_border(e, sm)) {
        continue;
      }

      halfedge_descriptor h = halfedge(e, sm);
      CGAL_assertion(CGAL::is_triangle(h, sm) &&
                     CGAL::is_triangle(opposite(h, sm), sm));

      // the current edge is v0-v1
      vertex_descriptor v0 = target(h, sm);
      vertex_descriptor v1 = source(h, sm);
      vertex_descriptor v2 = target(next(h, sm), sm);
      vertex_descriptor v3 = target(next(opposite(h, sm), sm), sm);
      CGAL_assertion(v0 != v1 && v0 != v2 && v0 != v3 && v1 != v2 && v1 != v3 && v2 != v3);

      int degree_v0 = static_cast<int>(degree_map[v0]);
      int degree_v1 = static_cast<int>(degree_map[v1]);
      int degree_v2 = static_cast<int>(degree_map[v2]);
      int degree_v3 = static_cast<int>(degree_map[v3]);

      int residual = degree_map[v0] + degree_map[v1] - degree_map[v2] - degree_map[v3];

      std::cout << v0 << " (" << degree_v0 << ") "
                << v1 << " (" << degree_v1 << ") "
                << v2 << " (" << degree_v2 << ") "
                << v3 << " (" << degree_v3 << ") "
                << " residual: " << residual << std::endl;

      // if the edge is not in a planar region, ignore it
      if (region_ids[face(h, sm)] != region_ids[face(opposite(h, sm), sm)]) {
        std::cout << " don't flip (not in the same region)" << std::endl;
        continue;
      }

      // check that all vertices are in convex position


      if (halfedge(v2, v3, sm).second) {
        std::cout << " don't flip (flipped edge exists)" << std::endl;
        continue;
      }

      if (residual <= residual_threshold) {
        std::cout << " don't flip (below threshold)" << std::endl;
        continue;
      }

      // check that the areas would be the same before and after flip
      FT area_b1 = CGAL::approximate_sqrt(CGAL::squared_area(sm.point(v0), sm.point(v1), sm.point(v2)));
      FT area_b2 = CGAL::approximate_sqrt(CGAL::squared_area(sm.point(v0), sm.point(v3), sm.point(v1)));
      FT area_a1 = CGAL::approximate_sqrt(CGAL::squared_area(sm.point(v0), sm.point(v3), sm.point(v2)));
      FT area_a2 = CGAL::approximate_sqrt(CGAL::squared_area(sm.point(v2), sm.point(v3), sm.point(v1)));
      if (CGAL::abs((area_b1 + area_b2) - (area_a1 + area_a2)) > 0.001) {
        std::cout << " don't flip (areas differ)" << std::endl;
        continue;
      }

      // reject if one of the faces becomes really thin
      if (area_a1 < 0.0000001 || area_a2 < 0.0000001) {
        std::cout << " don't flip (area too small)" << std::endl;
        continue;
      }

      std::cout << " flip" << std::endl;

      // Flip the edge.
      CGAL::Euler::flip_edge(h, sm);

      // Update the degree of each vertex
      degree_map[v0] = degree(v0, sm);
      degree_map[v1] = degree(v1, sm);
      degree_map[v2] = degree(v2, sm);
      degree_map[v3] = degree(v3, sm);

      // Update the edges to flip
      for (halfedge_descriptor h : halfedges_around_face(halfedge(v0, sm), sm)) {
        edges_to_test.insert(edge(h, sm));
      }
      for (halfedge_descriptor h : halfedges_around_face(halfedge(v1, sm), sm)) {
        edges_to_test.insert(edge(h, sm));
      }
      for (halfedge_descriptor h : halfedges_around_face(halfedge(v2, sm), sm)) {
        edges_to_test.insert(edge(h, sm));
      }
      for (halfedge_descriptor h : halfedges_around_face(halfedge(v3, sm), sm)) {
        edges_to_test.insert(edge(h, sm));
      }

      edges_to_test.erase(e);

      CGAL::IO::write_polygon_mesh("results/last_flip.off", sm, CGAL::parameters::stream_precision(17));
      CGAL_assertion(!PMP::does_self_intersect(sm));
    }

    std::cout << "residual threshold: " << residual_threshold << std::endl;

    degree_size_type max_degree = 0;
    for (vertex_descriptor v : vertices(sm)) {
      max_degree = (std::max)(max_degree, degree(v, sm));
    }
    std::cout << "Max degree: " << max_degree << std::endl;

    CGAL::IO::write_polygon_mesh("results/flipped_" + std::to_string(residual_threshold) + ".off",
                                 sm, CGAL::parameters::stream_precision(17));

    --residual_threshold;
  }
#endif

  std::cout << "== FINAL ==" << std::endl;
  for (vertex_descriptor v : vertices(sm)) {
    std::cout << sm.point(v) << " (" << degree(v, sm) << ")" << std::endl;
  }

  CGAL::IO::write_polygon_mesh("results/out_remeshed_degree.off", sm, CGAL::parameters::stream_precision(17));

  std::cout << "Done" << std::endl;

  return 0;
}
