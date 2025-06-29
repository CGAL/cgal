#define CGAL_SS3_DO_NOT_WARN_ABOUT_APPROXIMATE_SQRT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_traits_3.h>

#include <fstream>
#include <iostream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Plane_3 = K::Plane_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
using degree_size_type = boost::graph_traits<Mesh>::degree_size_type;

namespace PMP = CGAL::Polygon_mesh_processing;

// the point of this function is to apply flips to lower the maximal degree
// of any vertex in the mesh
void remesh_with_flips(int argc, char** argv)
{
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
}

void remesh_with_insertion(int argc, char** argv)
{
  std::cout << "===== remesh_with_insertion =====" << std::endl;

  Mesh sm;
  CGAL::IO::read_polygon_mesh(argv[1], sm);
  CGAL::IO::write_polygon_mesh("results/input.obj", sm, CGAL::parameters::stream_precision(17));

  CGAL::Bbox_3 bbox = PMP::bbox(sm);
  const FT diag_length = CGAL::approximate_sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                                      CGAL::square(bbox.ymax() - bbox.ymin()) +
                                                      CGAL::square(bbox.zmax() - bbox.zmin()));

  const FT cos_of_max_angle = 0.98;
  const FT max_distance = 0.0001 * diag_length;

  // Do one round of remeshing with no trickery to get rid of any nonsense
  std::vector<std::size_t> region_ids(num_faces(sm));
  boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected
  std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
  std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

  // detect planar regions in the mesh
  std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(sm,
                                          CGAL::make_random_access_property_map(region_ids),
                                          CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                            .region_primitive_map(plane_map)
                                                            .maximum_distance(max_distance));

  // detect corner vertices on the boundary of planar regions
  std::size_t nb_corners =
    PMP::detect_corners_of_regions(sm,
                                  CGAL::make_random_access_property_map(region_ids),
                                  nb_regions,
                                  CGAL::make_random_access_property_map(corner_id_map),
                                  CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                    maximum_distance(max_distance).
                                                    edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

  // run the remeshing algorithm using filled properties
  std::vector<Vector_3> normal_map(num_faces(sm));
  for(face_descriptor f : faces(sm))
    normal_map[f] = plane_map[f].orthogonal_vector();

  Mesh out;
  PMP::remesh_almost_planar_patches(sm,
                                    out,
                                    nb_regions, nb_corners,
                                    CGAL::make_random_access_property_map(region_ids),
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::make_random_access_property_map(ecm),
                                    CGAL::parameters::patch_normal_map(CGAL::make_random_access_property_map(normal_map)),
                                    CGAL::parameters::do_not_triangulate_faces(false));

  CGAL::IO::write_polygon_mesh("results/iter_-1.obj", out, CGAL::parameters::stream_precision(17));
  sm = out;

  // Main loop for remeshing with insertion
  for (;;) {
    static int iter = -1;
    CGAL::IO::write_polygon_mesh("results/iter_" + std::to_string(++iter) + ".obj",
                                 sm, CGAL::parameters::stream_precision(17));

    std::cout << "-- Iteration: " << iter << std::endl;
    std::cout << "Number of vertices: " << num_vertices(sm) << std::endl;
    std::cout << "Number of edges: " << num_edges(sm) << std::endl;
    std::cout << "Number of faces: " << num_faces(sm) << std::endl;

    degree_size_type max_degree = 0;
    for (vertex_descriptor v : vertices(sm)) {
      max_degree = (std::max)(max_degree, degree(v, sm));
    }
    std::cout << "Max degree: " << max_degree << std::endl;

    std::vector<std::size_t> region_ids(num_faces(sm));
    boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected

    // detect planar regions in the mesh
    std::size_t nb_regions =
      PMP::region_growing_of_planes_on_faces(sm,
                                            CGAL::make_random_access_property_map(region_ids),
                                            CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                             .region_primitive_map(plane_map)
                                                             .maximum_distance(max_distance));

    // detect edges incident to high degree vertices, and split them
    std::set<edge_descriptor> edges_to_split;
    for (vertex_descriptor v : vertices(sm)) {
      if (degree(v, sm) <= 10) {
         continue;
      }

      // find the longest edge incident to the vertex
      // only consider edges that are incident to two different regions
      edge_descriptor longest_edge;
      FT longest_sq_length = 0;
      for (halfedge_descriptor hh : halfedges_around_target(halfedge(v, sm), sm)) {
        edge_descriptor e = edge(hh, sm);
        if (region_ids[face(hh, sm)] == region_ids[face(opposite(hh, sm), sm)]) {
          continue; // skip edges that are not between two different regions
        }

        FT sq_length = CGAL::squared_distance(sm.point(source(hh, sm)), sm.point(target(hh, sm)));
        if (sq_length > longest_sq_length) {
          longest_sq_length = sq_length;
          longest_edge = e;
        }
      }

      edges_to_split.insert(longest_edge);
      std::cout << "Vertex " << v << " has degree " << degree(v, sm)
                << ", longest edge: " << longest_edge
                << " with length " << CGAL::sqrt(longest_sq_length) << std::endl;
    }

    std::cout << "Edges to split: " << edges_to_split.size() << std::endl;
    if (edges_to_split.empty()) {
      std::cout << "No edges to split, stopping remeshing." << std::endl;
      break; // no edges to split, exit the loop
    }

    for (edge_descriptor e : edges_to_split) {
      halfedge_descriptor h = halfedge(e, sm);
      Point_3 mp = CGAL::midpoint(sm.point(source(h, sm)), sm.point(target(h, sm)));
      halfedge_descriptor new_h = CGAL::Euler::split_edge_and_incident_faces(h, sm);
      sm.point(target(new_h, sm)) = mp;
      std::cout << "Inserted vertex " << target(new_h, sm) << " at midpoint of edge " << e << std::endl;
    }

    region_ids.assign(num_faces(sm), -1);


    // detect planar regions in the mesh
    nb_regions =
      PMP::region_growing_of_planes_on_faces(sm,
                                            CGAL::make_random_access_property_map(region_ids),
                                            CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                             .region_primitive_map(plane_map)
                                                             .maximum_distance(max_distance));

    // detect corner vertices on the boundary of planar regions
    std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

    std::size_t nb_corners =
      PMP::detect_corners_of_regions(sm,
                                    CGAL::make_random_access_property_map(region_ids),
                                    nb_regions,
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                      maximum_distance(max_distance).
                                                      edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    // run the remeshing algorithm using filled properties
    std::vector<Vector_3> normal_map(num_faces(sm));
    for(face_descriptor f : faces(sm))
      normal_map[f] = plane_map[f].orthogonal_vector();

    // mark all vertices as corners to keep them during remeshing
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
  }
}

// the point of this function is:
// - start from a triangle mesh
// - simplify it to the maximum
// - retriangulate it minimally, i.e. only triangulate faces if there is a face
//   with more than 2 high degree vertices
void minimal_triangulation(int argc, char** argv)
{
  Mesh sm;
  CGAL::IO::read_polygon_mesh(argv[1], sm);

  std::cout << "===== minimal_triangulation =====" << std::endl;
  std::cout << "Input:" << std::endl;
  std::cout << "Number of vertices: " << num_vertices(sm) << std::endl;
  std::cout << "Number of edges: " << num_edges(sm) << std::endl;
  std::cout << "Number of faces: " << num_faces(sm) << std::endl;

  CGAL::Bbox_3 bbox = PMP::bbox(sm);
  const FT diag_length = CGAL::approximate_sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                                CGAL::square(bbox.ymax() - bbox.ymin()) +
                                                CGAL::square(bbox.zmax() - bbox.zmin()));

  const FT cos_of_max_angle = 0.98;
  const FT max_distance = 0.0001 * diag_length;

  std::vector<std::size_t> region_ids(num_faces(sm));
  boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected
  std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
  std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

  // detect planar regions in the mesh
  std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(sm,
                                          CGAL::make_random_access_property_map(region_ids),
                                          CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                            .region_primitive_map(plane_map)
                                                            .maximum_distance(max_distance));

  // detect corner vertices on the boundary of planar regions
  std::size_t nb_corners =
    PMP::detect_corners_of_regions(sm,
                                  CGAL::make_random_access_property_map(region_ids),
                                  nb_regions,
                                  CGAL::make_random_access_property_map(corner_id_map),
                                  CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                    .maximum_distance(max_distance)
                                                    .edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

  // run the remeshing algorithm using filled properties
  std::vector<Vector_3> normal_map(num_faces(sm));
  for(face_descriptor f : faces(sm))
    normal_map[f] = plane_map[f].orthogonal_vector();

  Mesh out;
  PMP::remesh_almost_planar_patches(sm,
                                    out,
                                    nb_regions, nb_corners,
                                    CGAL::make_random_access_property_map(region_ids),
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::make_random_access_property_map(ecm),
                                    CGAL::parameters::patch_normal_map(CGAL::make_random_access_property_map(normal_map)),
                                    CGAL::parameters::do_not_triangulate_faces(true));
  sm = out;

  CGAL::IO::write_polygon_mesh("results/iter_-1.obj", sm, CGAL::parameters::stream_precision(17));

  std::cout << "Simplified:" << std::endl;
  std::cout << "Number of vertices: " << vertices(sm).size() << std::endl;
  std::cout << "Number of edges: " << edges(sm).size() << std::endl;
  std::cout << "Number of faces: " << faces(sm).size() << std::endl;

  // add a color map: red for faces with high degree vertices, green for others
  // do not use faces directly but use regions (because with a Surface_mesh, we cannot
  // represent some faces like faces with holes)
  std::vector<bool> region_with_high_degree_corner(nb_regions, false);
  for (face_descriptor f : faces(sm)) {
      for (vertex_descriptor v : CGAL::vertices_around_face(halfedge(f, sm), sm)) {
          std::set<std::size_t> incident_regions;
          for (face_descriptor f : CGAL::faces_around_target(halfedge(v, sm), sm)) {
              incident_regions.insert(region_ids[f]);
          }

          if (incident_regions.size() > 3) {
              for (std::size_t rid : incident_regions) {
                  region_with_high_degree_corner.at(rid) = true;
              }
          }
      }
  }

  auto color_map = sm.add_property_map<face_descriptor, CGAL::Color>("f:color", CGAL::Color(0, 0, 255)).first;
  for (face_descriptor f : faces(sm)) {
      std::size_t rid = region_ids[f];
      if (region_with_high_degree_corner.at(rid)) {
          color_map[f] = CGAL::Color(255, 0, 0);
      } else {
          color_map[f] = CGAL::Color(0, 255, 0);
      }
  }

  CGAL::IO::write_polygon_mesh("results/out_high_degree_faces.ply",
                               sm, CGAL::parameters::stream_precision(17).face_color_map(color_map));

  std::list<face_descriptor> faces_to_consider;
  for (face_descriptor f : faces(sm)) {
      faces_to_consider.push_back(f);
  }

  unsigned int iter = -1;
  while (!faces_to_consider.empty()) {
      face_descriptor f = faces_to_consider.front();
      faces_to_consider.pop_front();

      halfedge_descriptor h = halfedge(f, sm);
      if (CGAL::is_triangle(h, sm)) {
          continue;
      }

      unsigned int high_degree_vertices_n = 0;
      for (vertex_descriptor v : CGAL::vertices_around_face(halfedge(f, sm), sm)) {
          if (degree(v, sm) > 3) {
              ++high_degree_vertices_n;
          }
      }

      if (high_degree_vertices_n <= 2) {
          continue;
      }

      // std::cout << "face " << f << " has " << high_degree_vertices_n << " high degree vertices" << std::endl;

      struct Enlister_visitor
          : public PMP::Triangulate_polygons::Default_visitor
      {
          Enlister_visitor(std::list<face_descriptor>& r) : r_(r) { }

          void after_subface_created(face_descriptor f_new) {
              r_.push_back(f_new);
          }

          std::list<face_descriptor>& r_;
      };

      std::set<face_descriptor> neighbors;
      for (vertex_descriptor v : CGAL::vertices_around_face(h, sm)) {
          for (halfedge_descriptor hh : halfedges_around_target(halfedge(v, sm), sm)) {
              face_descriptor nf = face(hh, sm);
              if (nf != f) {
                  neighbors.insert(nf);
              }
          }
      }

      std::ofstream out_f("results/face_" + std::to_string(++iter) + ".off");
      out_f << "OFF" << std::endl;
      out_f << degree(f, sm) << " 1 0" << std::endl;
      for (vertex_descriptor v : CGAL::vertices_around_face(h, sm)) {
        out_f << sm.point(v) << "\n";
      }
      out_f << degree(f, sm);
      for (degree_size_type i=0; i<degree(f, sm); ++i) {
        out_f << " " << i;
      }
      out_f << std::endl;
      out_f.close();

      std::list<face_descriptor> new_faces;
      Enlister_visitor visitor(new_faces);
      PMP::triangulate_face(f, sm, CGAL::parameters::visitor(visitor));

      faces_to_consider.insert(faces_to_consider.end(),
                               neighbors.begin(), neighbors.end());
      faces_to_consider.insert(faces_to_consider.end(),
                               new_faces.begin(), new_faces.end());
  }

  std::cout << "Output:" << std::endl;
  std::cout << "Number of vertices: " << vertices(sm).size() << std::endl;
  std::cout << "Number of edges: " << edges(sm).size() << std::endl;
  std::cout << "Number of faces: " << faces(sm).size() << std::endl;

  degree_size_type max_degree = 0;
  for (vertex_descriptor v : vertices(sm)) {
    max_degree = (std::max)(max_degree, degree(v, sm));
  }
  std::cout << "Max degree: " << max_degree << std::endl;

  // add a color property map for faces: red for triangles, blue for others
  for (face_descriptor f : faces(sm)) {
    if (CGAL::is_triangle(halfedge(f, sm), sm)) {
      color_map[f] = CGAL::Color(255, 0, 0); // red for triangles
    } else {
      color_map[f] = CGAL::Color(0, 0, 255); // blue for other polygons
    }
  }

  CGAL::IO::write_polygon_mesh("results/out_minimal_triangulation.ply",
                               sm, CGAL::parameters::stream_precision(17).face_color_map(color_map));
}

// This particular functions aims to maximize the number of polygonal faces in the mesh,
// and specifically large polygonal faces.
// The nice formulation would be a system with nonlinear constraints, but since we can't really solve
// that exactly, let's try a greedy method.
void remesh_with_constrained_polygonal_faces(int argc, char** argv)
{
  Mesh sm;
  CGAL::IO::read_polygon_mesh(argv[1], sm);
  CGAL::IO::write_polygon_mesh("results/input.obj", sm, CGAL::parameters::stream_precision(17));

  CGAL::Bbox_3 bbox = PMP::bbox(sm);
  const FT diag_length = CGAL::approximate_sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                                      CGAL::square(bbox.ymax() - bbox.ymin()) +
                                                      CGAL::square(bbox.zmax() - bbox.zmin()));

  const FT cos_of_max_angle = 0.98;
  const FT max_distance = 0.0001 * diag_length;

  // Do one round of remeshing with no trickery to get rid of any nonsense
  std::vector<std::size_t> region_ids(num_faces(sm));
  boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected
  std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
  std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

  // detect planar regions in the mesh
  std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(sm,
                                          CGAL::make_random_access_property_map(region_ids),
                                          CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                            .region_primitive_map(plane_map)
                                                            .maximum_distance(max_distance));

  // detect corner vertices on the boundary of planar regions
  std::size_t nb_corners =
    PMP::detect_corners_of_regions(sm,
                                  CGAL::make_random_access_property_map(region_ids),
                                  nb_regions,
                                  CGAL::make_random_access_property_map(corner_id_map),
                                  CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                    maximum_distance(max_distance).
                                                    edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

  // run the remeshing algorithm using filled properties
  std::vector<Vector_3> normal_map(num_faces(sm));
  for(face_descriptor f : faces(sm))
    normal_map[f] = plane_map[f].orthogonal_vector();

  Mesh out;
  PMP::remesh_almost_planar_patches(sm,
                                    out,
                                    nb_regions, nb_corners,
                                    CGAL::make_random_access_property_map(region_ids),
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::make_random_access_property_map(ecm),
                                    CGAL::parameters::patch_normal_map(CGAL::make_random_access_property_map(normal_map)),
                                    CGAL::parameters::do_not_triangulate_faces(true));

  CGAL::IO::write_polygon_mesh("results/iter_-1.obj", out, CGAL::parameters::stream_precision(17));
  sm = out;

  // dump with colors for face status:
  // - blue for default faces
  // - green for fixed faces
  // - red for not fixable faces
  auto color_map = sm.add_property_map<face_descriptor, CGAL::Color>("f:color", CGAL::Color(0, 0, 255)).first;

  enum class Facet_status
  {
    DEFAULT = 0,
    FIXED,
    NOT_FIXABLE
  };

  // High degree vertex := vertex with strictly more than 3 incident edges.
  // Over constrained face := face with strictly more than 3 high degree vertices
  std::vector<Facet_status> face_status;

  auto flood = [&](Mesh& sm)
  {
    face_status.assign(num_faces(sm), Facet_status::DEFAULT);

    auto vertex_constraints = [&](const vertex_descriptor v) -> unsigned int
    {
      unsigned int n = 0;
      for (halfedge_descriptor h : CGAL::halfedges_around_target(v, sm)) {
        face_descriptor f = face(h, sm);
        if (face_status[f] == Facet_status::FIXED)
          ++n;
      }
      return n;
    };

    for (;;) {
      for (face_descriptor f : faces(sm)) {
        if (face_status[f] == Facet_status::FIXED) {
          color_map[f] = CGAL::Color(0, 255, 0); // green for fixed faces
        } else if (face_status[f] == Facet_status::NOT_FIXABLE) {
          color_map[f] = CGAL::Color(255, 0, 0); // red for not fixable faces
        } else {
          color_map[f] = CGAL::Color(0, 0, 255); // default blue
        }
      }

      static int iter = -1;
      CGAL::IO::write_polygon_mesh("results/iter_" + std::to_string(++iter) + ".ply",
                                   sm, CGAL::parameters::stream_precision(17).face_color_map(color_map));


      // Start of the loop: find the high degree vertex with the largest incident face (#vertices)
      // brute force, for now:
      vertex_descriptor v_max = boost::graph_traits<Mesh>::null_vertex();
      face_descriptor f_max = boost::graph_traits<Mesh>::null_face();
      degree_size_type max_f_degree = 0;
      for (vertex_descriptor v : vertices(sm)) {
        if (degree(v, sm) <= 3) {
          continue;
        }

        for (halfedge_descriptor h : CGAL::halfedges_around_target(v, sm)) {
          face_descriptor f = face(h, sm);
          if (face_status[f] != Facet_status::DEFAULT) {
            continue;
          }
          if (degree(f, sm) > max_f_degree) {
            max_f_degree = degree(f, sm);
            f_max = f;
            v_max = v;
          }
        }
      }

      // once everything incident to high degree vertices is fixed or not fixable, let's end things for now
      if (f_max == boost::graph_traits<Mesh>::null_face()) {
        break;
      }

      // We can fix the maximal face as long as this would not create more than 3 fixed points
      std::cout << "iter " << iter << std::endl;
      std::cout << "  deal with " << f_max << " (" << degree(f_max, sm) << ")" << std::endl;

      // triangles are not very interesting
      if (CGAL::is_triangle(halfedge(f_max, sm), sm)) {
        face_status[f_max] = Facet_status::FIXED;
        continue;
      }

      // check how many high degree vertices already have 2 incident fixed faces
      std::cout << "  ";
      unsigned int c2_vertices = 0;
      for (vertex_descriptor iv : CGAL::vertices_around_face(halfedge(f_max, sm), sm)) {
          if (degree(iv, sm) > 3) {
            unsigned int fc = vertex_constraints(iv);
            std::cout << iv << "(" << fc << ") ";
            CGAL_assertion(fc <= 3);
            if (fc >= 2) {
              ++c2_vertices;
            }
          }
      }
      std::cout << std::endl;

      // we can't guarantee that we will be able to compute a plane through these 3 fixed planes
      if (c2_vertices > 3) {
        face_status[f_max] = Facet_status::NOT_FIXABLE;
      }

      face_status[f_max] = Facet_status::FIXED;
    }

    std::cout << "Statuses:" << std::endl;
    unsigned int fixed_faces_n = 0, not_fixable_faces_n = 0, unexplored_n = 0;
    for (face_descriptor f : faces(sm)) {
      if (face_status[f] == Facet_status::FIXED) {
        ++fixed_faces_n;
      } else if(face_status[f] == Facet_status::NOT_FIXABLE) {
        ++not_fixable_faces_n;
      } else {
        ++unexplored_n;
      }
    }

    std::cout << fixed_faces_n << " fixed faces" << std::endl;
    std::cout << not_fixable_faces_n << " not fixable faces" << std::endl;
    std::cout << unexplored_n << " unexplored faces" << std::endl;
  };

  flood(sm);

  CGAL::IO::write_polygon_mesh("results/out_constrained_polygonal_faces.ply",
                               sm, CGAL::parameters::stream_precision(17).face_color_map(color_map));
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // remesh_with_flip(argc, argv);
  // remesh_with_insertion(argc, argv);
  minimal_triangulation(argc, argv);
  // remesh_with_constrained_polygonal_faces(argc, argv);

  std::cout << "Done" << std::endl;

  return 0;
}
