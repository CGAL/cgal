

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>

#include <boost/property_map/vector_property_map.hpp>

#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Plane_3 = K::Plane_3;
using Vector_3 = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

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

  // run the remeshing algorithm using filled properties
  std::vector<Vector_3> normal_map(num_faces(sm));
  for(face_descriptor f:faces(sm))
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

  CGAL::IO::write_polygon_mesh("out_remeshed.off", out, CGAL::parameters::stream_precision(17));

  // check the degree of each vertex
  for(vertex_descriptor v : vertices(out))
  {
    // @todo this should count incident face _ID_ and not incident faces:
    // incident to 'v' could be a face that is not simply connected, but
    // if all the vertices of this face are still of degree 3 (with IDs)
    // then we can still tilt the face
    std::size_t d = CGAL::halfedges_around_target(v, out).size();
    if(d != 3) {
      std::cout << out.point(v) << " is of degree " << d << "\n";
    }
  }

  // perturbate planes
  auto nudge = [](double v)
  {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> rdist(1e-10, 1e-9);

    double step = rdist(gen);
    double nv = v + step;
    return nv;
  };

  // @todo with incident face _IDS_ instead of incident faces,
  // maybe we can avoid having to triangulate EVERYTHING as soon as one vertex
  // does not have degree 3 by marking the vertex that needs to be moved manually
  // and tilts are possible for the incident faces as long as they have less than 3 vertices
  // that need to be manually moved (if 3 out of n vertices must be moved manually, we can still
  // tilt the face such that it matches the 3 moved vertices and other vertices will move
  // with the tilted plane)
  for(face_descriptor f : faces(out))
  {
    // @todo see other todo
    double na = nudge(plane_map[f].a()),
           nb = nudge(plane_map[f].b()),
           nc = nudge(plane_map[f].c()),
           d = plane_map[f].d(); // pointless to nudge d
    double n = std::sqrt(CGAL::square(na) + CGAL::square(nb) + CGAL::square(nc));

    // @todo if it's zero, nudge differently
    // @todo also, it shouldn't invert the face (bounded normal change)
    CGAL_assertion(n != 0);

    plane_map[f] = Plane_3(na/n, nb/n, nc/n, d/n);
  }

  // now recompute the coordinates of the vertices
  for(vertex_descriptor v : vertices(out))
  {
    std::array<Plane_3, 3> planes;
    halfedge_descriptor h = halfedge(v, out);
    for(int i=0; i<3;++i)
    {
      planes[i] = plane_map[face(h, out)];
      h = opposite(next(h, out), out);
    }

    auto inter = CGAL::intersection(planes[0], planes[1], planes[2]);
    if(inter.has_value())
    {
      Point_3 p;
      if(!CGAL::assign(p, inter.value()))
      {
        std::cout << "Could not compute new coordinates for " << v << "(" << out.point(v) << ") (intersection is not a point)\n";
      }
      else
      {
        CGAL_warning(p != out.point(v));
        put(CGAL::vertex_point, out, v, p);
      }
    }
    else
    {
      std::cout << "Could not compute new coordinates for " << v << "(" << out.point(v) << ") (intersection is empty)\n";
    }
  }

  CGAL::IO::write_polygon_mesh("out_remeshed_perturbed.off", out, CGAL::parameters::stream_precision(17));

  return 0;
}
