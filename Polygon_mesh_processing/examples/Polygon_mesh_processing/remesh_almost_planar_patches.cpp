#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>

#include <boost/property_map/vector_property_map.hpp>

#include <iostream>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main()
{
  Surface_mesh sm;
  CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/fandisk.off"), sm);

  //apply a perturbation to input vertices so that points are no longer coplanar
  PMP::random_perturbation(sm, 0.001);

  // declare vectors to store mesh properties
  std::vector<std::size_t> region_ids(num_faces(sm));
  std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
  std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions
  boost::vector_property_map<CGAL::Epick::Vector_3> normal_map; // normal of the supporting planes of the regions detected

  // detect planar regions in the mesh
  std::size_t nb_regions =
    PMP::region_growing_of_planes_on_faces(sm,
                                           CGAL::make_random_access_property_map(region_ids),
                                           CGAL::parameters::cosine_of_maximum_angle(0.98).
                                                             region_primitive_map(normal_map).
                                                             maximum_distance(0.011));

  // detect corner vertices on the boundary of planar regions
  std::size_t nb_corners =
    PMP::detect_corners_of_regions(sm,
                                   CGAL::make_random_access_property_map(region_ids),
                                   nb_regions,
                                   CGAL::make_random_access_property_map(corner_id_map),
                                   CGAL::parameters::cosine_of_maximum_angle(0.98).
                                                     maximum_distance(0.011).
                                                     edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

  // run the remeshing algorithm using filled properties
  Surface_mesh out;
  PMP::remesh_almost_planar_patches(sm,
                                    out,
                                    nb_regions, nb_corners,
                                    CGAL::make_random_access_property_map(region_ids),
                                    CGAL::make_random_access_property_map(corner_id_map),
                                    CGAL::make_random_access_property_map(ecm),
                                    CGAL::parameters::patch_normal_map(normal_map));

  CGAL::IO::write_polygon_mesh("fandisk_remeshed.off", out);

  return 0;
}
