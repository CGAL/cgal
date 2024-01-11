#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_ratio_stop_predicate.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>

// deprecated
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/internal/helpers.h>

#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

typedef SMS::Edge_count_stop_predicate<Mesh> Edge_count_stop;
typedef SMS::Face_count_stop_predicate<Mesh> Face_count_stop;
typedef SMS::Edge_count_ratio_stop_predicate<Mesh> Edge_count_ratio_stop;
typedef SMS::Face_count_ratio_stop_predicate<Mesh> Face_count_ratio_stop;

typedef SMS::Count_stop_predicate<Mesh> Count_stop;
typedef SMS::Count_ratio_stop_predicate<Mesh> Count_ratio_stop;

typedef SMS::Count_ratio_stop_predicate<Mesh> Count_ratio_stop;

typedef SMS::Edge_length_cost<Mesh> Cost;
typedef SMS::Midpoint_placement<Mesh> Placement;

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube-meshed.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh has " << num_vertices(mesh) << " nv "
                                 << num_edges(mesh) << " ne "
                                 << num_faces(mesh) << " nf" << std::endl;

  Cost cost;
  Placement placement;

  // Edge_count_stop
  {
    Mesh mesh_cpy = mesh;
    const std::size_t expected_ne = num_edges(mesh_cpy) / 2;
    Edge_count_stop stop(expected_ne);
    SMS::edge_collapse(mesh_cpy, stop,
                       CGAL::parameters::get_cost(cost)
                                        .get_placement(placement));

    std::cout << "Output mesh has " << CGAL::internal::exact_num_vertices(mesh_cpy) << " nv "
                                    << CGAL::internal::exact_num_edges(mesh_cpy) << " ne "
                                    << CGAL::internal::exact_num_faces(mesh_cpy) << " nf" << std::endl;

    assert(CGAL::internal::exact_num_edges(mesh_cpy) < expected_ne);
  }

  // Count_stop
  {
    Mesh mesh_cpy = mesh;
    const std::size_t expected_ne = num_edges(mesh_cpy) + 1;
    Count_stop stop(expected_ne);
    SMS::edge_collapse(mesh_cpy, stop,
                       CGAL::parameters::get_cost(cost)
                                        .get_placement(placement));

    std::cout << "Output mesh has " << CGAL::internal::exact_num_vertices(mesh_cpy) << " nv "
                                    << CGAL::internal::exact_num_edges(mesh_cpy) << " ne "
                                    << CGAL::internal::exact_num_faces(mesh_cpy) << " nf" << std::endl;

    assert(CGAL::internal::exact_num_edges(mesh_cpy) < expected_ne);
  }

  // Face_count_stop
  {
    Mesh mesh_cpy = mesh;
    const std::size_t expected_nf = num_faces(mesh_cpy) / 4;
    Face_count_stop stop(expected_nf);
    SMS::edge_collapse(mesh_cpy, stop,
                       CGAL::parameters::get_cost(cost)
                                        .get_placement(placement));

    std::cout << "Output mesh has " << CGAL::internal::exact_num_vertices(mesh_cpy) << " nv "
                                    << CGAL::internal::exact_num_edges(mesh_cpy) << " ne "
                                    << CGAL::internal::exact_num_faces(mesh_cpy) << " nf" << std::endl;

    assert(CGAL::internal::exact_num_faces(mesh_cpy) < expected_nf);
  }

  /// RATIO

  // Edge_count_ratio_stop
  {
    Mesh mesh_cpy = mesh;
    const double ratio = 0.5;
    const std::size_t initial_ne = num_edges(mesh_cpy);
    Edge_count_ratio_stop stop(ratio);
    SMS::edge_collapse(mesh_cpy, stop,
                       CGAL::parameters::get_cost(cost)
                                        .get_placement(placement));

    std::cout << "Output mesh has " << CGAL::internal::exact_num_vertices(mesh_cpy) << " nv "
                                    << CGAL::internal::exact_num_edges(mesh_cpy) << " ne "
                                    << CGAL::internal::exact_num_faces(mesh_cpy) << " nf" << std::endl;

    assert(CGAL::internal::exact_num_edges(mesh_cpy) / initial_ne < ratio);
  }

  // Count_ratio_stop
  {
    Mesh mesh_cpy = mesh;
    const double ratio = 1.;
    const std::size_t initial_ne = num_edges(mesh_cpy);
    Count_ratio_stop stop(ratio);
    SMS::edge_collapse(mesh_cpy, stop,
                       CGAL::parameters::get_cost(cost)
                                        .get_placement(placement));

    std::cout << "Output mesh has " << CGAL::internal::exact_num_vertices(mesh_cpy) << " nv "
                                    << CGAL::internal::exact_num_edges(mesh_cpy) << " ne "
                                    << CGAL::internal::exact_num_faces(mesh_cpy) << " nf" << std::endl;

    assert(CGAL::internal::exact_num_edges(mesh_cpy) / initial_ne < ratio);
  }

  // Face_count_ratio_stop
  {
    Mesh mesh_cpy = mesh;
    const double ratio = 0.7;
    const std::size_t initial_nf = num_faces(mesh_cpy);
    Face_count_ratio_stop stop(ratio, mesh_cpy);
    SMS::edge_collapse(mesh_cpy, stop,
                       CGAL::parameters::get_cost(cost)
                                        .get_placement(placement));

    std::cout << "Output mesh has " << CGAL::internal::exact_num_vertices(mesh_cpy) << " nv "
                                    << CGAL::internal::exact_num_edges(mesh_cpy) << " ne "
                                    << CGAL::internal::exact_num_faces(mesh_cpy) << " nf" << std::endl;

    assert(CGAL::internal::exact_num_faces(mesh_cpy) / initial_nf < ratio);
  }

  return 0;
}

