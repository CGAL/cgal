#include <CGAL/config.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Default.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/io.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/iterator.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_3/parameters.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Sizing_field_with_aabb_tree.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/tags.h>

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <thread>
#include <utility>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Surface_mesh = CGAL::Surface_mesh<K::Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K, Surface_mesh>;


#ifdef CGAL_LINKED_WITH_TBB
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index>;

// Criteria
using Features_sizing_field = CGAL::Sizing_field_with_aabb_tree<K, Mesh_domain>;
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

std::atomic<bool> stop_meshing{false};
void signal_handler(int signal)
{
  if(signal == SIGINT) {
    stop_meshing.store(true, std::memory_order_release);
  }
}

int main(int argc, char*argv[])
{
  std::cout.precision(17);
  std::cerr.precision(17);
  std::clog.precision(17);
  const std::string input_file_name = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fandisk.off");
  const std::string output_file_name = (argc > 2) ? argv[2] : "out-tetmesh.vtu";
  Surface_mesh surface_mesh;
  if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(input_file_name, surface_mesh)) {
    std::cerr << "Error: Cannot read file " << input_file_name << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Read surface mesh from file " << input_file_name << std::endl;
    std::cout << "  - number of vertices: " << num_vertices(surface_mesh) << std::endl;
    std::cout << "  - number of faces:    " << num_faces(surface_mesh) << std::endl;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh)) {
    std::cerr << "ERROR: the input surface mesh is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  const bool is_closed = CGAL::is_closed(surface_mesh);
  if(is_closed) {
    std::cout << "  - the input surface mesh is closed." << std::endl;
  } else {
    std::cout << "  - WARNING: the input surface mesh is not closed." << std::endl;
  }

  // check if the input surface mesh is manifold
  std::size_t nb_of_non_manifold_vertices{0};
  CGAL::Counting_output_iterator count_it{&nb_of_non_manifold_vertices};
  count_it = CGAL::Polygon_mesh_processing::non_manifold_vertices(surface_mesh, count_it);
  if(nb_of_non_manifold_vertices > 0) {
    std::cout << "  - WARNING: the input surface mesh has " << nb_of_non_manifold_vertices
              << " non-manifold vertices."<< std::endl;
  }
  const bool self_intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(surface_mesh);
  if(self_intersecting) {
    std::cout << "  - WARNING: the input surface mesh is self-intersecting." << std::endl;
  }

  const bool manifold = !self_intersecting && nb_of_non_manifold_vertices == 0;

  if(manifold && is_closed) {
    std::cout << "  - the input surface mesh is manifold." << std::endl;
  }

  // Create domain
  Mesh_domain domain(std::move(surface_mesh));

  // Get sharp features
  domain.detect_features();

  namespace params = CGAL::parameters;

  // Mesh criteria
  auto bbox = domain.bbox();
  const double diagonal_length =
      CGAL::sqrt(CGAL::square(bbox.x_span()) + CGAL::square(bbox.y_span()) + CGAL::square(bbox.z_span()));
  const double size_bound = diagonal_length * 0.05;
  const Features_sizing_field edge_sizing_field(0.07, domain);
  const Mesh_criteria criteria(params::edge_size(edge_sizing_field)
                                      .edge_min_size(size_bound / 1.e6)
                                      .facet_angle(25)
                                      .facet_size(size_bound)
                                      .facet_distance(size_bound * 0.1)
                                      .cell_radius_edge_ratio(3)
                                      .cell_size(size_bound));

  const auto manifold_criteria = manifold ? (is_closed ? params::manifold()
                                                       : params::manifold_with_boundary())
                                          : (params::non_manifold());

  auto mesh_options = params::mesh_3_options(params::nonlinear_growth_of_balls(true)
                                                    .pointer_to_stop_atomic_boolean(&stop_meshing));
  // Mesh generation
  C3t3 c3t3;

  std::cout << "Start meshing (" << (Concurrency_tag::is_parallel ? "parallel" : "sequential")
            << ")..." << std::endl;

  std::atomic<std::size_t> nb_of_protecting_balls{0};
  std::atomic<bool> meshing_finished{false};

  namespace chr = std::chrono;
  auto start_time = chr::system_clock::now();

  auto monitoring_task = [&] {
    std::size_t nb_of_vertices = 0;
    bool protection_done = false;
    while(!meshing_finished.load(std::memory_order_acquire))
    {
      std::this_thread::sleep_for(chr::seconds(1));

      if(!protection_done) {
        auto nb = nb_of_protecting_balls.load(std::memory_order_acquire);
        if(nb > 0) {
          std::cout << "Protection phase completed (" << nb << " protecting balls). Starting meshing phase..." << std::endl;
          protection_done = true;
        }
      }
      auto current_time = chr::system_clock::now();
      auto elapsed_time = chr::duration_cast<chr::seconds>(current_time - start_time).count();
      std::size_t current_nb_of_vertices = c3t3.triangulation().number_of_vertices();
      std::cout << "  " <<  std::to_string(elapsed_time) + " seconds, ";
      if(stop_meshing.load(std::memory_order_acquire)) {
        std::cout << "interrupted, ";
      } else if(protection_done) {
        std::cout << "meshing,     ";
      } else {
        std::cout << "protection,  ";
      }
      std::cout << "number of vertices: " << current_nb_of_vertices
                << "   (" << (current_nb_of_vertices - nb_of_vertices) << " new vertices)"
                << std::endl;
      nb_of_vertices = current_nb_of_vertices;
    }
  };

  std::signal(SIGINT, signal_handler);

  // Start monitoring thread
  std::thread monitoring_thread(monitoring_task);

  // Run meshing on main thread
  // c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, manifold_criteria, mesh_options);
  CGAL::Mesh_3::internal::C3t3_initializer<C3t3, Mesh_domain, Mesh_criteria> mesh_initializer;
  mesh_initializer(c3t3, domain, criteria, true, mesh_options.v);

  if(!stop_meshing) {
    nb_of_protecting_balls.store(c3t3.triangulation().number_of_vertices(), std::memory_order_release);
    CGAL::refine_mesh_3(c3t3, domain, criteria, manifold_criteria, mesh_options);
  }

  // Signal monitoring thread to stop and wait for it
  meshing_finished.store(true, std::memory_order_release);
  monitoring_thread.join();

  int exit_code = EXIT_SUCCESS;
  if(stop_meshing) {
    std::cout << "ERROR: Meshing interrupted.\n";
    exit_code = EXIT_FAILURE;
  } else {
    std::cout << "Meshing completed.\n";
  }
  auto end_time = chr::system_clock::now();
  auto elapsed_time = chr::duration_cast<chr::seconds>(end_time - start_time).count();
  std::cout << "  elapsed time:             " << static_cast<double>(elapsed_time) << " seconds." << std::endl;
  std::cout << "  number of vertices:       " << c3t3.triangulation().number_of_vertices() << std::endl;
  std::cout << "  number of surface facets: " << c3t3.number_of_facets() << std::endl;
  std::cout << "  number of cells:          " << c3t3.number_of_cells() << std::endl;

  // Output
  std::ofstream ofs(output_file_name);
  ofs.precision(17);
  CGAL::IO::output_to_vtu(ofs, c3t3, CGAL::IO::ASCII);
  if(ofs.fail()) {
    std::cerr << "Error: Cannot write the surface mesh to file " << output_file_name << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Tet mesh written to file " << output_file_name << std::endl;
  }

  return exit_code;
}
