#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <CGAL/boost/graph/partition.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <tbb/task_group.h>

#include <iostream>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;

typedef Kernel::Point_3                                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                                         Triangle_mesh;


namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef SMS::GarlandHeckbert_plane_and_line_policies<Triangle_mesh, Kernel> Plane_and_line;
typedef Plane_and_line::Get_cost                                            GH_cost;
typedef Plane_and_line::Get_placement                                       GH_placement;
typedef SMS::Bounded_normal_change_filter<>                                 Filter;

struct Border_is_constrained_edge_map{
  const Triangle_mesh& sm;
  typedef Triangle_mesh::Edge_index key_type;
  typedef bool value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  Border_is_constrained_edge_map(const Triangle_mesh& sm)
    : sm(sm)
  {}
  friend bool get(Border_is_constrained_edge_map m, const key_type& edge) {
    return  CGAL::is_border(edge, m.sm);
  }
};

struct Simplify
{
  Triangle_mesh* tm_ptr;
  std::size_t thread_expected_nb_edges;

  Simplify()
    :tm_ptr(NULL), thread_expected_nb_edges(0)
  {}

  Simplify(Triangle_mesh& tm, double ratio)
    : tm_ptr(&tm)
  {
    std::size_t nb_edges = num_edges(tm);
    thread_expected_nb_edges=nb_edges*ratio*1.01; // we take more on purpose to let the last pass achieve expected nb
  }

  void operator()() const
  {
    SMS::Edge_count_stop_predicate<Triangle_mesh> stop(thread_expected_nb_edges);

    Plane_and_line gh_policies(*tm_ptr);
    const GH_cost& gh_cost = gh_policies.get_cost();
    const GH_placement& gh_placement = gh_policies.get_placement();
    Filter filter;
    Border_is_constrained_edge_map eicm(*tm_ptr);

    SMS::edge_collapse(*tm_ptr, stop,
                       CGAL::parameters::get_cost(gh_cost)
                                        .edge_is_constrained_map(eicm)
                                        .filter(filter)
                                        .get_placement(gh_placement));
  }
};

int main(int argc, char** argv)
{
  typedef boost::graph_traits<Triangle_mesh>::face_descriptor face_descriptor;

  std::ifstream in((argc>1) ? std::string(argv[1]) : CGAL::data_file_path("meshes/refined_elephant.off"));
  double ratio = (argc > 2) ? std::stod(argv[2]) : 0.2;
  int number_of_parts = (argc>3) ? std::atoi(argv[3]) : 8;

  Triangle_mesh tm;
  in >> tm;

  std::cerr << "Input: #V = "<< num_vertices(tm) << " #E = "<< num_edges(tm)
            << " #F = " << num_faces(tm) << std::endl;

  std::size_t nb_in_edges=num_edges(tm);
  std::size_t expected_nb_edges=nb_in_edges*ratio;

  std::cerr << "expected # output edges " << expected_nb_edges << "\n";

  // we give on purpose a higher target per thread to let the last pass achieve the expected number
  CGAL::Real_timer t;
  t.start();

  // Partition ID map
  Triangle_mesh::Property_map<face_descriptor, std::size_t> fpmap
    = tm.add_property_map<face_descriptor, std::size_t>("f:partition").first;

  // Set some custom options for METIS
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_CONTIG] = 1; // need to have contiguous subdomains
  options[METIS_OPTION_UFACTOR] = 1;

  CGAL::METIS::partition_dual_graph(tm, number_of_parts, CGAL::parameters::METIS_options(&options)
                                    // .vertex_index_map(get(boost::vertex_index,tm))
                                                                          .face_partition_id_map(fpmap));


  CGAL::Face_filtered_graph<Triangle_mesh> fg(tm);
  std::vector<Triangle_mesh> meshes(number_of_parts);
  for (int i=0; i< number_of_parts; ++i)
  {
    fg.set_selected_faces(i, fpmap);
    CGAL::copy_face_graph(fg, meshes[i]);
  }
  std::cout << "Partition computed in " << t.time() << "\n";


  tbb::task_group tasks;
  for(int id=0; id<number_of_parts; ++id)
  {
    tasks.run(Simplify(meshes[id], ratio));
  }
  tasks.wait();

  std::cout << "Parallel decimation done in " << t.time() << "\n";

  // reassemble the mesh
  Triangle_mesh final_mesh;
  Triangle_mesh::size_type nv = 0, nf = 0, ne = 0;
  for (int i=0; i< number_of_parts; ++i)
  {
    meshes[i].collect_garbage();
    nv += num_vertices(meshes[i]);
    nf += num_faces(meshes[i]);
    ne += num_edges(meshes[i]);
  }
  final_mesh.reserve(nv, ne, nf);
  for (int i=0; i< number_of_parts; ++i)
  {
    CGAL::copy_face_graph(meshes[i], final_mesh);
  }
  PMP::stitch_borders(final_mesh);

  // final pass (could be applied only on interface boundaries)
  SMS::Edge_count_stop_predicate<Triangle_mesh> stop(expected_nb_edges);

  Plane_and_line gh_policies(final_mesh);
  const GH_cost& gh_cost = gh_policies.get_cost();
  const GH_placement& gh_placement = gh_policies.get_placement();
  Filter filter;

  SMS::edge_collapse(final_mesh, stop,
                     CGAL::parameters::get_cost(gh_cost)
                                      .filter(filter)
                                      .get_placement(gh_placement));


  final_mesh.collect_garbage();
  t.stop();

  std::cerr << "  Total time: " << t.time() << " sec.\n"
          << "  Result: #V = " << num_vertices(final_mesh)
          << " #E = " << num_edges(final_mesh)
          << " #F = " << num_faces(final_mesh) << std::endl;

  std::ofstream("decimated_mesh_parallel.off") << std::setprecision(17) << final_mesh;

  return EXIT_SUCCESS;
}
