#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <CGAL/boost/graph/partition.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <tbb/task_group.h>

#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>                                      Kernel;

typedef Kernel::Point_3                                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                                         Triangle_mesh;


namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

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
  double edge_length;

  Simplify()
    :tm_ptr(NULL), edge_length(0)
  {}

  Simplify(Triangle_mesh& tm, double edge_length)
    : tm_ptr(&tm), edge_length(edge_length)
  {}

  void operator()() const
  {

    SMS::Edge_length_cost<Triangle_mesh> cost;
    SMS::Edge_length_stop_predicate<double> stop(edge_length);
    Border_is_constrained_edge_map eicm(*tm_ptr);
    SMS::Constrained_placement<SMS::LindstromTurk_placement<Triangle_mesh>,
                               Border_is_constrained_edge_map > placement(eicm);

    SMS::edge_collapse(*tm_ptr, stop,
                       CGAL::parameters::edge_is_constrained_map(eicm)
                        .get_placement(placement)
                        .get_cost(cost));
  }
};

int main(int argc, char** argv)
{
  typedef boost::graph_traits<Triangle_mesh>::face_descriptor face_descriptor;

  std::ifstream in((argc>1) ? std::string(argv[1]) : CGAL::data_file_path("meshes/elephant.off"));
  double edge_length = (argc>2) ? std::atof(argv[2]) : 0.02;
  int number_of_parts = (argc>3) ? std::atoi(argv[3]) : 8;

  Triangle_mesh tm;
  in >> tm;

  std::cerr << "Input: #V = "<< num_vertices(tm) << " #E = "<< num_edges(tm)
            << " #F = " << num_faces(tm) << std::endl;

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

  SMS::LindstromTurk_placement<Triangle_mesh> placement;
  SMS::Edge_length_cost<Triangle_mesh> cost;
  SMS::Edge_length_stop_predicate<double> stop(edge_length);

  CGAL::Face_filtered_graph<Triangle_mesh> fg(tm);
  std::vector<Triangle_mesh> meshes(number_of_parts);
  for (int i=0; i< number_of_parts; ++i)
  {
    fg.set_selected_faces(i, fpmap);
    CGAL::copy_face_graph(fg, meshes[i]);
  }

  tbb::task_group tasks;
  for(int id=0; id<number_of_parts; ++id)
  {
    tasks.run(Simplify(meshes[id], edge_length));
  }
  tasks.wait();

  // reassemble the mesh
  Triangle_mesh final_mesh;
  std::size_t nv=0, nf=0, ne=0;
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
  SMS::edge_collapse(final_mesh, stop,
                     CGAL::parameters::get_placement(placement)
                     .get_cost(cost));

  final_mesh.collect_garbage();
  t.stop();

  std::cerr << "  Total time: " << t.time() << " sec.\n"
          << "  Result: #V = " << num_vertices(final_mesh)
          << " #E = " << num_edges(final_mesh)
          << " #F = " << num_faces(final_mesh) << std::endl;

  std::ofstream("simplified_mesh.off") << std::setprecision(17) << final_mesh;

  return EXIT_SUCCESS;
}
