#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_placement.h>

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;

  const char* filename = argv[1];
  std::ifstream is(argv[1]);

  if(!is)
  {
    std::cerr << "Filename provided is invalid\n";
    return EXIT_FAILURE;
  }

  is >> surface_mesh;
  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh has " << num_vertices(surface_mesh) << " nv "
                                 << num_edges(surface_mesh) << " ne "
                                 << num_faces(surface_mesh) << " nf" << std::endl;

  const double stop_threshold = (argc > 2) ? std::stod(argv[2]) : 0.1;

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_threshold);

  // Garland&Heckbert simplification maintains an error matrix at each vertex,
  // which must be accessible for the cost and placement evaluations.
  typedef typename SMS::GarlandHeckbert_cost_matrix<Kernel>::type               Cost_matrix;
  typedef CGAL::dynamic_vertex_property_t<Cost_matrix>                          Cost_property;
  typedef typename boost::property_map<Surface_mesh, Cost_property>::type       Vertex_cost_map;
  Vertex_cost_map vcm = get(Cost_property(), surface_mesh);

  typedef SMS::GarlandHeckbert_cost<Vertex_cost_map>                            GH_cost;
  typedef SMS::GarlandHeckbert_placement<Vertex_cost_map>                       GH_placement;
  typedef SMS::Bounded_normal_change_placement<GH_placement>                    Bounded_GH_placement;

  GH_cost cost(vcm);
  GH_placement gh_placement(vcm);
  Bounded_GH_placement placement(gh_placement);

  int r = SMS::edge_collapse(surface_mesh, stop, CGAL::parameters::get_cost(cost)
                                                                  .get_placement(placement));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "Time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "\nFinished...\n" << r << " edges removed.\n" << surface_mesh.number_of_edges() << " final edges.\n";

  std::ofstream os(argc > 3 ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
