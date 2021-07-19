#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_probabilistic_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_probabilistic_tri_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_triangle_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <boost/histogram.hpp>
#include <boost/format.hpp>
#include <CGAL/Polygon_mesh_processing/measure.h>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

typedef typename boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

using hist = boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double, 
      boost::use_default, boost::use_default, boost::use_default>>>; 

// taken from a boost example
template<typename histo>
void print_hist(histo h)
{
  using namespace boost::histogram;
  
  std::ostringstream os;
  for (auto&& x : indexed(h, coverage::all)) 
  {
    // alternative output formatting
    //os << boost::format("bin %2i [%4.1f, %4.1f): %i\n")
    //      % x.index() % x.bin().lower() % x.bin().upper() % *x;
    os << boost::format("%i, ")
          % *x;
  }

  std::cout << os.str() << std::endl;
}

hist generate_statistics(Surface_mesh mesh, hist histo, std::function<FT(edge_descriptor, 
      const Surface_mesh&)> f)
{
  for (auto e : mesh.edges())
  {
    FT value = f(e, mesh); //PMP::edge_length(mesh.halfedge(e), mesh); 
    histo(value);
  }

  return histo;
}

template<typename Policy>
hist generate_statistics(Surface_mesh mesh, Policy p, 
    hist histo, std::function<FT(edge_descriptor, const Surface_mesh&)> f)
{
  typedef typename Policy::Get_cost Cost; 
  typedef typename Policy::Get_placement Placement;
  
  typedef SMS::Bounded_normal_change_placement<Placement> Bounded_placement;
  
  const Cost& cost = p.get_cost();
  const Placement& unbounded_placement = p.get_placement();

  Bounded_placement placement(unbounded_placement);

  //TODO don't hardcode ratio 
  const double ratio = 0.2;
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);
  
  // collapse edges ignoring result code
  SMS::edge_collapse(mesh, stop,
      CGAL::parameters::get_cost(cost).get_placement(placement));
  
  return generate_statistics(mesh, histo, f);
}

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/cube-meshed.off";
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  using namespace boost::histogram;
  hist edge_histo = make_histogram(axis::regular<>(14, 0.1, 4.0, "length"));
  
  auto edge_length = [] (edge_descriptor e, const Surface_mesh& mesh) { 
    return PMP::edge_length(mesh.halfedge(e), mesh);
  };
  
  std::cout << "Original mesh histogram:\n";
  print_hist(generate_statistics(surface_mesh, edge_histo, edge_length));

  std::cout << "Decimated probabilistic mesh histogram:\n";
  print_hist(generate_statistics<SMS::GarlandHeckbert_probabilistic_policies<Surface_mesh, Kernel>>
      (surface_mesh, {surface_mesh, 100}, edge_histo, edge_length));

  std::cout << "Decimated classic mesh histogram:\n";
  print_hist(generate_statistics<SMS::GarlandHeckbert_policies<Surface_mesh, Kernel>>
      (surface_mesh, {surface_mesh, 100}, edge_histo, edge_length));

  std::cout << "Decimated classic tri mesh histogram:\n";
  print_hist(generate_statistics<SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel>>
      (surface_mesh, {surface_mesh, 100}, edge_histo, edge_length));

  return EXIT_SUCCESS;
}
