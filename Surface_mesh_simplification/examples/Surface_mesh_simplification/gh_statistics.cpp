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
#include <boost/histogram.hpp>
#include <boost/format.hpp>
#include <CGAL/Polygon_mesh_processing/measure.h>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

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

hist edge_length_histo(Surface_mesh mesh)
{
  using namespace boost::histogram;
  hist histo = make_histogram(axis::regular<>(14, 0.1, 4.0, "length"));
  
  for (auto e : mesh.edges())
  {
    FT length = PMP::edge_length(mesh.halfedge(e), mesh);
    histo(length);
  }
  
  return histo;
}

template<typename Policy>
hist edge_length_histo(Surface_mesh mesh, Policy p)
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
  
  return edge_length_histo(mesh);
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

  std::cout << "Original mesh histogram:\n";
  print_hist(edge_length_histo(surface_mesh));
  
  std::cout << "Decimated probabilistic mesh histogram:\n";
  print_hist(edge_length_histo<SMS::GarlandHeckbert_probabilistic_policies<Surface_mesh, Kernel>>
    (surface_mesh, {surface_mesh, 100}));
  
  std::cout << "Decimated probabilistic tri mesh histogram:\n";
  print_hist(edge_length_histo<SMS::GarlandHeckbert_probabilistic_tri_policies<Surface_mesh, Kernel>>
    (surface_mesh, {surface_mesh, 100}));

  std::cout << "Decimated classic mesh histogram:\n";
  print_hist(edge_length_histo<SMS::GarlandHeckbert_policies<Surface_mesh, Kernel>>
    (surface_mesh, {surface_mesh, 100}));
  
  std::cout << "Decimated classic tri mesh histogram:\n";
  print_hist(edge_length_histo<SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel>>
    (surface_mesh, {surface_mesh, 100}));
  
  return EXIT_SUCCESS;
}
