#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_probabilistic_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

#define CGAL_USE_BASIC_VIEWER

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;
typedef CGAL::Bbox_3 Bbox_3;

namespace SMS = CGAL::Surface_mesh_simplification;

template<typename Cost, typename Placement>
Surface_mesh collapse(Surface_mesh surface_mesh, std::string outname, 
    FT ratio, Cost cost, Placement placement)
{
  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return { } ; 
  }

  std::cout << "Input mesh has " << num_vertices(surface_mesh) << " nv "
                                 << num_edges(surface_mesh) << " ne "
                                 << num_faces(surface_mesh) << " nf" << std::endl;


  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
  
  std::cout << "Aiming for " << 100 * ratio << "% of the input edges..." << std::endl;
  
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);
  
  int r = SMS::edge_collapse(surface_mesh, stop,
                             CGAL::parameters::get_cost(cost).get_placement(placement));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "Time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "\nFinished!\n" << r 
            << " edges removed.\n" << surface_mesh.number_of_edges() 
            << " final edges.\n";

  CGAL::IO::write_polygon_mesh(outname, surface_mesh, CGAL::parameters::stream_precision(17));
  
  return surface_mesh;
}

void normalize(Surface_mesh& mesh) 
{
  Bbox_3 bbox { };

  for (auto v : mesh.vertices())
  {
    bbox += mesh.point(v).bbox();
  }

  FT max = std::max(std::max(bbox.xmax(), bbox.ymax()), bbox.zmax());
  
  for (auto v : mesh.vertices())
  {
    Point_3& p = mesh.point(v);
    p = {p.x() / max, p.y() / max, p.z() / max};
  }
}

// takes 4 arguments: input mesh, target reduction percent, name of output file (will be formated
// as name_classic, name_prob_1 etc.) and starting variance
//
// then do classic decimation (with bounded placement) and 5 iterations of probabilistic placement
// where we halve the variance each step
int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  
  // command line arguments
  const char* filename = (argc > 1) ? argv[1] : "data/cube-meshed.off";
  const double ratio = (argc > 2) ? std::stod(argv[2]) : 0.2;
  const std::string outname = (argc > 3) ? argv[3] : "test";
  
  FT variance = (argc > 4) ? std::stod(argv[4]) : 0.05;
  
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  // Garland&Heckbert simplification policies
  typedef typename SMS::GarlandHeckbert_policies<Surface_mesh, Kernel> Classic_policies;
  typedef typename SMS::GarlandHeckbert_probabilistic_policies<Surface_mesh, Kernel> Prob_policies;
  typedef typename Classic_policies::Get_cost Classic_cost;
  typedef typename Classic_policies::Get_placement Classic_placement_unbounded;
  typedef typename CGAL::Surface_mesh_simplification::Bounded_normal_change_placement<
    Classic_placement_unbounded> Classic_placement;
  typedef typename Prob_policies::Get_cost Prob_cost;
  typedef typename Prob_policies::Get_placement Prob_placement;
  
  Classic_policies classic_policies(surface_mesh, 100);
  const Classic_cost& classic_cost = classic_policies.get_cost();
  const Classic_placement& classic_placement { classic_policies.get_placement() };
  
  
  normalize(surface_mesh);
  
  collapse(surface_mesh, outname + "_classic.off", ratio, classic_cost, classic_placement);

  for (int i = 0; i < 5; ++i)
  {
    Prob_policies prob_policies(surface_mesh, 100, variance, variance);
    
    const Prob_cost& prob_cost = prob_policies.get_cost();
    const Prob_placement& prob_placement = prob_policies.get_placement();
    
    std::string filename = outname + "_prob_" + std::to_string(i) + ".off";

    std::cout << "Using variance of " << variance << " for both parameters\n";
    collapse(surface_mesh, filename, ratio, prob_cost, prob_placement);

    variance = variance / 2;
  }
  
  return EXIT_SUCCESS;
}
