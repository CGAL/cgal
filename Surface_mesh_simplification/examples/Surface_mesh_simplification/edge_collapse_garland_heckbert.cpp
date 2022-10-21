#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

typedef SMS::GarlandHeckbert_plane_policies<Surface_mesh, Kernel>                  Classic_plane;
typedef SMS::GarlandHeckbert_probabilistic_plane_policies<Surface_mesh, Kernel>    Prob_plane;
typedef SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel>               Classic_tri;
typedef SMS::GarlandHeckbert_probabilistic_triangle_policies<Surface_mesh, Kernel> Prob_tri;

template <typename GHPolicies>
void collapse_gh(Surface_mesh& mesh,
                 const double ratio)
{
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  // Garland&Heckbert simplification policies

  typedef typename GHPolicies::Get_cost                                        GH_cost;
  typedef typename GHPolicies::Get_placement                                   GH_placement;
  typedef SMS::Bounded_normal_change_placement<GH_placement>                    Bounded_GH_placement;

  GHPolicies gh_policies(mesh);
  const GH_cost& gh_cost = gh_policies.get_cost();
  const GH_placement& gh_placement = gh_policies.get_placement();
  Bounded_GH_placement placement(gh_placement);

  int r = SMS::edge_collapse(mesh, stop,
                             CGAL::parameters::get_cost(gh_cost)
                                              .get_placement(placement));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "Time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "\nFinished!\n" << r << " edges removed.\n" << edges(mesh).size() << " final edges.\n";
}

// Usage:
// ./command [input] [ratio] [policy] [outpout]
// policy can be "cp" (classic plane), "ct" (classic triangle), "pp" (probabilistic plane), "pt" (probabilistic triangle)
int main(int argc, char** argv)
{
  Surface_mesh mesh;
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube-meshed.off");

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

  const double ratio = (argc > 2) ? std::stod(argv[2]) : 0.2;
  std::cout << "Collapsing edges of mesh: " << filename << ", aiming for " << 100 * ratio << "% of the input edges..." << std::endl;

  const std::string policy = (argc > 3) ? argv[3] : "cp";
  if(policy == "cp")
    collapse_gh<Classic_plane>(mesh, ratio);
  else if(policy == "ct")
    collapse_gh<Classic_tri>(mesh, ratio);
  else if(policy == "pp")
    collapse_gh<Prob_plane>(mesh, ratio);
  else
    collapse_gh<Prob_tri>(mesh, ratio);

  CGAL::IO::write_polygon_mesh((argc > 4) ? argv[4] : "out.off", mesh, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
