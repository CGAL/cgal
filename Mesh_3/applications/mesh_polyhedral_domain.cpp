#include "debug.h"
#include <CGAL/AABB_intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Domain 
// (we use exact intersection computation with Robust_intersection_traits_3)
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
typedef CGAL::Polyhedron_3<Geom_traits> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Geom_traits> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

template <typename T>
T set_arg(const std::string& param_name,
          const std::string& param_string,
          const po::variables_map& vm)
{
  if ( vm.count(param_name) )
  {
    T param_value = vm[param_name].as<T>();
    std::cout << param_string << ": " << param_value << "\n";
    return param_value;
  }
  else
  {
    std::cout << param_string << " ignored.\n";
    return T();
  }
}

int main(int argc, char** argv)
{
  po::options_description generic("Generic options");
  generic.add_options() ("help", "Produce help message");
  generic.add_options()("file", po::value<std::string>(), "Mesh polyhedron contained in that file");
  
  po::options_description mesh("Mesh generation parameters");
  mesh.add_options()("facet_angle", po::value<double>(), "Set facet angle bound")
  ("facet_size", po::value<double>(), "Set facet size bound")
  ("facet_error", po::value<double>(), "Set facet approximation error bound")
  ("tet_shape", po::value<double>(), "Set tet radius-edge bound")
  ("tet_size", po::value<double>(), "Set tet size bound");
  
  po::options_description desc("Options");
  desc.add_options()
  ("exude", "Exude mesh after refinement")
  ("off_vertices", "Use polyhedron vertices as initialization step")
  ("no_label_rebind", "Don't rebind cell labels in medit output")
  ("show_patches", "Show surface patches in medit output");
  
  
  po::options_description cmdline_options("Usage",1);
  cmdline_options.add(generic).add(mesh).add(desc);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if (vm.count("help") || argc < 2)
  {    
    std::cout << cmdline_options << std::endl;
//    std::cout << "Generic options:\n"
//    << "  --help                Produce help message\n"
//    << "  --file arg            Mesh polyhedron contained in that file\n\n"
//    << "Mesh generation parameters:\n"
//    << "  --facet_angle arg     Set facet angle bound\n"
//    << "  --facet_size arg      Set facet size bound\n"
//    << "  --facet_error arg     Set facet approximation error bound\n"
//    << "  --tet_shape arg       Set tet radius-edge bound\n"
//    << "  --tet_size arg        Set tet size bound\n\n"
//    << "Options:\n"
//    << "  --exude               Exude mesh after refinement\n"
//    << "  --off_vertices        Use polyhedron vertices as initialization step\n
//    << "  --no_label_rebind     Don't rebind cell labels in medit output\n"
//    << "  --show_patches        Show surface patches in medit output\n";
    return 1;
  }
  
  std::cout << "=========== Params ==========="<< std::endl;
  
  double facet_angle = set_arg<double>("facet_angle","Facet angle",vm);
  double facet_size = set_arg<double>("facet_size","Facet size",vm);
  double facet_error = set_arg<double>("facet_error","Facet approximation error",vm);
  
  double tet_shape = set_arg<double>("tet_shape","Tet shape (radius-edge)",vm);
  double tet_size = set_arg<double>("tet_size","Tet size",vm);
  
  std::cout << std::endl;
  std::string polyhedron_filename = set_arg<std::string>("file", "Filename", vm);
  
  std::cout << "=============================="<< std::endl;
  std::cout << std::endl;
  
  if ( polyhedron_filename.empty() )
  {
    std::cout << "No file selected. Exit.\n";
    return 0;
  }

  // Loads polyhedron
  Polyhedron polyhedron;
  std::ifstream input(polyhedron_filename.c_str());
  input >> polyhedron;
  
  using std::atof;
  
  // Domain
  Mesh_domain domain(polyhedron);
  
  // Mesh criteria
  Facet_criteria facet_criteria(facet_angle,
                                facet_size,
                                facet_error); // angle, size, approximation
  Cell_criteria cell_criteria(tet_shape,
                              tet_size); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);
  
  // Meshing
  C3t3 c3t3;
  if ( !vm.count("off_vertices") )
  {
    c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, vm.count("exude"));
  }
  else
  {
    c3t3.insert_surface_points(polyhedron.points_begin(),
                               polyhedron.points_end(),
                               domain.make_surface_index());    
    
    CGAL::refine_mesh_3<C3t3>(c3t3, domain, criteria, vm.count("exude"));
  }
  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file,!vm.count("no_label_rebind"), vm.count("show_patches"));

  return 0;
}
