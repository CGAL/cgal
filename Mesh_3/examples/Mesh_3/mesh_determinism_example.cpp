#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/lloyd_optimize_mesh_3.h>
#include <CGAL/odt_optimize_mesh_3.h>
#include <CGAL/perturb_mesh_3.h>
#include <CGAL/exude_mesh_3.h>
#include <CGAL/Image_3.h>

#include <sstream>
#include <cstring>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  // Collect options
  std::size_t nb_runs = 1;
  bool do_lloyd = false;
  bool do_odt = false;
  bool do_perturb = false;
  bool do_exude = false;
  for(int i = 1; i < argc; ++i) 
  {
    std::string arg = argv[i];
    if(arg == "-n")             nb_runs = atoi(argv[i+1]);
    else if(arg == "-lloyd")    do_lloyd = true;
    else if(arg == "-odt")      do_odt = true;
    else if(arg == "-perturb")  do_perturb = true;
    else if(arg == "-exude")    do_exude = true;
  }

  // Domain
  CGAL::Image_3 image;
  image.read("data/liver.inr.gz");
  Mesh_domain domain(image);

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, 
                         facet_size=5, 
                         facet_distance=1.5,
                         cell_radius_edge_ratio=2, 
                         cell_size=7);

  for(std::size_t i = 0; i < nb_runs; ++i)
  {
    std::ostringstream oss;
    oss << i << "_";
    char* num_str = (char*)(oss.str().data());

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        no_perturb(),
                                        no_exude());

    std::ofstream medit_file(strcat(num_str,"out0-refinement.mesh"));
    c3t3.output_to_medit(medit_file);
  
    //LLOYD
    CGAL::lloyd_optimize_mesh_3(c3t3, domain, max_iteration_number = 10);
    std::ofstream medit_file1(strcat(num_str,"out1-lloyd.mesh"));
    c3t3.output_to_medit(medit_file1);

    //ODT
    CGAL::odt_optimize_mesh_3(c3t3, domain, max_iteration_number = 10);
    std::ofstream medit_file2(strcat(num_str,"out2-odt.mesh"));
    c3t3.output_to_medit(medit_file2);

    //PERTURB
    CGAL::perturb_mesh_3(c3t3, domain, sliver_bound=10);
    std::ofstream medit_file3(strcat(num_str,"out3-perturb.mesh"));
    c3t3.output_to_medit(medit_file3);

    //EXUDE
    CGAL::exude_mesh_3(c3t3, domain, sliver_bound=12);
    std::ofstream medit_file4(strcat(num_str,"out4-exude.mesh"));
    c3t3.output_to_medit(medit_file4);
  }

  return 0;
}
