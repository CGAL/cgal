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
#include <CGAL/Timer.h>

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
  char* filename = "run";
  bool do_lloyd = false;
  unsigned int nb_lloyd = 1;
  bool do_odt = false;
  unsigned int nb_odt = 1;
  bool do_perturb = false;
  double perturb_bound = 10.;
  bool do_exude = false;
  double exude_bound = 15.;

  for(int i = 1; i < argc; ++i) 
  {
    std::string arg = argv[i];
    if(arg == "-n")             nb_runs = atoi(argv[i+1]);
    else if(arg == "-name")     filename = argv[i+1];
    else if(arg == "-lloyd")
    {
      do_lloyd = true;
      nb_lloyd = atoi(argv[i+1]);
      ++i;
    }
    else if(arg == "-odt")
    {
      do_odt = true;
      nb_odt = atoi(argv[i+1]);
      ++i;
    }
    else if(arg == "-perturb")
    {
      do_perturb = true;
      perturb_bound = atof(argv[i+1]);
      ++i;
    }
    else if(arg == "-exude")
    {
      do_exude = true;
      exude_bound = atof(argv[i+1]);
      ++i;
    }
  }

  // Domain
  CGAL::Image_3 image;
  image.read("data/liver.inr.gz");
  Mesh_domain domain(image);

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, 
                         facet_size=5,//3, 
                         facet_distance=1.5,
                         cell_radius_edge_ratio=2, 
                         cell_size=7);//3);

  CGAL::Timer time;
  double total_time = 0.;

  for(std::size_t i = 0; i < nb_runs; ++i)
  {
    CGAL::default_random = CGAL::Random(0);

    std::cout << "------- Iteration " << (i+1) << " -------" << std::endl;

    std::ostringstream oss;
    oss << filename << (i+1) << "_";
    std::string num_str(oss.str().data());

    time.start();
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        no_perturb(),
                                        no_exude());
    time.stop();

    std::ofstream medit_file(num_str + std::string("out0-refinement.mesh"));
    c3t3.output_to_medit(medit_file);
  
    //LLOYD
    if(do_lloyd)
    {
      time.start();
      CGAL::lloyd_optimize_mesh_3(c3t3, domain, max_iteration_number = nb_lloyd);
      time.stop();
      std::ofstream medit_file1(num_str + std::string("out1-lloyd.mesh"));
      c3t3.output_to_medit(medit_file1);
    }
    //ODT
    if(do_odt)
    {
      time.start();
      CGAL::odt_optimize_mesh_3(c3t3, domain, max_iteration_number = nb_odt);
      time.stop();
      std::ofstream medit_file2(num_str + std::string("out2-odt.mesh"));
      c3t3.output_to_medit(medit_file2);
    }
    //PERTURB
    if(do_perturb)
    {
      time.start();
      CGAL::perturb_mesh_3(c3t3, domain, sliver_bound=perturb_bound);
      time.stop();
      std::ofstream medit_file3(num_str + std::string("out3-perturb.mesh"));
      c3t3.output_to_medit(medit_file3);
    }
    //EXUDE
    if(do_exude)
    {
      time.start();
      CGAL::exude_mesh_3(c3t3, sliver_bound=exude_bound);
      time.stop();
      std::ofstream medit_file4(num_str + std::string("out4-exude.mesh"));
      c3t3.output_to_medit(medit_file4);
    }
    std::cout << "[Timer at " << time.time() << " sec]" << std::endl;
  }

  std::cout << "Total time :         " << time.time() << std::endl;
  std::cout << "Time per iteration : " << (time.time() / nb_runs) << std::endl;

  return 0;
}
