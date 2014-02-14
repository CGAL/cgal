#define CGAL_MESH_3_EXAMPLE_POLYHEDRAL_DOMAIN

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//std::size_t TS;

//template <typename C3T3>
//void incident(const C3T3& c3t3)
//{
//  typedef typename C3T3::Triangulation Tr;
//
//  const Tr& tr_ = c3t3.triangulation();
//  for(typename Tr::Finite_vertices_iterator it = tr_.finite_vertices_begin();
//    it!= tr_.finite_vertices_end();
//    ++it)
//  {
//    typename C3T3::Vertex_handle v = it;
//    incident(c3t3,v);
//  }
//}
//
//template <typename C3T3>
//void incident(const C3T3& c3t3,
//  typename C3T3::Vertex_handle vh)
//{
//  typedef typename C3T3::Triangulation Tr;
//
//  const Tr& tr_ = c3t3.triangulation();
//  typedef std::vector<typename Tr::Facet> Facet_vector;
//  Facet_vector facets;
//
//  tr_.finite_incident_facets(vh, std::back_inserter(facets));
//  std::cout << "vertex " << vh->ts << std::endl;
//  for(typename Facet_vector::iterator fit2 = facets.begin() ;
//    fit2 != facets.end() ;
//    ++fit2 )
//  {
//    typename Tr::Cell_handle c = fit2->first;
//    int ii = fit2->second;
//    std::cout << "  f  " << c->ts << " " << ii;
//
//    typename Tr::Facet mf = tr_.mirror_facet(*fit2);
//    c = mf.first;
//    ii = mf.second;
//    std::cout << "  n  " << c->ts << " " << ii << std::endl;
//  }
//}

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#ifdef CGAL_MESH_3_EXAMPLE_POLYHEDRAL_DOMAIN
  #include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#else
  #include <CGAL/Labeled_image_mesh_domain_3.h>
#endif
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
#ifdef CGAL_MESH_3_EXAMPLE_POLYHEDRAL_DOMAIN
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
#else
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#ifdef CGAL_MESH_3_EXAMPLE_POLYHEDRAL_DOMAIN
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;
#else
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
#endif

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;



int main(int argc, char* argv[])
{
  std::cout.precision(17);

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

#ifdef CGAL_MESH_3_EXAMPLE_POLYHEDRAL_DOMAIN
  // Domain
  Mesh_domain domain("data/fandisk.off");
  
  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.02,
                         facet_angle = 30, facet_size = 0.02, facet_distance = 0.002,
                         cell_radius_edge_ratio = 3, cell_size = 0.02);
  //Mesh_criteria criteria(edge_size = 0.025,
  //                       facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
  //                       cell_radius_edge_ratio = 3, cell_size = 0.05);
#else
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
#endif
  
  CGAL::Timer time;
  double total_time = 0.;

  for(std::size_t i = 0; i < nb_runs; ++i)
  {
//    TS = 0;

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
//    std::cerr << "TS = " << TS << std::endl;
  }

  std::cout << "Total time :         " << time.time() << std::endl;
  std::cout << "Time per iteration : " << (time.time() / nb_runs) << std::endl;

  return 0;
}
