
// Poor man's profile counters to see which is the failure that leads to exact computation
// in do_intersect(Bbox_3, Segment_3)
int EXIT1, EXIT2, EXIT3, BASE1, BASE2, BASE3, BASE4, BASE5, BASE6, BASE7, BASE8,  BASE9, BASE10, CALLS;

#define ADD_BBOX_POINTS
#define DOUBLE_FILTER
#define DUMP_FAILURES

bool add_bbox_points, double_filter;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

void
bbox_3_to_off(std::ostream& os, const CGAL::Bbox_3& bb)
{
  CGAL::Simple_cartesian<double>::Iso_cuboid_3 ic(bb);
  os << "OFF\n8 6 0\n"; 
  for(int i=0; i <8; i++){
    os << ic[i] << std::endl;
  }
  os << "4 0 1 2 3\n";
  os << "4 0 1 6 5\n";
  os << "4 1 2 7 6\n";
  os << "4 3 2 7 4\n";
  os << "4 0 3 4 5\n";
  os << "4 4 5 6 7" << std::endl;
}


#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Timer.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

#include <string>
#include <boost/lexical_cast.hpp>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  std::cerr << "Usage: mesh_3 [Add_bbox_points (0 or 1)] "
    << "[Double filtering (0 or 1)] [input_filename] [facet_size]" 
    << std::endl << std::endl;
  
  add_bbox_points = double_filter = false;

  if(argc>1){
    add_bbox_points = boost::lexical_cast<int>(argv[1]);
    if(add_bbox_points){
      std::cerr << "Add bbox points" << std::endl;
    }
  }

  if(argc>2){
    double_filter = boost::lexical_cast<int>(argv[2]);
    if(double_filter){
      std::cerr << "Do double filtering" << std::endl;
    }
  }

  char *input_filename = NULL;
  if(argc > 3)
    input_filename = argv[3];
  else
    input_filename = "data/rocker-arm.off";
  
  double facetsize = 0.3034 *0.0045;
  if(argc > 4)
  {
    facetsize = boost::lexical_cast<double>(argv[4]);
  }

  std::cerr.precision(20);
  std::cout.precision(20);
  CGAL::default_random = CGAL::Random(0);

  EXIT1 = EXIT2 = EXIT3 = CALLS = BASE1 = BASE2 = BASE3 = BASE4 = BASE5 = BASE6 = BASE7 = BASE8 = BASE9 = BASE10 = 0;
  // Create input polyhedron
  Polyhedron polyhedron;
  std::cerr << "Loading " << input_filename << "." << std::endl;
  std::ifstream input(input_filename);
  input >> polyhedron;
  if (input.fail())
  {
    std::cerr << "Unable to open file" << std::endl;
    return EXIT_FAILURE;
  }
   
  // Create domain
  Mesh_domain domain(polyhedron);
 
  CGAL::Bbox_3 bb = bbox_3(polyhedron.points_begin(), polyhedron.points_end());
  std::cerr << "bbox(polyhedron) = " << bb << std::endl;

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_size=facetsize);
 
  std::cerr << "Criteria:" << std::endl
    << "  * Facet size = " << facetsize << std::endl;
  std::cerr << "Meshing... ";
  CGAL::Timer t;
  t.start();
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
  t.stop();
  std::cerr << t.time() << std::endl;
  std::cerr <<  BASE1 << "  " <<  BASE2 << "  " <<  BASE3 << "  " <<  BASE4 << "  " <<  BASE5 << "  " <<  BASE6 << "  " <<  BASE7 << "  " <<  BASE8 << "  " <<  BASE9 << "  " << BASE10 << std::endl;
  std::cerr << " "<< CALLS << std::endl;
  std::cerr << " "<< EXIT1 << " "<< EXIT2 << " " << EXIT3 << std::endl;

#ifdef DOUBLE_FILTER
  if (double_filter)
    std::cerr << "Used double filter" << std::endl;
#endif 

#ifdef ADD_BBOX_POINTS
  if (add_bbox_points)
    std::cerr << "Added bbox points" << std::endl;
#endif

  // Output
  // std::ofstream medit_file("fandisk_CGAL.mesh");
  // c3t3.output_to_medit(medit_file);
  // medit_file.close();

  // Set tetrahedron size (keep cell_radius_edge), ignore facets
  //Mesh_criteria new_criteria(cell_radius_edge=3, cell_size=0.03);

  // Mesh refinement
  //CGAL::refine_mesh_3(c3t3, domain, new_criteria);

  // Output
  //medit_file.open("out_2.mesh");
  //c3t3.output_to_medit(medit_file);

  return 0;
}
