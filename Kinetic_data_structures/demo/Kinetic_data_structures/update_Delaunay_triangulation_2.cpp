#define CGAL_CHECK_EXPENSIVE
//#define CGAL_CHECK_EXACTNESS

#include <CGAL/basic.h>

#ifdef CGAL_USE_QT

#include <CGAL/Kinetic/IO/Qt_moving_points_2.h>
#include <CGAL/Kinetic/IO/Qt_triangulation_2.h>
#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Random.h>

#include <CGAL/Updatable_Delaunay_triangulation_2.h>


#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

#define UD_DEBUG(x) 


int main(int argc, char *argv[]) {
  int n=10;
  int d=2;
  int seed=std::time(NULL);
  bool print_help=false;
  bool exact=false;
  std::string ifile, ffile;
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("exact", boost::program_options::bool_switch(&exact), "Run an exact simulation")
    ("num-points,n", boost::program_options::value<int>(&n), "Number of points to use.")
    ("random-seed,s", boost::program_options::value<int>(&seed), "The value to use for the random seed.")
    ("degree,d", boost::program_options::value<int>(&d), "The degree of the motions to use.")
    ("ifile,i", boost::program_options::value<std::string>(&ifile), "The inital coordinates.")
    ("ffile,f", boost::program_options::value<std::string>(&ffile), "The final coordinates.")
    ;

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(desc).run(), vm);
  boost::program_options::notify(vm);

  if (print_help) {
    std::cout << desc << "\n";
    return EXIT_FAILURE;
  }
#endif
  
  //if (exact) {
  typedef CGAL::Updatable_Delaunay_triangulation_2<CGAL::Indirect_point_2_kernel<CGAL::Exact_predicates_inexact_constructions_kernel> > UD;
  return UD::run(argc, argv, n,d,seed, ifile, ffile);
    /*} else {
    Update_Delaunay_triangulation_2<CGAL::Kinetic::Inexact_simulation_traits_2> ud;
    return ud.run(argc, argv, n,d,seed, ifile, ffile);
    }*/
};
#else
int main(int, char *[]){
  return 0;
}
#endif
