//#define MOVE_ALL
//#define HYBRID

#ifndef NDEBUG
#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS
//#define CGAL_KINETIC_DISABLE_AUDITING
#endif

#include <CGAL/basic.h>

#include <CGAL/Updatable_Delaunay_triangulation_2.h>

#include <CGAL/Kinetic/IO/Qt_moving_points_2.h>
#include <CGAL/Kinetic/IO/Qt_triangulation_2.h>
#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Random.h>

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif



int main(int argc, char *argv[]) {
  int n=10;
  int d=2;
  int seed=std::time(NULL);
  bool print_help=false;
  bool exact=false;
  std::string ifile="data/before002", ffile="data/after002";
  std::cout.precision(18); // std::ios::precision(18);
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("disable-0,0", boost::program_options::bool_switch(&CGAL::disable_filter_0_), "Disable filter 0")
    ("disable-1,1", boost::program_options::bool_switch(&CGAL::disable_filter_1_), "Disable filter 1")
    ("disable-2,2", boost::program_options::bool_switch(&CGAL::disable_filter_2_), "Disable filter 2")
    ("disable-3,3", boost::program_options::bool_switch(&CGAL::disable_filter_3_), "Disable filter 3")
    ("exact", boost::program_options::bool_switch(&exact), "Run an exact simulation")
    ("num-points,n", boost::program_options::value<int>(&n), "Number of points to use.")
    ("random-seed,s", boost::program_options::value<int>(&seed), "The value to use for the random seed.")
    //("degree,d", boost::program_options::value<int>(&d), "The degree of the motions to use.")
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
  typedef CGAL::Updatable_Delaunay_triangulation_2<CGAL::Indirect_point_2_kernel<CGAL::Exact_predicates_inexact_constructions_kernel> > UD;
  return UD::run(argc, argv, n,d,seed, ifile, ffile);   
};
