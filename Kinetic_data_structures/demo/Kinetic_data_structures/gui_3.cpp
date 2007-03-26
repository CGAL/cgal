#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Regular_triangulation_exact_simulation_traits.h>
#include <boost/program_options.hpp>
#include <CGAL/Kinetic/Enclosing_box_3.h>
#include <CGAL/Random.h>
#include <vector>

#ifdef CGAL_USE_COIN
#include "include/SoQt_widget_3.h"
#include "include/SoQt_moving_points_3.h"
#include "include/SoQt_moving_weighted_points_3.h"
#endif

#include <CGAL/Kinetic/Insert_event.h>

/*!
  \file coin_check.cc A simple example using the coin GUI.
*/

int main(int argc, char *argv[])
{
#ifdef CGAL_USE_COIN
  int n=10;
  int d=2;
  bool print_help=false;
  std::string file;
  bool verbose=false;
  bool weighted=false;
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("verbose,v", boost::program_options::bool_switch(&verbose), "produce lots of output")
    ("weighted,w", boost::program_options::bool_switch(&weighted), "Use weighted points instead of unweighted")
    ("num-points,n", boost::program_options::value<int>(&n), "Number of points to use.")
    ("degree,d", boost::program_options::value<int>(&d), "The degree of the motions to use.")
    ("file,f", boost::program_options::value<std::string>(&file), "Read points from a file.");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(desc).run(), vm);
  boost::program_options::notify(vm);

  if (print_help) {
    std::cout << desc << "\n";
    return EXIT_FAILURE;
  }

  if (weighted) {
    typedef CGAL::Kinetic::Regular_triangulation_exact_simulation_traits Traits;
    typedef CGAL::Kinetic::SoQt_widget_3<Traits::Simulator> Qt_gui;
    typedef CGAL::Kinetic::SoQt_moving_weighted_points_3<Traits, Qt_gui> Qt_mpt;

    Traits tr(0,1000000);
    Qt_gui::Handle qtsim= new Qt_gui(argc, argv, tr.simulator_handle());
    Qt_mpt::Handle qtmpt= new Qt_mpt(tr, qtsim);

    Traits::Simulator::Handle sim= tr.simulator_handle();
    Traits::Active_points_3_table::Handle mpt= tr.active_points_3_table_handle();

    typedef Traits::Kinetic_kernel::Motion_function Fn;

    if (file.empty()) {
      CGAL::Random rand;
      for (int i=0; i< n; ++i) {
	std::vector<double> coefsx, coefsy, coefsz, coefsw;
	for (int j=0; j< d; ++j) {
	  coefsx.push_back((rand.get_double()*10-5)/(j+1));
	  coefsy.push_back((rand.get_double()*10-5)/(j+1));
	  coefsz.push_back((rand.get_double()*10-5)/(j+1));
	  coefsw.push_back((rand.get_double())/(j+1));
	}
	Traits::Kinetic_kernel::Weighted_point_3 mp(Traits::Kinetic_kernel::Point_3(Fn(coefsx.begin(), coefsx.end()),
										    Fn(coefsy.begin(), coefsy.end()),
										    Fn(coefsz.begin(), coefsz.end())),
						    Fn(coefsw.begin(), coefsw.end()));
	tr.active_points_3_table_handle()->insert(mp);
      }
    }
    else {
      std::ifstream in(file.c_str());
      if (!in) {
	std::cerr << "Error opening input file: " << file << std::endl;
	return EXIT_FAILURE;
      }
      char buf[1000];
      int nread=0;
      while (true ) {
	in.getline(buf, 1000);
	if (!in) break;
	std::istringstream il(buf);
	Traits::Kinetic_kernel::Weighted_point_3 p;
	il >> p;
	tr.active_points_3_table_handle()->insert(p);
	++nread;
      }
      std::cout << nread << " points read.\n";
    }

    return qtsim->begin_event_loop();
  }
  else {
    typedef CGAL::Kinetic::Exact_simulation_traits Traits;
    typedef CGAL::Kinetic::SoQt_widget_3<Traits::Simulator> Qt_gui;
    typedef CGAL::Kinetic::SoQt_moving_points_3<Traits, Qt_gui> Qt_mpt;

    Traits tr(0,100000);
    Qt_gui::Handle qtsim= new Qt_gui(argc, argv, tr.simulator_handle());
    Qt_mpt::Handle qtmpt= new Qt_mpt(tr, qtsim);

    typedef Traits::Kinetic_kernel::Motion_function Fn;
    Traits::Simulator::Handle sim= tr.simulator_handle();
    Traits::Active_points_3_table::Handle mpt= tr.active_points_3_table_handle();
    CGAL::Kinetic::Enclosing_box_3<Traits> eb(tr,-10,10,-10,10,-10,10);

    if (file.empty()) {
      CGAL::Random rand;
      for (int i=0; i< n; ++i) {
	std::vector<double> coefsx, coefsy, coefsz;
	for (int j=0; j< d; ++j) {
	  coefsx.push_back((rand.get_double()*10-5)/(j+1));
	  coefsy.push_back((rand.get_double()*10-5)/(j+1));
	  coefsz.push_back((rand.get_double()*10-5)/(j+1));
	}
	Traits::Kinetic_kernel::Point_3 mp(Fn(coefsx.begin(), coefsx.end()),
					   Fn(coefsy.begin(), coefsy.end()),
					   Fn(coefsz.begin(), coefsz.end()));
	tr.active_points_3_table_handle()->insert(mp);
      }
    }
    else {
      std::ifstream in(file.c_str());
      if (!in) {
	std::cerr << "Error opening input file: " << file << std::endl;
	return EXIT_FAILURE;
      }
      char buf[1000];
      int nread=0;
      while (true ) {
	in.getline(buf, 1000);
	if (!in) break;
	std::istringstream il(buf);

	Traits::Kinetic_kernel::Point_3 p;
	il >> p;
	tr.active_points_3_table_handle()->insert(p);
	++nread;
      }
      std::cout << nread << " points read.\n";
    }

    return qtsim->begin_event_loop();
  }
#else
    std::cout << "An install of Inventor and SoQt are required for this demo.  "
    "Please make sure they are installed and then compile "
    "using the makefile 'makefile.soqt'.\n"
    "They can be found at http://www.coin3d.org or as an rpm from "
    "your linux distribution (they are part of Fedora extras, for example).\n";
    return EXIT_FAILURE;
#endif

}
