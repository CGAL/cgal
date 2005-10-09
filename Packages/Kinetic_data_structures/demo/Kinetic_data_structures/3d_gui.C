#include <CGAL/KDS/Exact_simulation_traits_3.h>
#include <CGAL/KDS/Regular_triangulation_exact_simulation_traits_3.h>
#include <boost/program_options.hpp>
#include <CGAL/KDS/Enclosing_box_3.h>
#include <CGAL/Random.h>
#include <vector>


#ifdef CGAL_USE_COIN
#include <CGAL/KDS/IO/Qt_widget_3.h>
#include <CGAL/KDS/IO/Qt_moving_points_3.h>
#include <CGAL/KDS/IO/Qt_moving_weighted_points_3.h>
#endif

#include <CGAL/KDS/Insert_event.h>

/*!
  \file coin_check.cc A simple example using the coin GUI.
*/

int main(int argc, char *argv[]){
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
    typedef CGAL::KDS::Regular_triangulation_exact_simulation_traits_3 Traits;
    typedef CGAL::KDS::Qt_widget_3<Traits::Simulator> Qt_gui;
    typedef CGAL::KDS::Qt_moving_weighted_points_3<Traits, Qt_gui> Qt_mpt;
    
    Traits tr;
    Qt_gui::Pointer qtsim= new Qt_gui(argc, argv, tr.simulator_pointer());
    Qt_mpt::Pointer qtmpt= new Qt_mpt(tr, qtsim);

    Traits::Function_kernel::Construct_function cf= tr.function_kernel_object().construct_function_object();
    Traits::Simulator::Pointer sim= tr.simulator_pointer();
    Traits::Moving_point_table::Pointer mpt= tr.moving_point_table_pointer();
    
    if (file.empty()) {
      CGAL::Random rand;
      for (int i=0; i< n; ++i){
	std::vector<double> coefsx, coefsy, coefsz, coefsw;
	for (int j=0; j< d; ++j){
	  coefsx.push_back((rand.get_double()*10-5)/(j+1));
	  coefsy.push_back((rand.get_double()*10-5)/(j+1));
	  coefsz.push_back((rand.get_double()*10-5)/(j+1));
	  coefsw.push_back((rand.get_double())/(j+1));
	}
	Traits::Kinetic_kernel::Weighted_point_3 mp(Traits::Kinetic_kernel::Point_3(cf(coefsx.begin(), coefsx.end()),
										    cf(coefsy.begin(), coefsy.end()),
										    cf(coefsz.begin(), coefsz.end())),
						     cf(coefsw.begin(), coefsw.end()));
	tr.moving_point_table_pointer()->insert(mp);
      }
    } else {
      std::ifstream in(file.c_str());
      if (!in) {
	std::cerr << "Error opening input file: " << file << std::endl;
	return EXIT_FAILURE;
      }
      char buf[1000];
      int nread=0;
      while (true ){
	in.getline(buf, 1000);
	if (!in) break;
	std::istringstream il(buf);
	Traits::Kinetic_kernel::Weighted_point_3 p;
	il >> p;
	tr.moving_point_table_pointer()->insert(p);
	++nread;
      }
      std::cout << nread << " points read.\n";
    }

    return qtsim->begin_event_loop();
  } else {
    typedef CGAL::KDS::Exact_simulation_traits_3 Traits;
    typedef CGAL::KDS::Qt_widget_3<Traits::Simulator> Qt_gui;
    typedef CGAL::KDS::Qt_moving_points_3<Traits, Qt_gui> Qt_mpt;
    
    Traits tr;
    Qt_gui::Pointer qtsim= new Qt_gui(argc, argv, tr.simulator_pointer());
    Qt_mpt::Pointer qtmpt= new Qt_mpt(tr, qtsim);  
  
    Traits::Function_kernel::Construct_function cf= tr.function_kernel_object().construct_function_object();
    
    Traits::Simulator::Pointer sim= tr.simulator_pointer();
    Traits::Moving_point_table::Pointer mpt= tr.moving_point_table_pointer();
    CGAL::KDS::Enclosing_box_3<Traits> eb(tr,-10,10,-10,10,-10,10);

    if (file.empty()) {
      CGAL::Random rand;
      for (int i=0; i< n; ++i){
	std::vector<double> coefsx, coefsy, coefsz;
	for (int j=0; j< d; ++j){
	  coefsx.push_back((rand.get_double()*10-5)/(j+1));
	  coefsy.push_back((rand.get_double()*10-5)/(j+1));
	  coefsz.push_back((rand.get_double()*10-5)/(j+1));
	}
	Traits::Kinetic_kernel::Point_3 mp(cf(coefsx.begin(), coefsx.end()),
					   cf(coefsy.begin(), coefsy.end()),
					   cf(coefsz.begin(), coefsz.end()));
	tr.moving_point_table_pointer()->insert(mp);
      }
    } else {
      std::ifstream in(file.c_str());
      if (!in) {
	std::cerr << "Error opening input file: " << file << std::endl;
	return EXIT_FAILURE;
      }
      char buf[1000];
      int nread=0;
      while (true ){
	in.getline(buf, 1000);
	if (!in) break;
	std::istringstream il(buf);
      
	Traits::Kinetic_kernel::Point_3 p;
	il >> p;
	tr.moving_point_table_pointer()->insert(p);
	++nread;
      }
      std::cout << nread << " points read.\n";
    }

    return qtsim->begin_event_loop();
  }
#else
    bool warning_coin_is_not_available_therefore_no_3d_visualization;
    return 0;
#endif

}
