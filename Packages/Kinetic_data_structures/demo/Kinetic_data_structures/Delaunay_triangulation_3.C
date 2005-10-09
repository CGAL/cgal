#include <CGAL/KDS/Insert_event.h>
#include <CGAL/KDS/Exact_simulation_traits_3.h>
#include <algorithm>
#include <CGAL/KDS/Delaunay_triangulation_3.h>
#include <CGAL/KDS/IO/Qt_widget_3.h>
#include <CGAL/KDS/IO/Qt_moving_points_3.h>
#include <CGAL/KDS/IO/Qt_triangulation_3.h>
#include <boost/program_options.hpp>
#include <CGAL/KDS/Enclosing_box_3.h>



int main(int argc, char *argv[]){
  int n=10;
  int d=2;
  bool print_help=false;
  std::string file;
  bool verbose=false;
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("verbose,v", boost::program_options::bool_switch(&verbose), "produce lots of output")
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


  typedef CGAL::KDS::Exact_simulation_traits_3 Traits;
  typedef CGAL::KDS::Delaunay_triangulation_3<Traits> KDel;
  typedef Traits::Simulator::Time Time;
  typedef CGAL::KDS::Insert_event<Traits::Moving_point_table> IE;
  typedef Traits::Moving_point_table::Data MP;
  typedef Traits::Kinetic_kernel::Motion_function MF;
  typedef Traits::Kinetic_kernel::Point_3 MP;
  typedef CGAL::KDS::Qt_widget_3<Traits::Simulator> Qt_gui;
  typedef CGAL::KDS::Qt_moving_points_3<Traits, Qt_gui> Qt_mpt;
  typedef CGAL::KDS::Qt_triangulation_3<KDel, Qt_gui, Qt_mpt> Qt_del;

  Traits tr;
  Qt_gui::Pointer qtsim= new Qt_gui(argc, argv, tr.simulator_pointer());
  Qt_mpt::Pointer qtmpt= new Qt_mpt(tr, qtsim);
  KDel::Pointer kdel= new KDel(tr);
  Qt_del::Pointer cd= new Qt_del(kdel, qtsim, qtmpt);

  
 


  
  
  typedef Traits::Simulator::Time Time;
  typedef CGAL::KDS::Insert_event<Traits::Moving_point_table> MOI;
  
  Traits::Function_kernel::Construct_function cf= tr.function_kernel_object().construct_function_object();
  //sim->end_time();
  Traits::Simulator::Pointer sim= tr.simulator_pointer();
  Traits::Moving_point_table::Pointer mpt= tr.moving_point_table_pointer();
  CGAL::KDS::Enclosing_box_3<Traits> eb(tr,-10,10,-10,10,-10,10);

  if (verbose) {
    CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_LOTS);
  } else {
    CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_NONE);
  }

  if (file.empty()) {
    CGAL::Random rand;
    for (int i=0; i< n; ++i){
      std::vector<double> coefsx, coefsy, coefsz;
      for (int j=0; j< d; ++j){
	coefsx.push_back((rand.get_double()*10-5)/(j+1));
	coefsy.push_back((rand.get_double()*10-5)/(j+1));
	coefsz.push_back((rand.get_double()*10-5)/(j+1));
      }
      Traits::Kinetic_kernel::Point_3 mp(Traits::Function_kernel::Function(coefsx.begin(), coefsx.end()),
					 Traits::Function_kernel::Function(coefsy.begin(), coefsy.end()),
					 Traits::Function_kernel::Function(coefsz.begin(), coefsz.end()));
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

  kdel->set_has_certificates(true);
  return qtsim->begin_event_loop();
};
