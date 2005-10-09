#include <CGAL/basic.h>

#include <CGAL/KDS/Exact_simulation_traits_2.h>
#include <CGAL/KDS/Delaunay_triangulation_2.h>
#include <CGAL/KDS/Delaunay_triangulation_recent_edges_visitor_2.h>
#include <CGAL/KDS/Enclosing_box_2.h>
#include <CGAL/KDS/IO/Qt_moving_points_2.h>
#include <CGAL/KDS/IO/Qt_triangulation_2.h>
#include <CGAL/KDS/IO/Qt_widget_2.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/Random.h>
#include <CGAL/Random.h>
#include <algorithm>
#include <boost/program_options.hpp>





int main(int argc, char *argv[]){
  int n=10;
  int d=2;
  bool print_help=false;
  std::string file;
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
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

  typedef CGAL::KDS::Exact_simulation_traits_2 Traits;
  typedef CGAL::Triangulation_data_structure_2<
    CGAL::Triangulation_vertex_base_2<Traits::Instantaneous_kernel>, 
    CGAL::KDS::Delaunay_triangulation_face_base_2<Traits> > TDS;
  typedef CGAL::Delaunay_triangulation_2<Traits::Instantaneous_kernel, TDS > Del;
  typedef CGAL::KDS::Delaunay_triangulation_recent_edges_visitor_2<TDS> Visitor;
  typedef CGAL::KDS::Delaunay_triangulation_2<Traits, Visitor, Del> KDel;
  typedef CGAL::KDS::Qt_widget_2<Traits::Simulator> Qt_gui;
  typedef CGAL::KDS::Qt_moving_points_2<Traits, Qt_gui> Qt_mps;
  typedef CGAL::KDS::Qt_triangulation_2<KDel, Qt_gui, Qt_mps> Qt_triangulation;
  typedef CGAL::KDS::Enclosing_box_2<Traits> Box;

  //CGAL_KDS_SET_LOG_LEVEL(CGAL::KDS::LOG_LOTS);
  


  Traits tr;
  Box::Pointer box= new Box(tr);
  KDel::Pointer kdel= new KDel(tr);
  
  Qt_gui::Pointer qtsim= new Qt_gui(argc, argv, tr.simulator_pointer());
  Qt_mps::Pointer qtmps= new Qt_mps(qtsim, tr);

  Qt_triangulation::Pointer qtdel= new Qt_triangulation(kdel, qtsim, qtmps);

  if (file.empty()) {
    CGAL::Random rand;
    for (int i=0; i< n; ++i){
      std::vector<double> coefsx, coefsy;
      for (int j=0; j< d; ++j){
	coefsx.push_back((rand.get_double()*10-5)/(j+1));
	coefsy.push_back((rand.get_double()*10-5)/(j+1));
      }
      Traits::Kinetic_kernel::Point_2 mp(Traits::Function_kernel::Function(coefsx.begin(), coefsx.end()),
					 Traits::Function_kernel::Function(coefsy.begin(), coefsy.end()));
      tr.moving_point_table_pointer()->insert(mp);
      //std::cout << mp << std::endl;
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
      Traits::Kinetic_kernel::Point_2 p;
      il >> p;
      tr.moving_point_table_pointer()->insert(p);
      ++nread;
    }
    std::cout << nread << " points read.\n";
  }

  return qtsim->begin_event_loop();
};
