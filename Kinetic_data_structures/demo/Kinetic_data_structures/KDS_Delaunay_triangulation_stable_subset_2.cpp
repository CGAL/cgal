#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Random.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Inexact_simulation_traits.h>
#include "Qt_Delaunay_stable_subset_2.h"
#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Kinetic/IO/Qt_moving_points_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_recent_edges_visitor_2.h>
#include <CGAL/Kinetic/Enclosing_box_2.h>

#include <CGAL/Random.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif


int main(int argc, char *argv[])
{
  double threshold=.9;
  int n=10;
  int d=2;
  bool print_help=false;
  std::string file;
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("num-points,n", boost::program_options::value<int>(&n), "Number of points to use.")
    ("degree,d", boost::program_options::value<int>(&d), "The degree of the motions to use.")
    ("file,f", boost::program_options::value<std::string>(&file), "Read points from a file.")
    ("threshold,t", boost::program_options::value<double>(&threshold), "The threshold for displaying the edges.");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(desc).run(), vm);
  boost::program_options::notify(vm);

  if (print_help) {
    std::cout << desc << "\n";
    return EXIT_FAILURE;
  }
#endif

  typedef CGAL::Kinetic::Inexact_simulation_traits Traits;
  typedef CGAL::Triangulation_data_structure_2<
  CGAL::Kinetic::Delaunay_triangulation_vertex_base_2<Traits::Instantaneous_kernel>,
    CGAL::Kinetic::Delaunay_triangulation_face_base_2<Traits> > TDS;
  typedef CGAL::Delaunay_triangulation_2<Traits::Instantaneous_kernel, TDS > Del;
  typedef CGAL::Kinetic::Enclosing_box_2<Traits> EB;
  typedef CGAL::Kinetic::Delaunay_triangulation_recent_edges_visitor_2<TDS> Visitor;
  typedef CGAL::Kinetic::Delaunay_triangulation_2<Traits, Visitor, Del> KDel;
  typedef CGAL::Kinetic::Qt_widget_2<Traits::Simulator> Qt_gui;
  typedef CGAL::Kinetic::Qt_moving_points_2<Traits, Qt_gui> Qt_mps;
  typedef CGAL::Kinetic::Qt_Delaunay_stable_subset_2<KDel, Qt_gui, Qt_mps> Qt_triangulation;

  Traits tr(0,100000.0);

  KDel::Handle kdel= new KDel(tr);
  EB::Handle eb= new EB(tr,-10,10,-10,10);
  Qt_gui::Handle qtsim= new Qt_gui(argc, argv, tr.simulator_handle(), -10,10,-10,10);
  Qt_mps::Handle qtmps= new Qt_mps(qtsim, tr);

  Qt_triangulation::Handle qtdel= new Qt_triangulation(qtsim, qtmps, kdel, threshold);

  if (file.empty()) {
    CGAL::Random rand;
    for (int i=0; i< n; ++i) {
      std::vector<double> coefsx, coefsy;
      for (int j=0; j< d; ++j) {
	coefsx.push_back((rand.get_double()*10-5)/(j+1));
	coefsy.push_back((rand.get_double()*10-5)/(j+1));
      }
      Traits::Kinetic_kernel::Point_2 mp(Traits::Kinetic_kernel::Motion_function(coefsx.begin(), coefsx.end()),
					 Traits::Kinetic_kernel::Motion_function(coefsy.begin(), coefsy.end()));
      tr.active_points_2_table_handle()->insert(mp);
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
      Traits::Kinetic_kernel::Point_2 p;
      il >> p;
      tr.active_points_2_table_handle()->insert(p);
      ++nread;
    }
    std::cout << nread << " points read.\n";
  }

  kdel->set_has_certificates(true);
  return qtsim->begin_event_loop();
}
