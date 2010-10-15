#include <CGAL/basic.h>
#include <CGAL/Kinetic/Inexact_simulation_traits.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_vertex_base_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_recent_edges_visitor_2.h>
#include <CGAL/Kinetic/Enclosing_box_2.h>
#include <CGAL/Kinetic/IO/Qt_moving_points_2.h>
#include <CGAL/Kinetic/IO/Qt_triangulation_2.h>
#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Random.h>

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

template <class Traits>
int run(int argc, char *argv[], int n, int d, int seed, std::string file) {
  // typedefs to set everything up
  typedef CGAL::Triangulation_data_structure_2<
  CGAL::Kinetic::Delaunay_triangulation_vertex_base_2<typename Traits::Instantaneous_kernel>,
    CGAL::Kinetic::Delaunay_triangulation_face_base_2<Traits> > TDS;
  typedef CGAL::Delaunay_triangulation_2<typename Traits::Instantaneous_kernel, TDS > Del;
  typedef CGAL::Kinetic::Delaunay_triangulation_recent_edges_visitor_2<Del> Visitor;
  typedef CGAL::Kinetic::Delaunay_triangulation_2<Traits, Visitor, Del> KDel;
  typedef CGAL::Kinetic::Qt_widget_2<typename Traits::Simulator> Qt_gui;
  typedef CGAL::Kinetic::Qt_moving_points_2<Traits, Qt_gui> Qt_mps;
  typedef CGAL::Kinetic::Qt_triangulation_2<KDel, typename Traits::Instantaneous_kernel, Qt_gui> Qt_triangulation;
  typedef CGAL::Kinetic::Enclosing_box_2<Traits> Box;

  CGAL_SET_LOG_LEVEL(CGAL::Log::LOTS);

  Traits tr(0, 10000000);
  typename Box::Handle box= new Box(tr);
  typename KDel::Handle kdel= new KDel(tr);

  typename Qt_gui::Handle qtsim= new Qt_gui(argc, argv, tr.simulator_handle());

  typename Qt_mps::Handle qtmps= new Qt_mps(qtsim, tr);
  typename Qt_triangulation::Handle qtdel= new Qt_triangulation(kdel, tr.instantaneous_kernel_object(), qtsim);


  if (file.empty()) {
    // Generate some random points
    typename CGAL::Random rand= CGAL::Random(seed);
    typename Traits::Active_points_2_table::Key lk;
    std::vector<typename Traits::Kinetic_kernel::Point_2> pts;

    for (int i=0; i< n; ++i) {
      std::vector<double> coefsx, coefsy;
      for (int j=0; j< d; ++j) {
	coefsx.push_back((rand.get_double()*10-5.0)/(j+1));
	coefsy.push_back((rand.get_double()*10-5.0)/(j+1));
	std::cout << coefsx.back() << " " << coefsy.back() << std::endl;
      }
      typename Traits::Kinetic_kernel::Point_2 mp(typename Traits::Kinetic_kernel::Motion_function(coefsx.begin(),
												   coefsx.end()),
						  typename Traits::Kinetic_kernel::Motion_function(coefsy.begin(),
												   coefsy.end()));
      std::cout << "Adding point " << mp << std::endl;
      pts.push_back(mp);
      //std::cout << mp << std::endl;
    }
    for (unsigned int i=0; i< pts.size(); ++i) {
      lk=tr.active_points_2_table_handle()->insert(pts[i]);
    }
    tr.active_points_2_table_handle()->erase(lk);
  } else {
    // read from a file
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
      typename Traits::Kinetic_kernel::Point_2 p;
      il >> p;
      tr.active_points_2_table_handle()->insert(p);
      ++nread;
    }
    std::cout << nread << " points read.\n";
  }

  std::cout << "Green edges just flipped, grey edges will not flip until"
	    << " their certificate changes and black edges will flip." << std::endl;

  return qtsim->begin_event_loop();
}


int main(int argc, char *argv[])
{
  int n=10;
  int d=2;
  int seed=std::time(NULL);
  std::string file;
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  bool print_help=false;
  bool exact=false;
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("exact", boost::program_options::bool_switch(&exact), "Run an exact simulation")
    ("num-points,n", boost::program_options::value<int>(&n), "Number of points to use.")
    ("random-seed,s", boost::program_options::value<int>(&seed), "The value to use for the random seed.")
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
#else
  std::cerr << "Warning, this demo is not very functional without "
            << "boost program options. You probably need to modify "
            << "the code directly.\n";
#endif

  if (true) {
    return run<CGAL::Kinetic::Exact_simulation_traits>(argc, argv, n,d,seed, file);
  } else {
    //return run<CGAL::Kinetic::Inexact_simulation_traits_2>(argc, argv, n,d,seed, file);
  }
}
