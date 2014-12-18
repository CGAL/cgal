#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>
#include <CGAL/Random.h>
#include <CGAL/Mesh_3/Profiling_tools.h>

#include "../../test/Tangential_complex/test_utilities.h"


#include <fstream>
#include <math.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif
#include "XML_exporter.h"
#define CGAL_TC_EXPORT_PERFORMANCE_DATA
#define CGAL_TC_SET_PERFORMANCE_DATA(value_name, value) \
        XML_perf_data::set(value_name, value);

class XML_perf_data
{
public:
  typedef Streaming_XML_exporter<std::string> XML_exporter;

  XML_perf_data(const std::string &filename)
    : m_xml(filename, "ContainerPerformance", "Perf",
            construct_subelements_names())
  {}

  virtual ~XML_perf_data()
  {
  }

  static XML_perf_data &get()
  {
    static XML_perf_data singleton(build_filename());
    return singleton;
  }

  template <typename Value_type>
  static void set(const std::string &name, Value_type value)
  {
    get().set_data(name, value);
  }

  static void commit()
  {
    get().commit_current_element();
  }

protected:
  static std::string build_filename()
  {
    std::stringstream sstr;
    sstr << "Performance_log_" << time(0) << ".xml";
    return sstr.str();
  }

  static std::vector<std::string> construct_subelements_names()
  {
    std::vector<std::string> subelements;
    subelements.push_back("Name");
    subelements.push_back("Intrinsic_dim");
    subelements.push_back("Ambient_dim");
    subelements.push_back("Num_points");
    subelements.push_back("Noise");
    subelements.push_back("Init_time");
    subelements.push_back("Comput_time");
    subelements.push_back("Fix_time");

    return subelements;
  }

  void set_data(const std::string &name, const std::string &value)
  {
    m_current_element[name] = value;
  }

  template <typename Value_type>
  void set_data(const std::string &name, Value_type value)
  {
    std::stringstream sstr;
    sstr << value;
    set_data(name, sstr.str());
  }

  void commit_current_element()
  {
    m_xml.add_element(m_current_element);
    m_current_element.clear();
  }

  XML_exporter m_xml;
  XML_exporter::Element_with_map m_current_element;
};

#ifdef _DEBUG
  const int NUM_POINTS = 50;
#else
  const int NUM_POINTS = 30000;
#endif

int main()
{
#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  tbb::task_scheduler_init init(1);
# else
  tbb::task_scheduler_init init(10);
# endif
#endif

  const int INTRINSIC_DIMENSION = 2;
  const int AMBIENT_DIMENSION   = 4;

  typedef CGAL::Epick_d<CGAL::Dimension_tag<AMBIENT_DIMENSION> > Kernel;
  typedef Kernel::FT                                             FT;
  typedef Kernel::Point_d                                        Point;
 
  int i = 0;
  bool stop = false;
  //for ( ; !stop ; ++i)
  {
    Kernel k;
    Wall_clock_timer t;
    CGAL::default_random = CGAL::Random(i);
    std::cerr << "Random seed = " << i << std::endl;
    
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t_gen;
#endif
    //std::vector<Point> points;
    //load_points_from_file<Point>(
    //  "data/points_10_10k.cin", std::back_inserter(points)/*, 600*/);

    std::vector<Point> points =
      //generate_points_on_circle_2<Kernel>(NUM_POINTS, 3.);
      //generate_points_on_moment_curve<Kernel>(NUM_POINTS, AMBIENT_DIMENSION, 0., 1.);
      //generate_points_on_plane<Kernel>(NUM_POINTS);
      //generate_points_on_sphere_3<Kernel>(NUM_POINTS, 3.0);
      //generate_points_on_sphere_d<Kernel>(NUM_POINTS, AMBIENT_DIMENSION, 3.0);
      //generate_points_on_klein_bottle_3D<Kernel>(NUM_POINTS, 4., 3.);
      generate_points_on_klein_bottle_4D<Kernel>(NUM_POINTS, 4., 3.);
      //generate_points_on_klein_bottle_variant_5D<Kernel>(NUM_POINTS, 4., 3.);
    

    CGAL_TC_SET_PERFORMANCE_DATA("Name", "Klein bottle 4D");
    CGAL_TC_SET_PERFORMANCE_DATA("Intrinsic_dim", INTRINSIC_DIMENSION);
    CGAL_TC_SET_PERFORMANCE_DATA("Ambient_dim", AMBIENT_DIMENSION);
    CGAL_TC_SET_PERFORMANCE_DATA("Noise", 0.01);
    
#ifdef CGAL_TC_PROFILING
    std::cerr << "Point set generated in " << t_gen.elapsed()
              << " seconds." << std::endl;
#endif
    
    std::cerr << "Number of points before sparsification: "
      << points.size() << std::endl;
    //points = sparsify_point_set(
    //  k, points, FT(INPUT_SPARSITY)*FT(INPUT_SPARSITY));
    std::cerr << "Number of points after sparsification: "
      << points.size() << std::endl;

    CGAL::Tangential_complex<
      Kernel,
      INTRINSIC_DIMENSION,
      CGAL::Parallel_tag> tc(points.begin(), points.end(), k);
    double init_time = t.elapsed(); t.reset();

    //tc.estimate_intrinsic_dimension();
    tc.compute_tangential_complex();
    double computation_time = t.elapsed(); t.reset();

    std::set<std::set<std::size_t> > incorrect_simplices;
    //stop = !tc.check_if_all_simplices_are_in_the_ambient_delaunay(&incorrect_simplices);

    t.reset();
    tc.fix_inconsistencies();
    double fix_time = t.elapsed(); t.reset();

    double export_time = -1.;
    if (INTRINSIC_DIMENSION <= 3)
    {
      t.reset();
      std::stringstream output_filename;
      output_filename << "data/test_tc_" << INTRINSIC_DIMENSION
        << "_in_R" << AMBIENT_DIMENSION << ".off";
      std::ofstream off_stream(output_filename.str().c_str());
      tc.export_to_off(off_stream, true, &incorrect_simplices, true);
      double export_time = t.elapsed(); t.reset();
    }
    /*else
      tc.number_of_inconsistent_simplices();*/


    std::cerr << std::endl
      << "================================================" << std::endl
      << "Number of vertices: " << tc.number_of_vertices() << std::endl
      << "Computation times (seconds): " << std::endl
      << "  * Tangential complex: " << init_time + computation_time
      << std::endl
      << "    - Init + kd-tree = " << init_time << std::endl
      << "    - TC computation = " << computation_time << std::endl
      << "  * Fix inconsistencies: " << fix_time << std::endl
      << "  * Export to OFF: " << export_time << std::endl
      << "================================================" << std::endl
      << std::endl;
    
    CGAL_TC_SET_PERFORMANCE_DATA("Num_points", tc.number_of_vertices());
    CGAL_TC_SET_PERFORMANCE_DATA("Init_time", init_time);
    CGAL_TC_SET_PERFORMANCE_DATA("Comput_time", computation_time);
    CGAL_TC_SET_PERFORMANCE_DATA("Fix_time", fix_time);
  }

  return 0;
}
