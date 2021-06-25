#include <CGAL/basic.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Benchmark/config.hpp>
#if defined(CGAL_TRAITS_COUNTING)
#include <CGAL/Arr_counting_traits_2.h>
#endif
#if defined(CGAL_TRAITS_TRACING)
#include <CGAL/Arr_tracing_traits_2.h>
#endif

#include "bench_config.hpp"
#include "number_type.hpp"
#include "kernel_type.hpp"

CGAL_BENCHMARK_BEGIN_NAMESPACE
CGAL_BENCHMARK_END_NAMESPACE

namespace cb = CGAL::benchmark;

enum MaxFilesNumber {
  MAX_FILES_NUMBER = 20
};

// PostScript support:
#if BENCH_TRAITS != LEDA_CONIC_TRAITS && BENCH_TRAITS != CORE_CONIC_TRAITS && \
    BENCH_TRAITS != EXACUS_CONIC_TRAITS && BENCH_TRAITS != CK_CIRCLE_TRAITS && \
    BENCH_TRAITS != CK_CONIC_TRAITS
#define POSTSCRIPT_SUPPORTED 1
#endif
#undef POSTSCRIPT_SUPPORTED

// Landmark point-location supported:
#if BENCH_TRAITS != LEDA_CONIC_TRAITS && BENCH_TRAITS != CORE_CONIC_TRAITS && \
    BENCH_TRAITS != POLYLINE_TRAITS && \
    BENCH_TRAITS != NON_CACHING_POLYLINE_TRAITS
#define LANDMARK_SUPPORTED 1
#endif
#undef LANDMARK_SUPPORTED

// Triangle point-location supported:
#if BENCH_TRAITS != LEDA_CONIC_TRAITS && BENCH_TRAITS != CORE_CONIC_TRAITS && \
    BENCH_TRAITS != POLYLINE_TRAITS && \
    BENCH_TRAITS != NON_CACHING_POLYLINE_TRAITS
#define TRIANGLE_SUPPORTED 1
#endif
#undef TRIANGLE_SUPPORTED

// Stream:
#if BENCH_KERNEL == LEDA_KERNEL || BENCH_KERNEL == MY_KERNEL
#if defined(USE_CGAL_WINDOW)
#include <LEDA/graphics/rat_window.h>
#else
#include <CGAL/IO/Qt_widget_Leda_rat.h>
#endif
#endif

// Traits:
#if BENCH_TRAITS == SEGMENT_TRAITS
#include <CGAL/Arr_segment_traits_2.h>

#elif BENCH_TRAITS == NON_CACHING_SEGMENT_TRAITS
#include <CGAL/Arr_non_caching_segment_traits_2.h>

#elif BENCH_TRAITS == LEDA_SEGMENT_TRAITS
#include <CGAL/Arr_leda_segment_traits_2.h>

#elif BENCH_TRAITS == POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/IO/Qt_widget_Polyline_2.h>

#elif BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#include <CGAL/IO/Qt_widget_Polyline_2.h>

#elif BENCH_TRAITS == LEDA_CONIC_TRAITS
#include <CGAL/Arr_conic_traits_2.h>
#if defined(USE_CGAL_WINDOW)
#include <CGAL/IO/Conic_arc_2_Window_stream.h>
#else
#include <CGAL/IO/Qt_widget_Conic_arc_2.h>
#endif

#elif BENCH_TRAITS == CORE_CONIC_TRAITS
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#if defined(USE_CGAL_WINDOW)

//#include <CGAL/IO/Conic_arc_2_Window_stream_core.h>
//#include <CGAL/Arrangement_2/Conic_arc_2_core.h>

#include <CGAL/IO/Conic_arc_2_Window_stream.h>

#else
#include <CGAL/IO/Qt_widget_Conic_arc_2.h>
#endif

// Exacus Conics:
#elif BENCH_TRAITS == EXACUS_CONIC_TRAITS
#include <CnX/Conic_sweep_traits_2.h>
#include <CnX/Conic_segment_v2_2.h>
#include <SoX/GAPS/CGAL_Arr_2_for_GAPS_traits.h>

// Curved-kernel Conics:
#elif BENCH_TRAITS == CK_CIRCLE_TRAITS
#elif BENCH_TRAITS == CK_CONIC_TRAITS
#else
#error No traits (TRAITS) specified!
#endif

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_default_dcel.h>

#include <CGAL/IO/Arr_iostream.h>
// #include <CGAL/IO/Arr_Window_stream.h>
#if defined(POSTSCRIPT_SUPPORTED)
#include <CGAL/IO/Arr_Postscript_file_stream.h>
#endif

#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_naive_point_location.h>

#if defined(LANDMARK_SUPPORTED)
#include <CGAL/Arr_nearest_neighbor.h>
#include <CGAL/Arr_landmarks_point_location.h>
#endif

#if defined(TRIANGLE_SUPPORTED)
#include <CGAL/Arr_triangle_point_location.h>
#endif

#include <CGAL/Benchmark/Benchmark.hpp>

#if defined(USE_CGAL_WINDOW)
#include <CGAL/IO/Window_stream.h>
#else
#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h>
#endif

#include <cstdlib>
#include <iostream>
#include <list>

#include "Option_parser.hpp"

// Traits:
#if BENCH_TRAITS == SEGMENT_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Base_traits;
#define TRAITS_TYPE "Segments"

#elif BENCH_TRAITS == LEDA_SEGMENT_TRAITS
typedef CGAL::Arr_leda_segment_traits_2<Kernel>         Base_traits;
#define TRAITS_TYPE "Leda Segments"

#elif BENCH_TRAITS == NON_CACHING_SEGMENT_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Base_traits;
#define TRAITS_TYPE "Non Caching Segments"

#elif BENCH_TRAITS == POLYLINE_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Base_traits;
#define TRAITS_TYPE "Polylines"

#elif BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Base_traits;
#define TRAITS_TYPE "Non Caching Polylines"

#elif BENCH_TRAITS == LEDA_CONIC_TRAITS
typedef CGAL::Arr_conic_traits_2<Kernel>                Base_traits;
#define TRAITS_TYPE "Conics"

#elif BENCH_TRAITS == CORE_CONIC_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;

#if BENCH_KERNEL == CARTESIAN_KERNEL
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
#elif BENCH_KERNEL == SIMPLE_CARTESIAN_KERNEL
typedef CGAL::Simple_cartesian<Rational>                Rat_kernel;
typedef CGAL::Simple_cartesian<Algebraic>               Alg_kernel;
#endif

typedef CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>
                                                        Base_traits;
#define TRAITS_TYPE "Conics"

// Exacus Conics:
#elif BENCH_TRAITS == EXACUS_CONIC_TRAITS
typedef CnX::Conic_sweep_traits_2<Arithmetic_traits>    CST;
typedef CnX::Conic_segment_2< CST>                      Input_segment;
typedef SoX::CGAL_Arr_2_for_GAPS_traits< Input_segment, CST> Base_traits;
#define TRAITS_TYPE "Exacus Conics"

#else
#error No traits (TRAITS) specified!
#endif

// Arrangement types:
#if defined(CGAL_TRAITS_COUNTING)
typedef CGAL::Arr_counting_traits_2<Base_traits>        Traits;
#elif defined(CGAL_TRAITS_TRACING)
typedef CGAL::Arr_tracing_traits_2<Base_traits>         Traits;
#else
typedef  Base_traits                                    Traits;
#endif
typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arr;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef std::list<Curve_2>                              Curve_list;
typedef std::list<Point_2>                              Point_list;
typedef Arr::Halfedge_iterator                          Arr_halfedge_iterator;

// Point location strategies:
typedef CGAL::Arr_trapezoid_ric_point_location<Arr>     Trap_point_location;
typedef CGAL::Arr_naive_point_location<Arr>             Naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Arr>   Walk_point_location;

#if defined(LANDMARK_SUPPORTED)
typedef CGAL::Arr_nearest_neighbor<Arr>                 Nearest_neighbor;
typedef CGAL::Arr_landmarks_point_location<Arr,Nearest_neighbor>
  Landmarks_point_location;
#endif

#if defined(TRIANGLE_SUPPORTED)
typedef CGAL::Arr_triangle_point_location<Arr>          Triangle_point_location;
#endif

typedef Option_parser::Type_code        Type_code;
typedef Option_parser::Strategy_code    Strategy_code;

// Window stream:
#if defined(USE_CGAL_WINDOW)
typedef CGAL::Window_stream Window_stream;
#else
typedef CGAL::Qt_widget Window_stream;
QApplication * App;
#endif

// Readers:
// Conic reader:
#if BENCH_TRAITS == SEGMENT_TRAITS || BENCH_TRAITS == NON_CACHING_SEGMENT_TRAITS
#include "Segment_reader.hpp"

#elif BENCH_TRAITS == LEDA_CONIC_TRAITS || BENCH_TRAITS == CORE_CONIC_TRAITS ||\
      BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS || \
      BENCH_TRAITS == EXACUS_CONIC_TRAITS
#include "Conic_reader.hpp"

// Polyline reader:
#elif BENCH_TRAITS == POLYLINE_TRAITS || \
      BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
#include "Polyline_reader.hpp"

#else
#error No traits (TRAITS) specified
#endif

namespace CGAL {

inline std::ostream & operator<<(std::ostream & os, const Arr::Vertex & vertex)
{
  os << vertex.point();
  return os;
}

} //namespace CGAL

/*! Exporter of an arrangement to a window stream */
inline Window_stream & operator<<(Window_stream & ws, Arr & arr)
{
  Arr::Edge_iterator ei;
  ws << CGAL::IO::blue();
  for (ei = arr.edges_begin(); ei != arr.edges_end(); ++ei)
    ws << (*ei).curve();
  Arr::Vertex_iterator vi;
  ws << CGAL::IO::red();
  for (vi = arr.vertices_begin(); vi != arr.vertices_end(); ++vi)
    ws << (*vi).point();
  return ws;
}

/*! */
class Basic_arr {
public:
  /*! Constructor */
  Basic_arr() :
    m_filename(0), m_verbose_level(0), m_postscript(false),
    m_bbox(0.0, 0.0, 0.0, 0.0),
    m_width(DEF_WIN_WIDTH), m_height(DEF_WIN_HEIGHT)
  { }

  virtual ~Basic_arr() {}

  /*! */
  //virtual void op() = 0;

  /*! */
  int init()
  {
#if BENCH_TRAITS == SEGMENT_TRAITS || BENCH_TRAITS == NON_CACHING_SEGMENT_TRAITS
    Segment_reader<Traits> reader;

#elif BENCH_TRAITS == POLYLINE_TRAITS || \
      BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
    Polyline_reader<Traits> reader;

#elif BENCH_TRAITS == LEDA_CONIC_TRAITS || BENCH_TRAITS == CORE_CONIC_TRAITS ||\
      BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS || \
      BENCH_TRAITS == EXACUS_CONIC_TRAITS
    Conic_reader<Traits> reader;

#else
#error "Run out of options!"
#endif

    int rc = reader.read_data(m_filename, std::back_inserter(m_curve_list),
                              m_bbox);
    if (rc < 0) return rc;
    if (m_verbose_level > 0)
      std::cout << m_curve_list.size() << " curves" << std::endl;

    return 0;
  }

  /* */
  void clean()
  {
    m_curve_list.clear();
#if defined(CGAL_TRAITS_COUNTING)
    if (m_verbose_level > 0) std::cout << m_traits;
#endif
  }
  void sync(){}

  void set_file_name(const char * filename, int file_index=0)
  { m_filename = filename; }
  void set_verbose_level(const unsigned int level) { m_verbose_level = level; }
  void set_postscript(const bool postscript) { m_postscript = postscript; }

  void set_width(unsigned int width) { m_width = width; }
  void set_height(unsigned int height) { m_height = height; }

protected:
  const char * m_filename;
  Curve_list m_curve_list;
  unsigned int m_verbose_level;
  bool m_postscript;
  CGAL::Bbox_2 m_bbox;

  int m_width, m_height;
  float m_x0, m_x1, m_y0, m_y1;

  Traits m_traits;
};

/*!  This class is supplied to the bench as a templatwe parameter.
 * it defines the operations we want to measure.
 * in this case we measure only the incremental construction
 * of the arrangement.
*/
template <class Strategy>
class Increment_arr : public Basic_arr {
public:
  // op is the main function that must be implemented in each benchable class
  // that is supplied to the Bench class.
  // the other functions that should be implemented are inherited from Basic_arr
  // This (op) is the main method, that is different from all other classes
  // inherite from Basic_arr
  void op()
  {
    Arr arr(&m_traits);
    Curve_list::const_iterator i;
    Strategy strategy(arr);
    for (i = m_curve_list.begin(); i != m_curve_list.end(); ++i)
      insert(arr, *i, strategy);

    if (m_verbose_level > 0) {      //print to cout
      if (m_verbose_level > 1) {
        if (!arr.is_valid()) std::cerr << "map invalid!" << std::endl;
      }
      std::cout << "# of vertices: " << arr.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << arr.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << arr.number_of_faces() << std::endl;
    }

#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
    // Clear the exacus cache:
    CST::Algebraic_curve_2::clear_cache();
    CST::Point_2::clear_cache();
#endif
  }
};

typedef Increment_arr<Trap_point_location>      Trap_inc_arr;
typedef cb::Benchmark<Trap_inc_arr>             Trap_inc_arr_bench;

typedef Increment_arr<Naive_point_location>     Naive_inc_arr;
typedef cb::Benchmark<Naive_inc_arr>            Naive_inc_arr_bench;

typedef Increment_arr<Walk_point_location>      Walk_inc_arr;
typedef cb::Benchmark<Walk_inc_arr>             Walk_inc_arr_bench;

#if defined(LANDMARK_SUPPORTED)
typedef Increment_arr<Landmarks_point_location> Landmarks_inc_arr;
typedef cb::Benchmark<Landmarks_inc_arr>        Landmarks_inc_arr_bench;
#endif

#if defined(TRIANGLE_SUPPORTED)
typedef Increment_arr<Triangle_point_location>  Triangle_inc_arr;
typedef cb::Benchmark<Triangle_inc_arr>         Triangle_inc_arr_bench;
#endif

/*! */
class Aggregate_arr : public Basic_arr {
public:
  void op()
  {
    if (m_verbose_level > 0) std::cout << "Inserting Aggregate" << std::endl;
    Arr arr(&m_traits);
    insert(arr, m_curve_list.begin(), m_curve_list.end());
    if (m_verbose_level > 0) {
      if (m_verbose_level > 1) {
        if (!arr.is_valid()) std::cerr << "map invalid!" << std::endl;
      }
      std::cout << "# of vertices: " << arr.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << arr.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << arr.number_of_faces() << std::endl;
    }
#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
    // Clear the exacus cache:
    CST::Algebraic_curve_2::clear_cache();
    CST::Point_2::clear_cache();
#endif
  }
};

/*! */
typedef cb::Benchmark<Aggregate_arr>            Agg_arr_bench;

/*! */
class Display_arr : public Basic_arr {
public:
  /*! */
  void op()
  {
    Arr arr(&m_traits);
    insert(arr, m_curve_list.begin(), m_curve_list.end());
    if (m_verbose_level > 0) {
      if (m_verbose_level > 1) {
        if (!arr.is_valid()) std::cerr << "map invalid!" << std::endl;
      }
      std::cout << "# of vertices: " << arr.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << arr.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << arr.number_of_faces() << std::endl;

      if (m_verbose_level > 1) {
        std::copy(arr.vertices_begin(), arr.vertices_end(),
                  std::ostream_iterator<Arr::Vertex>(std::cout, "\n"));
      }
    }

#if defined(USE_CGAL_WINDOW)
    m_window->set_flush(0);
    (*m_window) << arr;
    m_window->set_flush(1);
    m_window->flush();
#else
    m_window->lock();
    *m_window << CGAL::BackgroundColor(CGAL::white()) << CGAL::red();
    (*m_window) << arr;
    m_window->unlock();
    App->flush();
#endif

#if defined(POSTSCRIPT_SUPPORTED)
    if (m_postscript) {
      // Print to Postscript file:
      std::cout << "Print to Postscript file" << std::endl;
      CGAL::Postscript_file_stream ps_stream(m_width, m_height ,"arr.ps");
      ps_stream.init(m_x0, m_x1, m_y0);
      // CGAL::cgalize(ps_stream);
      ps_stream.set_line_width(1);
      CGAL::Arr_drawer<Arr, CGAL::Postscript_file_stream> drawer(ps_stream);
      // drawer.draw_faces(arr.faces_begin(), arr.faces_end());
      ps_stream << CGAL::blue();
      drawer.draw_halfedges(arr.halfedges_begin(), arr.halfedges_end());
      ps_stream << CGAL::red();
      drawer.draw_vertices(arr.vertices_begin(), arr.vertices_end());

      // draw_arr(arr, drawer, ps_stream);
      // ps_stream << arr;
    }
#endif

#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
    // Clear the exacus cache:
    CST::Algebraic_curve_2::clear_cache();
    CST::Point_2::clear_cache();
#endif
}

  /*! */
  int init()
  {
    int rc = Basic_arr::init();
    if (rc < 0) return rc;

    float x_range = m_bbox.xmax() - m_bbox.xmin();
    float y_range = m_bbox.ymax() - m_bbox.ymin();
    float height = (y_range * m_width) / x_range;

    float min_range = (x_range < y_range) ? x_range : y_range;
    float x_margin = min_range / 4;
    float y_margin = (height * x_margin) / m_width;

    m_x0 = m_bbox.xmin() - x_margin;
    m_x1 = m_bbox.xmax() + x_margin;
    m_y0 = m_bbox.ymin() - y_margin;

    m_height = static_cast<int>(height);

#if defined(USE_CGAL_WINDOW)
    m_window = new Window_stream(m_width, m_height);
    if (!m_window) return -1;
    m_window->init(m_x0, m_x1, m_y0);   // logical window size

    m_window->set_redraw(&Display_arr::redraw);
    m_window->set_mode(leda_src_mode);
    m_window->set_node_width(3);
    m_window->set_point_style(leda_cross_point);
    m_window->set_line_width(1);
    m_window->display(leda_window::center, leda_window::center);
#else
    m_y1 = m_bbox.ymax() + y_margin;

    m_window = new Window_stream();
    if (!m_window) return -1;
    App->setMainWidget(m_window);
    m_window->resize(m_width, m_height);
    m_window->set_window(m_x0, m_x1, m_y0, m_y1);   // logical window size
    m_window->setLineWidth(1);
    m_window->setPointSize(3);
    m_window->show();
#endif
    return 0;
  }

  /*! */
  void clean()
  {
    Basic_arr::clean();
    delete m_window;
  }

private:
  /*! */
#if defined(USE_CGAL_WINDOW)
  static void
    redraw(leda_window * wp, double x0, double y0, double x1, double y1)
  { wp->flush_buffer(x0,y0,x1,y1); }
#endif

  Window_stream * m_window;
};

typedef cb::Benchmark<Display_arr>              Dis_arr_bench;

/*! */
template <class Bench_inst, class Benchable>
void run_bench(Bench_inst & bench_inst, Benchable & benchable,
               const char * fullname,
               int samples, int iterations, unsigned int verbose_level,
               bool postscript = false,
               const char * fullname2 = 0)
{
  //set some variable
  benchable.set_file_name(fullname);
  benchable.set_verbose_level(verbose_level);
  benchable.set_postscript(postscript);
  if (fullname2) benchable.set_file_name(fullname2, 1);

  if (samples > 0) bench_inst.set_samples(samples);
  else if (iterations > 0) bench_inst.set_iterations(iterations);

  //opertor () in the Bench - does all the work !
  bench_inst();
}

/* */
int main(int argc, char * argv[])
{
  Option_parser option_parser;
  try {
    option_parser(argc, argv);
  } catch(Option_parser::Generic_option_exception & /* e */) {
    return 0;
  } catch(std::exception & e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  unsigned int verbose_level = option_parser.get_verbose_level();
  unsigned int type_mask = option_parser.get_type_mask();
  unsigned int strategy_mask = option_parser.get_strategy_mask();
  unsigned int samples = option_parser.get_samples();
  int iterations = option_parser.get_iterations();
  int seconds = option_parser.get_seconds();
  bool print_header = option_parser.is_print_header();
  int nameLength = option_parser.get_name_length();
  unsigned int files_number = option_parser.get_number_files();
  for (unsigned int i = 0; i < files_number; ++i) {
    if (verbose_level > 0) {
      std::cout << "filename: " << option_parser.get_file_name(i).c_str()
                << std::endl
                << "fullname: " << option_parser.get_full_name(i).c_str()
                << std::endl;
    }
  }

  bool postscript = option_parser.get_postscript();
  if (postscript) samples = 1;

  const char * file_name = option_parser.get_file_name(0).c_str();
  const char * full_name = option_parser.get_full_name(0).c_str();
  if (!file_name || !full_name) return -1;

  cb::Benchmark_base::set_name_length(nameLength);
  if (print_header) cb::Benchmark_base::print_header();
  //end parameters from command line

  //start bench
  if (verbose_level > 0) std::cout << "start bench " << std::endl;

  //general definitions
  Type_code type_code;           // operation we want to apply
  Strategy_code strategy_code;   // point location strategy

  if (verbose_level > 0) {
    std::cout << "type_mask = "<< type_mask << std::endl;
    std::cout << "strategy_mask = " << strategy_mask  << std::endl;
  }

  // Construct Incrementaly  (only if type_code == incremental)
  type_code = Option_parser::TYPE_INCREMENT;
  if (type_mask & (0x1 << type_code)) {
    if (verbose_level > 0) std::cout << "TYPE_INCREMENT " << std::endl;

#if defined(RIC_SUPPORTED)
    // Ric point location:
    strategy_code = Option_parser::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_code)) {
      std::string name =
        std::string(option_parser.get_type_name(type_code)) + " " +
        std::string(option_parser.get_strategy_name(strategy_code)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + " (" + std::string(filename[0]) + ")";
      Trap_inc_arr_bench bench_inst(name, seconds, false);
      Trap_inc_arr & benchable = bench_inst.get_benchable();
      run_bench<Trap_inc_arr_bench, Trap_inc_arr>(bench_inst, benchable,
                                                  full_name,
                                                  samples, iterations,
                                                  verbose_level, postscript);
    }
#endif

    // Naive point location:
    strategy_code = Option_parser::STRATEGY_NAIVE;
    if (strategy_mask & (0x1 << strategy_code)) {
      std::string name =
        std::string(option_parser.get_type_name(type_code)) + " " +
        std::string(option_parser.get_strategy_name(strategy_code)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + " (" + std::string(file_name) + ")";
      Naive_inc_arr_bench bench_inst(name, seconds, false);
      Naive_inc_arr & benchable = bench_inst.get_benchable();
      run_bench<Naive_inc_arr_bench, Naive_inc_arr>(bench_inst, benchable,
                                                    full_name,
                                                    samples, iterations,
                                                    verbose_level, postscript);
    }

    // Walk point location:
    strategy_code = Option_parser::STRATEGY_WALK;
    if (strategy_mask & (0x1 << strategy_code)) {
      std::string name =
        std::string(option_parser.get_type_name(type_code)) + " " +
        std::string(option_parser.get_strategy_name(strategy_code)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + " (" + std::string(file_name) + ")";
      Walk_inc_arr_bench bench_inst(name, seconds, false);
      Walk_inc_arr & benchable = bench_inst.get_benchable();
      run_bench<Walk_inc_arr_bench, Walk_inc_arr>(bench_inst, benchable,
                                                  full_name,
                                                  samples, iterations,
                                                  verbose_level, postscript);
    }

#if defined(LANDMARK_SUPPORTED)
    // Landmarks point location:
    strategy_code = Option_parser::STRATEGY_LANDMARKS;
    if (strategy_mask & (0x1 << strategy_code)) {
      std::string name =
        std::string(option_parser.get_type_name(type_code)) + " " +
        std::string(option_parser.get_strategy_name(strategy_code)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + " (" + std::string(file_name) + ")";
      Landmarks_inc_arr_bench bench_inst(name, seconds, false);
      Landmarks_inc_arr & benchable = bench_inst.get_benchable();
      run_bench<Landmarks_inc_arr_bench, Landmarks_inc_arr>(bench_inst,
                                                            benchable,
                                                            full_name,
                                                            samples,
                                                            iterations,
                                                            verbose_level,
                                                            postscript);
    }
#endif

#if defined(TRIANGLE_SUPPORTED)
    // Triangle point location:
    strategy_code = Option_parser::STRATEGY_TRIANGLE;
    if (strategy_mask & (0x1 << strategy_code)) {
      std::string name =
        std::string(option_parser.get_type_name(type_code)) + " " +
        std::string(option_parser.get_strategy_name(strategy_code)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE +" " + " (" + std::string(file_name) + ")";
      Triangle_inc_arr_bench bench_inst(name, seconds, false);
      Triangle_inc_arr & benchable = bench_inst.get_benchable();
      run_bench<Triangle_inc_arr_bench,Triangle_inc_arr>(bench_inst, benchable,
                                                         full_name,
                                                         samples, iterations,
                                                         verbose_level,
                                                         postscript);
    }
#endif
  }

  // Construct Aggregately
  type_code = Option_parser::TYPE_AGGREGATE;
  if (type_mask & (0x1 << type_code)) {
    if (verbose_level > 0) std::cout << "TYPE_AGGREGATE " << std::endl;
    std::string name =
      std::string(option_parser.get_type_name(type_code)) + " " +
      TRAITS_TYPE + " " + KERNEL_TYPE + " " +
      NUMBER_TYPE + " " + " (" + std::string(file_name) + ")";
    Agg_arr_bench bench_inst(name, seconds, false);
    Aggregate_arr & benchable = bench_inst.get_benchable();
    run_bench<Agg_arr_bench,Aggregate_arr>(bench_inst, benchable, full_name,
                                           samples, iterations,
                                           verbose_level, postscript);
  }

  // Construct and Display
  type_code = Option_parser::TYPE_DISPLAY;

  if (type_mask & (0x1 << type_code)) {
    unsigned int width = option_parser.get_width();
    unsigned int height = option_parser.get_height();

#if !defined(USE_CGAL_WINDOW)
    QApplication app(argc, argv);
    App = &app;
#endif
    if (verbose_level > 0) std::cout << "TYPE_DISPLAY " << std::endl;
    std::string name =
      std::string(option_parser.get_type_name(type_code)) + " " +
      TRAITS_TYPE + " " + KERNEL_TYPE + " " +
      NUMBER_TYPE + " " + " (" + std::string(file_name) + ")";
    Dis_arr_bench bench_inst(name, seconds, false);
    Display_arr & benchable = bench_inst.get_benchable();
    benchable.set_width(width);
    benchable.set_height(height);
    run_bench<Dis_arr_bench,Display_arr>(bench_inst, benchable, full_name,
                                         samples, iterations,
                                         verbose_level, postscript);
  }

  return 0;
}
