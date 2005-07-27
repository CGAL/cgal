#include <CGAL/config.h>
#include <CGAL/basic.h>
#include <CGAL/IO/Color.h>

#include "bench_config.h"
#include "number_type.h"

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

// Kernel:
#if BENCH_KERNEL == LEDA_KERNEL
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>

#elif BENCH_KERNEL == MY_KERNEL
#include <CGAL/Arr_segment_traits_leda_kernel_2.h>

#elif BENCH_KERNEL == SIMPLE_CARTESIAN_KERNEL
#include <CGAL/Simple_cartesian.h>

#elif BENCH_KERNEL == CARTESIAN_KERNEL
#include <CGAL/Cartesian.h>

#else
#error No kernel (KERNEL) specified!
#endif

#if BENCH_KERNEL == LEDA_KERNEL || BENCH_KERNEL == MY_KERNEL
#if defined(USE_CGAL_WINDOW)
#include <LEDA/rat_window.h>
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

#elif BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>

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

// This does not belong here
#if BENCH_NT == DOUBLE_NT
#include <ECG/Root_of/double.h>
#elif BENCH_NT == LEDA_REAL_NT
#include <ECG/Root_of/leda_real.h>
#elif BENCH_NT == QUOTIENT_MP_FLOAT_NT
#include <ECG/Root_of/CGAL_Quotient.h>
#elif BENCH_NT == GMPZ_NT
#include <ECG/Root_of/gmpxx.h>
#elif BENCH_NT == GMPQ_NT
#include <ECG/Root_of/gmpxx.h>
#elif BENCH_NT == CGAL_GMPQ_NT
#include <ECG/Root_of/CGAL_Gmpq.h>
#elif BENCH_NT == LAZY_LEDA_RAT_NT
#include <ECG/Root_of/CGAL_Lazy_exact_nt.h>
#elif BENCH_NT == LAZY_CGAL_GMPQ_NT
#include <ECG/Root_of/CGAL_Lazy_exact_nt.h>
#include <ECG/Root_of/gmpxx.h>
#elif BENCH_NT == LAZY_GMPZ_NT
#include <ECG/Root_of/CGAL_Lazy_exact_nt.h>
#include <ECG/Root_of/gmpxx.h>
#elif BENCH_NT == LAZY_QUOTIENT_MP_FLOAT_NT
#include <ECG/Root_of/CGAL_Lazy_exact_nt.h>
#include <ECG/Root_of/CGAL_Quotient.h>
#elif BENCH_NT == CORE_EXPR_NT
#include <ECG/Root_of/CORE_Expr.h>
#endif

#include <ECG/Root_of/Root_of_2.h>

#include <ECG/Circular_kernel.h>
#include <ECG/Circular_arc_traits.h>
#include <ECG/Circular_arc_traits_tracer.h> 
#include <ECG/IO/Qt_widget_circular_arc_2.h>
#include <ECG/IO/Qt_widget_circular_arc_endpoint_2.h>

#elif BENCH_TRAITS == CK_CONIC_TRAITS
// #include <CGAL/Random.h>
#include <ECG/Synaps_kernel.h>
#include <ECG/Conic_kernel.h>
#include <ECG/Conic_arc_traits.h>
#include <ECG/Conic_arc_traits_tracer.h> 
#include <ECG/IO/Qt_widget_conic_arc_2.h>
#include <ECG/IO/Qt_widget_conic_arc_endpoint_2.h>

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

#include <CGAL/Bench.h>
#include <CGAL/Bench.C>
#include <CGAL/Bench_parse_args.h>
#include <CGAL/Bench_parse_args.C>

#if defined(USE_CGAL_WINDOW)
#include <CGAL/IO/Window_stream.h>
#else
#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h>
#endif

#include <stdlib.h>
#include <iostream>
#include <list>

#include "Point_reader.h"       //for locate and not insert

// Readers:
// Conic reader:
#if BENCH_TRAITS == SEGMENT_TRAITS || BENCH_TRAITS == NON_CACHING_SEGMENT_TRAITS
#include "Segment_reader.h"

#elif BENCH_TRAITS == LEDA_CONIC_TRAITS || BENCH_TRAITS == CORE_CONIC_TRAITS ||\
      BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS || \
      BENCH_TRAITS == EXACUS_CONIC_TRAITS
#include "Conic_reader.h"

// Polyline reader:
#elif BENCH_TRAITS == POLYLINE_TRAITS || \
      BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
#include "Polyline_reader.h"

#else
#error No traits (TRAITS) specified
#endif

// Kernel:
#if BENCH_KERNEL == LEDA_KERNEL
typedef CGAL::leda_rat_kernel_traits                    Kernel;
#define KERNEL_TYPE "Leda"

#elif BENCH_KERNEL == MY_KERNEL
typedef CGAL::Arr_segment_traits_leda_kernel_2          Kernel;
#define KERNEL_TYPE "Partial Leda"

#elif BENCH_KERNEL == SIMPLE_CARTESIAN_KERNEL
typedef CGAL::Simple_cartesian<WNT>                     Kernel;
#define KERNEL_TYPE "Simple Cartesian"

#elif BENCH_KERNEL == CARTESIAN_KERNEL
typedef CGAL::Cartesian<WNT>                            Kernel;
#define KERNEL_TYPE "Cartesian"

#else
#error No kernel (KERNEL) specified!
#endif

// Traits:
#if BENCH_TRAITS == SEGMENT_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#define TRAITS_TYPE "Segments"

#elif BENCH_TRAITS == LEDA_SEGMENT_TRAITS
typedef CGAL::Arr_leda_segment_traits_2<Kernel>         Traits;
#define TRAITS_TYPE "Leda Segments"

#elif BENCH_TRAITS == NON_CACHING_SEGMENT_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Traits;
#define TRAITS_TYPE "Non Caching Segments"

#elif BENCH_TRAITS == POLYLINE_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Traits;
#define TRAITS_TYPE "Polylines"

#elif BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Traits;
#define TRAITS_TYPE "Non Caching Polylines"

#elif BENCH_TRAITS == LEDA_CONIC_TRAITS
typedef CGAL::Arr_conic_traits_2<Kernel>                Traits;
#define TRAITS_TYPE "Conics"

#elif BENCH_TRAITS == CORE_CONIC_TRAITS
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>
                                                        Traits;
#define TRAITS_TYPE "Conics"

// Exacus Conics:
#elif BENCH_TRAITS == EXACUS_CONIC_TRAITS
typedef CnX::Conic_sweep_traits_2<Arithmetic_traits>    CST;
typedef CnX::Conic_segment_2< CST>                      Input_segment;
typedef SoX::CGAL_Arr_2_for_GAPS_traits< Input_segment, CST> Traits;
#define TRAITS_TYPE "Exacus Conics"

// Curved-kernel Circle:
#elif BENCH_TRAITS == CK_CIRCLE_TRAITS
typedef ECG::Curved_kernel<Kernel>                      Curved_k;
typedef ECG::Circular_arc_traits<Curved_k>              Traits;
#define TRAITS_TYPE "Curved Kernel Circles"

// Curved-kernel Conics:
#elif BENCH_TRAITS == CK_CONIC_TRAITS
typedef ECG::Algebraic::Synaps_kernel<RT,FT>            Algebraic_k;
typedef ECG::Curved_kernel<Kernel,Algebraic_k>          Curved_k;
typedef ECG::Conic_arc_traits<Curved_k>                 Traits;
#define TRAITS_TYPE "Curved Kernel Conics"

#else
#error No traits (TRAITS) specified!
#endif

// Arrangement types:
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

typedef CGAL::Bench_parse_args::TypeId                  Type_id;
typedef CGAL::Bench_parse_args::StrategyId              Strategy_id;
typedef CGAL::Bench_parse_args::FormatId                Format_id;

// Window stream:
#if defined(USE_CGAL_WINDOW)
typedef CGAL::Window_stream Window_stream;
#else
typedef CGAL::Qt_widget Window_stream;
QApplication * App;
#endif

CGAL_BEGIN_NAMESPACE

inline std::ostream & operator<<(std::ostream & os, const Arr::Vertex & vertex)
{
  os << vertex.point();
  return os;
}

CGAL_END_NAMESPACE

#if BENCH_TRAITS == POLYLINE_TRAITS || \
    BENCH_TRAITS == NON_CACHING_POLYLINE_TRAITS
/*! Exporter operator for a polyline to a window stream */
template <class T_SegmentTraits>
Window_stream & operator<<(Window_stream & ws,
                           const CGAL::Polyline_2<T_SegmentTraits> & cv)
{
  typename CGAL::Polyline_2<T_SegmentTraits>::const_iterator iter = cv.begin();
  for (; iter != cv.end(); ++iter) ws << *iter;
  return ws;
}
#endif

/*! Exporter of an arrangement to a window stream */
inline Window_stream & operator<<(Window_stream & os, Arr & arr)
{
  Arr::Edge_iterator ei;
  os << CGAL::BLUE;
  for (ei = arr.edges_begin(); ei != arr.edges_end(); ++ei)
    os << (*ei).curve();
  Arr::Vertex_iterator vi;
  os << CGAL::RED;
  for (vi = arr.vertices_begin(); vi != arr.vertices_end(); ++vi)
    os << (*vi).point();
  return os;
}

/*! */
class Basic_arr {
public:

  /*! Constructor */
  Basic_arr() :
    m_filename(0), m_verbose(false), m_postscript(false),
    m_bbox(0.0, 0.0, 0.0, 0.0),
    m_width(1024), m_height(1024)
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
                              m_format, m_bbox);
    if (rc < 0) return rc;
    if (m_verbose) std::cout << m_curve_list.size() << " curves" << std::endl;

    return 0;
  }
    
  /* */
  void clean() { m_curve_list.clear(); }
  void sync(){}

  void set_format(Format_id fmt) { m_format = fmt; }
  void set_file_name(const char * filename, int file_index=0) 
  { m_filename = filename; }
  void set_verbose(const bool verbose) { m_verbose = verbose; }
  void set_postscript(const bool postscript) { m_postscript = postscript; }

protected:
  const char * m_filename;
  Curve_list m_curve_list;
  bool m_verbose;
  bool m_postscript;
  Format_id m_format;
  CGAL::Bbox_2 m_bbox;

  int m_width, m_height;
  float m_x0, m_x1, m_y0, m_y1;
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
    Arr arr;
    Curve_list::const_iterator i;
    Strategy strategy(arr);
    for (i = m_curve_list.begin(); i != m_curve_list.end(); ++i)
      insert(arr, *i, strategy);
    
    if (m_verbose) {      //print to cout
      if (!arr.is_valid()) std::cerr << "map invalid!" << std::endl;
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
typedef CGAL::Bench<Trap_inc_arr>               Trap_inc_arr_bench;

typedef Increment_arr<Naive_point_location>     Naive_inc_arr;
typedef CGAL::Bench<Naive_inc_arr>              Naive_inc_arr_bench;

typedef Increment_arr<Walk_point_location>      Walk_inc_arr;
typedef CGAL::Bench<Walk_inc_arr>               Walk_inc_arr_bench;

#if defined(LANDMARK_SUPPORTED)
typedef Increment_arr<Landmarks_point_location> Landmarks_inc_arr;
typedef CGAL::Bench<Landmarks_inc_arr>          Landmarks_inc_arr_bench;
#endif

#if defined(TRIANGLE_SUPPORTED)
typedef Increment_arr<Triangle_point_location>  Triangle_inc_arr;
typedef CGAL::Bench<Triangle_inc_arr>           Triangle_inc_arr_bench;
#endif

/*! */
class Aggregate_arr : public Basic_arr {
public:
  void op()
  {
    if (m_verbose) std::cout << "Inserting Aggregate" << std::endl;
    Arr arr;
    insert(arr, m_curve_list.begin(), m_curve_list.end());
    if (m_verbose) {
      if (!arr.is_valid()) std::cerr << "map invalid!" << std::endl;
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
typedef CGAL::Bench<Aggregate_arr>               Agg_arr_bench;

/*! */
class Display_arr : public Basic_arr {
public:
  /*! */
  void op()
  {
    Arr arr;
    insert(arr, m_curve_list.begin(), m_curve_list.end());
    if (m_verbose) {
      if (!arr.is_valid()) std::cerr << "map invalid!" << std::endl;
      std::cout << "# of vertices: " << arr.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << arr.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << arr.number_of_faces() << std::endl;

      std::copy(arr.vertices_begin(), arr.vertices_end(),
                std::ostream_iterator<Arr::Vertex>(std::cout, "\n"));
    }

#if defined(USE_CGAL_WINDOW)
    m_window->set_flush(0);
    (*m_window) << arr;
    m_window->set_flush(1);
    m_window->flush();
#else
    m_window->lock();
    *m_window << CGAL::BackgroundColor(CGAL::WHITE) << CGAL::RED;
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
      ps_stream << CGAL::BLUE;
      drawer.draw_halfedges(arr.halfedges_begin(), arr.halfedges_end());
      ps_stream << CGAL::RED;
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
    m_window->setLineWidth(4);
    m_window->setPointSize(4);
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
  static void
    redraw(leda_window * wp, double x0, double y0, double x1, double y1) 
  { wp->flush_buffer(x0,y0,x1,y1); }

  Window_stream * m_window;
};

typedef CGAL::Bench<Display_arr>                Dis_arr_bench;

/*! */
template <class Bench_inst, class Benchable>
void run_bench(Bench_inst & benchInst, Benchable & benchable,
               const char * fullname, Format_id format,
               int samples, int iterations, bool verbose,
               bool postscript = false, 
               const char * fullname2 = 0)
{
  //set some variable
  benchable.set_format(format);
  benchable.set_file_name(fullname);
  benchable.set_verbose(verbose);
  benchable.set_postscript(postscript);
  if (fullname2) benchable.set_file_name(fullname2, 1);
  
  if (samples > 0) benchInst.set_samples(samples);
  else if (iterations > 0) benchInst.set_iterations(iterations);
  
  //opertor () in the Bench - does all the work !
  benchInst();

  //cout
  if (verbose) std::cout << "(" << benchInst.get_samples() << ") "
                         << std::endl;
}

/* */
int main(int argc, char * argv[])
{
  // get arguments:
  if ((argc < 2)) {
    std::cout << "Usage: "<< argv[0] <<" [OPTION] input_filename"<<std::endl;
    std::cout << "Try `"<<argv[0]<<" -h` for more information"<<std::endl;
    exit(-1);
  }
  
  //get arguments from command line
  CGAL::Bench_parse_args parse_args(argc, argv);
  int rc = parse_args.parse();
  if (rc > 0) return 0;
  if (rc < 0) return rc;
  
  bool verbose = parse_args.get_verbose();
  unsigned int type_mask = parse_args.get_type_mask();
  unsigned int strategy_mask = parse_args.get_strategy_mask();
  Format_id format = parse_args.get_input_format();
  int samples = parse_args.get_samples();
  int iterations = parse_args.get_iterations();
  int seconds = parse_args.get_seconds();
  bool printHeader = parse_args.get_print_header();
  int nameLength = parse_args.get_name_length();
  unsigned int files_number = parse_args.get_files_number();
  //  std::cout << "files_number: " << files_number << std::endl;
  const char * filename[MAX_FILES_NUMBER];
  const char * fullname[MAX_FILES_NUMBER];
  for (unsigned int i = 0; i < files_number; ++i) {
    filename[i] = parse_args.get_file_name(i);
    fullname[i] = parse_args.get_full_name(i);
    if (verbose) {
      std::cout << "filename: " << filename[i] << std::endl;
      std::cout << "fullname: " << fullname[i] << std::endl;
    }
  }

  bool postscript = parse_args.get_postscript();
  if (postscript) samples = 1;

  if (!filename[0] || !fullname[0]) return -1;

  CGAL::Bench_base::set_name_length(nameLength);
  if (printHeader) CGAL::Bench_base::print_header();
  //end parameters from command line

  //start bench
  if (verbose) std::cout << "start bench " << std::endl;

  //general definitions
  Type_id type_id;           //operation we want to apply
  Strategy_id strategy_id;   //point location strategy 
  
  if (verbose) {
    std::cout << "type_mask = "<< type_mask << std::endl;
    std::cout << "strategy_mask = " << strategy_mask  << std::endl;
  }

  // Construct Incrementaly  (only if type_id == incremental)
  type_id = CGAL::Bench_parse_args::TYPE_INCREMENT;
  if (type_mask & (0x1 << type_id)) {
    if (verbose) std::cout << "TYPE_INCREMENT " << std::endl;
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
        NUMBER_TYPE + " " + " (" + std::string(filename[0]) + ")";
      Trap_inc_arr_bench benchInst(name, seconds, false);
      Trap_inc_arr & benchable = benchInst.get_benchable();
      run_bench<Trap_inc_arr_bench,Trap_inc_arr>(benchInst, benchable,
                                                 fullname[0], format,
                                                 samples, iterations,
                                                 verbose, postscript);
    }
    
    // Naive point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + " (" + std::string(filename[0]) + ")";
      Naive_inc_arr_bench benchInst(name, seconds, false);
      Naive_inc_arr & benchable = benchInst.get_benchable();
      run_bench<Naive_inc_arr_bench,Naive_inc_arr>(benchInst, benchable,
                                                   fullname[0], format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }

    // Walk point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + " (" + std::string(filename[0]) + ")";
      Walk_inc_arr_bench benchInst(name, seconds, false);
      Walk_inc_arr & benchable = benchInst.get_benchable();
      run_bench<Walk_inc_arr_bench,Walk_inc_arr>(benchInst, benchable,
                                                 fullname[0], format,
                                                 samples, iterations,
                                                 verbose, postscript);
    }

#if defined(LANDMARK_SUPPORTED)
    // Landmarks point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_LANDMARKS;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + " (" + std::string(filename[0]) + ")";
      Landmarks_inc_arr_bench benchInst(name, seconds, false);
      Landmarks_inc_arr & benchable = benchInst.get_benchable();
      run_bench<Landmarks_inc_arr_bench,Landmarks_inc_arr>(benchInst, benchable,
                                                           fullname[0], format,
                                                           samples, iterations,
                                                           verbose, postscript);
    }
#endif
    
#if defined(TRIANGLE_SUPPORTED)
    // Triangle point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRIANGLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
        NUMBER_TYPE +" " + " (" + std::string(filename[0]) + ")";
      Triangle_inc_arr_bench benchInst(name, seconds, false);
      Triangle_inc_arr & benchable = benchInst.get_benchable();
      run_bench<Triangle_inc_arr_bench,Triangle_inc_arr>(benchInst, benchable,
                                                         fullname[0], format,
                                                         samples, iterations,
                                                         verbose, postscript);
    }
#endif
  }
  
  // Construct Aggregately
  type_id = CGAL::Bench_parse_args::TYPE_AGGREGATE;
  if (type_mask & (0x1 << type_id)) {
    if (verbose) std::cout << "TYPE_AGGREGATE " << std::endl;
    std::string name =
      std::string(parse_args.get_type_name(type_id)) + " " +
      TRAITS_TYPE + " " + KERNEL_TYPE + " " +
      NUMBER_TYPE + " " + " (" + std::string(filename[0]) + ")";
    Agg_arr_bench benchInst(name, seconds, false);
    Aggregate_arr & benchable = benchInst.get_benchable();
    run_bench<Agg_arr_bench,Aggregate_arr>(benchInst, benchable,
                                           fullname[0], format,
                                           samples, iterations,
                                           verbose, postscript);
  }
  
  // Construct and Display
  type_id = CGAL::Bench_parse_args::TYPE_DISPLAY;
  if (type_mask & (0x1 << type_id)) {
#if !defined(USE_CGAL_WINDOW)
    QApplication app(argc, argv);
    App = &app;
#endif
    if (verbose) std::cout << "TYPE_DISPLAY " << std::endl;
    std::string name =
      std::string(parse_args.get_type_name(type_id)) + " " +
      TRAITS_TYPE + " " + KERNEL_TYPE + " " +
      NUMBER_TYPE + " " + " (" + std::string(filename[0]) + ")";
    Dis_arr_bench benchInst(name, seconds, false);
    Display_arr & benchable = benchInst.get_benchable();
    run_bench<Dis_arr_bench,Display_arr>(benchInst, benchable,
                                         fullname[0], format,
                                         samples, iterations,
                                         verbose, postscript);
  }

  return 0;
}
