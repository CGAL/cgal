#include <CGAL/config.h>
#include "short_names.h"
#include <CGAL/basic.h>
#include <CGAL/IO/Color.h>
#include "bench_config.h"
#include "numberType.h"

// PostScript support:
#if BENCH_TRAITS != CONIC_TRAITS && BENCH_TRAITS != CORE_CONIC_TRAITS && \
    BENCH_TRAITS != EXACUS_CONIC_TRAITS && BENCH_TRAITS != CK_CIRCLE_TRAITS && \
    BENCH_TRAITS != CK_CONIC_TRAITS
#define POSTSCRIPT_SUPPORTED 1
#endif

// Kernel:
#if BENCH_KERNEL == LEDA_KERNEL
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>

#elif BENCH_KERNEL == MY_KERNEL
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>

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

#elif BENCH_TRAITS == SEGMENT_CACHED_TRAITS
#include <CGAL/Arr_segment_cached_traits_2.h>

#elif BENCH_TRAITS == LEDA_SEGMENT_TRAITS
#include <CGAL/Arr_leda_segment_traits_2.h>

#elif BENCH_TRAITS == POLYLINE_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#elif BENCH_TRAITS == POLYLINE_CACHED_TRAITS
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_cached_traits_2.h>

#elif BENCH_TRAITS == CONIC_TRAITS
#include <CGAL/Arr_conic_traits_2.h>
#if defined(USE_CGAL_WINDOW)
#include <CGAL/IO/Conic_arc_2_Window_stream.h>
#else
#include <CGAL/IO/Qt_widget_Conic_arc_2.h>
#endif

#elif BENCH_TRAITS == CORE_CONIC_TRAITS
#include <CORE/BigInt.h>
#include <CGAL/Arr_conic_traits_2_core.h>
#if defined(USE_CGAL_WINDOW)
#include <CGAL/IO/Conic_arc_2_Window_stream_core.h>
#else
#include <CGAL/IO/Qt_widget_Conic_arc_2_core.h>
#endif

// Exacus Conics:
#elif BENCH_TRAITS == EXACUS_CONIC_TRAITS
#include <CnX/Conic_sweep_traits_2.h>
#include <CnX/Conic_segment_v2_2.h>
#include <SoX/GAPS/CGAL_Pmwx_2_for_GAPS_traits.h>

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

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>
#if defined(POSTSCRIPT_SUPPORTED)
#include <CGAL/IO/Pm_Postscript_file_stream.h>
#endif

#include <CGAL/Pm_default_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_dummy_point_location.h>

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

// Readers:
// Conic reader:
#if BENCH_TRAITS == SEGMENT_TRAITS
#include "Segment_reader.h"

#elif BENCH_TRAITS == CONIC_TRAITS || BENCH_TRAITS == CORE_CONIC_TRAITS || \
      BENCH_TRAITS == CK_CIRCLE_TRAITS || BENCH_TRAITS == CK_CONIC_TRAITS || \
      BENCH_TRAITS == EXACUS_CONIC_TRAITS
#include "Conic_reader.h"

// Polyline reader:
#elif BENCH_TRAITS == POLYLINE_TRAITS || BENCH_TRAITS == POLYLINE_CACHED_TRAITS
#include "Polyline_reader.h"

#else
#error No traits (TRAITS) specified
#endif

// Kernel:
#if BENCH_KERNEL == LEDA_KERNEL
typedef CGAL::leda_rat_kernel_traits                    Kernel;
#define KERNEL_TYPE "Leda"

#elif BENCH_KERNEL == MY_KERNEL
typedef CGAL::Pm_segment_traits_leda_kernel_2           Kernel;
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

#if defined (USE_INSERT_OLD)
#define INSERT_TYPE "Old"
#else
#define INSERT_TYPE ""
#endif

// Traits:
#if BENCH_TRAITS == SEGMENT_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#define TRAITS_TYPE "Segments"

#elif BENCH_TRAITS == LEDA_SEGMENT_TRAITS
typedef CGAL::Arr_leda_segment_traits_2<Kernel>         Traits;
#define TRAITS_TYPE "Leda Segments"

#elif BENCH_TRAITS == SEGMENT_CACHED_TRAITS
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Traits;
#define TRAITS_TYPE "Cached Segments"

#elif BENCH_TRAITS == POLYLINE_TRAITS
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Traits;
#define TRAITS_TYPE "Polylines"

#elif BENCH_TRAITS == POLYLINE_CACHED_TRAITS
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Segment_traits;
typedef CGAL::Arr_polyline_traits_2<Segment_traits>     Traits;
#define TRAITS_TYPE "Cached Polylines"

#elif BENCH_TRAITS == CONIC_TRAITS
typedef CGAL::Arr_conic_traits_2<Kernel>                Traits;
#define TRAITS_TYPE "Conics"

#elif BENCH_TRAITS == CORE_CONIC_TRAITS
typedef CORE::BigInt                                    CfNT;
typedef CGAL::Arr_conic_traits_2<CfNT,Kernel>           Traits;
#define TRAITS_TYPE "Conics"

// Exacus Conics:
#elif BENCH_TRAITS == EXACUS_CONIC_TRAITS
typedef CnX::Conic_sweep_traits_2<Arithmetic_traits>    CST;
typedef CnX::Conic_segment_2< CST>                      Input_segment;
typedef SoX::CGAL_Pmwx_2_for_GAPS_traits< Input_segment, CST> Traits;
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

// Planar map types:
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;
typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef std::list<Curve_2>                              Curve_list;

// Point location strategies:
typedef CGAL::Pm_default_point_location<Pm>             Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm>               Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm>     Walk_point_location;
typedef CGAL::Pm_dummy_point_location<Pm>               Dummy_point_location;

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

/*! */
inline Window_stream & operator<<(Window_stream & os, Pmwx & pm)
{
  Pmwx::Edge_iterator ei;
  os << CGAL::BLUE;
  for (ei = pm.edges_begin(); ei != pm.edges_end(); ++ei)
    os << (*ei).curve();
  Pmwx::Vertex_iterator vi;
  os << CGAL::RED;
  for (vi = pm.vertices_begin(); vi != pm.vertices_end(); ++vi)
    os << (*vi).point();
  return os;
}

/*! */
class Basic_Pm {
public:
  /*! */
  Basic_Pm() :
    m_filename(0), m_verbose(false), m_postscript(false),
    m_bbox(0.0,0.0,0.0,0.0),
    m_width(1024), m_height(1024)
  { }

  virtual ~Basic_Pm() {}
  
  /*! */
  virtual void op() = 0;

  /*! */
  int init()
  {
#if BENCH_TRAITS == SEGMENT_TRAITS
    Segment_reader<Traits> reader;

#elif BENCH_TRAITS == POLYLINE_TRAITS || BENCH_TRAITS == POLYLINE_CACHED_TRAITS
    Polyline_reader<Traits> reader;

#elif BENCH_TRAITS == CONIC_TRAITS || BENCH_TRAITS == CORE_CONIC_TRAITS || \
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
  void set_file_name(const char * filename) { m_filename = filename; }
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

/*! */
template <class Strategy>
class Increment_pmwx : public Basic_Pm {
public:
  virtual void op()
  {
    Strategy strategy;
    Pmwx pm(&strategy);
    
    Curve_list::const_iterator i;
    for (i = m_curve_list.begin(); i != m_curve_list.end(); i++) pm.insert(*i);
    if (m_verbose) {
      if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
      std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
    }
  }
};

typedef Increment_pmwx<Trap_point_location>     Trap_inc_pmwx;
typedef CGAL::Bench<Trap_inc_pmwx>              Trap_inc_pmwx_bench;

typedef Increment_pmwx<Naive_point_location>    Naive_inc_pmwx;
typedef CGAL::Bench<Naive_inc_pmwx>             Naive_inc_pmwx_bench;

typedef Increment_pmwx<Walk_point_location>     Walk_inc_pmwx;
typedef CGAL::Bench<Walk_inc_pmwx>              Walk_inc_pmwx_bench;

/*! */
template <class Strategy>
class Aggregate_pmwx : public Basic_Pm {
public:
  virtual void op()
  {
    if (m_verbose) {
      std::cout << "Inserting Aggregate" << std::endl;
    }
    Strategy strategy;
    Pmwx pm(&strategy);
#if defined(USE_INSERT_OLD)
    pm.insert_old(m_curve_list.begin(), m_curve_list.end());
#else
    pm.insert(m_curve_list.begin(), m_curve_list.end());
#endif
    if (m_verbose) {
      if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
      std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
    }
  }
};

/*! */
typedef Aggregate_pmwx<Dummy_point_location>    Dummy_agg_pmwx;
typedef CGAL::Bench<Dummy_agg_pmwx>             Dummy_agg_pmwx_bench;

/*! */
template <class Strategy>
class Display_pmwx : public Basic_Pm {
public:
  /*! */
  virtual void op()
  {
    Strategy strategy;
    Pmwx pm(&strategy);
#if 1
#if defined(USE_INSERT_OLD)
    pm.insert_old(m_curve_list.begin(), m_curve_list.end());
#else
    pm.insert(m_curve_list.begin(), m_curve_list.end());
#endif
#else
    Curve_list::const_iterator i;
    for (i = m_curve_list.begin(); i != m_curve_list.end(); i++) pm.insert(*i);
#endif
    if (m_verbose) {
      if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
      std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
    }

#if defined(USE_CGAL_WINDOW)
    m_window->set_flush(0);
    (*m_window) << pm;
    m_window->set_flush(1);
    m_window->flush();
#else
    m_window->lock();
    *m_window << CGAL::BackgroundColor(CGAL::WHITE) << CGAL::RED;
    (*m_window) << pm;
    m_window->unlock();
    App->flush();
#endif

#if defined(POSTSCRIPT_SUPPORTED)
    if (m_postscript) {
      // Print to Postscript file:
      std::cout << "Print to Postscript file" << std::endl;
      CGAL::Postscript_file_stream ps_stream(m_width, m_height ,"pm.ps");
      ps_stream.init(m_x0, m_x1, m_y0);
      // CGAL::cgalize(ps_stream);
      ps_stream.set_line_width(1);
      CGAL::Pm_drawer<Pm, CGAL::Postscript_file_stream> drawer(ps_stream);
      // drawer.draw_faces(pm.faces_begin(), pm.faces_end());
      ps_stream << CGAL::BLUE;
      drawer.draw_halfedges(pm.halfedges_begin(), pm.halfedges_end());
      ps_stream << CGAL::RED;
      drawer.draw_vertices(pm.vertices_begin(), pm.vertices_end());

      // draw_pm(pm, drawer, ps_stream);
      // ps_stream << pm;
    }
#endif
  }
  
  /*! */
  int init()
  {
    int rc = Basic_Pm::init();
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

    m_window->set_redraw(&Display_pmwx::redraw);
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
    Basic_Pm::clean();
    delete m_window;
  }
  
private:
  /*! */
  static void
    redraw(leda_window * wp, double x0, double y0, double x1, double y1) 
  { wp->flush_buffer(x0,y0,x1,y1); }

  Window_stream * m_window;
};

typedef Display_pmwx<Trap_point_location>     Trap_dis_pmwx;
typedef CGAL::Bench<Trap_dis_pmwx>            Trap_dis_pmwx_bench;

typedef Display_pmwx<Naive_point_location>    Naive_dis_pmwx;
typedef CGAL::Bench<Naive_dis_pmwx>           Naive_dis_pmwx_bench;

typedef Display_pmwx<Walk_point_location>     Walk_dis_pmwx;
typedef CGAL::Bench<Walk_dis_pmwx>            Walk_dis_pmwx_bench;

typedef Display_pmwx<Dummy_point_location>    Dummy_dis_pmwx;
typedef CGAL::Bench<Dummy_dis_pmwx>           Dummy_dis_pmwx_bench;

/*! */
template <class Bench_inst, class Benchable>
void run_bench(Bench_inst & benchInst, Benchable & benchable,
               const char * fullname, Format_id format,
               int samples, int iterations, bool verbose,
               bool postscript = false)
{
    // Bench_inst benchInst(name, seconds, false);
    // Benchable & benchable = benchInst.get_benchable();
  benchable.set_format(format);
  benchable.set_file_name(fullname);
  benchable.set_verbose(verbose);
  benchable.set_postscript(postscript);
  
  if (samples > 0) benchInst.set_samples(samples);
  else if (iterations > 0) benchInst.set_iterations(iterations);

  benchInst();

  if (verbose) std::cout << "(" << benchInst.get_samples() << ") "
                         << std::endl;
}

/* */
int main(int argc, char * argv[])
{
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
  const char * filename = parse_args.get_file_name();
  const char * fullname = parse_args.get_full_name();
  bool postscript = parse_args.get_postscript();
  if (postscript) samples = 1;

  if (!filename || !fullname) return -1;

  CGAL::Bench_base::set_name_length(nameLength);
  if (printHeader) CGAL::Bench_base::print_header();

  Type_id type_id;
  Strategy_id strategy_id;
  
  // Construct Incrementaly
  type_id = CGAL::Bench_parse_args::TYPE_INCREMENT;
  if (type_mask & (0x1 << type_id)) {
    // Trapezoidal point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Trap_inc_pmwx_bench benchInst(name, seconds, false);
      Trap_inc_pmwx & benchable = benchInst.get_benchable();
      run_bench<Trap_inc_pmwx_bench,Trap_inc_pmwx>(benchInst, benchable,
                                                   fullname, format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }
    
    // Naive point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Naive_inc_pmwx_bench benchInst(name, seconds, false);
      Naive_inc_pmwx & benchable = benchInst.get_benchable();
      run_bench<Naive_inc_pmwx_bench,Naive_inc_pmwx>(benchInst, benchable,
                                                     fullname, format,
                                                     samples, iterations,
                                                     verbose, postscript);
    }

    // Walk point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Walk_inc_pmwx_bench benchInst(name, seconds, false);
      Walk_inc_pmwx & benchable = benchInst.get_benchable();
      run_bench<Walk_inc_pmwx_bench,Walk_inc_pmwx>(benchInst, benchable,
                                                   fullname, format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }
  }
  
  // Construct Aggregately
  type_id = CGAL::Bench_parse_args::TYPE_AGGREGATE;
  if (type_mask & (0x1 << type_id)) {
    // Dummy point location:
    Strategy_id strategy_id = CGAL::Bench_parse_args::STRATEGY_DUMMY;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Dummy_agg_pmwx_bench benchInst(name, seconds, false);
      Dummy_agg_pmwx & benchable = benchInst.get_benchable();
      run_bench<Dummy_agg_pmwx_bench,Dummy_agg_pmwx>(benchInst, benchable,
                                                     fullname, format,
                                                     samples, iterations,
                                                     verbose, postscript);
    }
  }
  
  // Construct and Display
  type_id = CGAL::Bench_parse_args::TYPE_DISPLAY;
  if (type_mask & (0x1 << type_id)) {
#if !defined(USE_CGAL_WINDOW)
    QApplication app(argc, argv);
    App = &app;
#endif

    // Trapezoidal point location:
    Strategy_id strategy_id = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Trap_dis_pmwx_bench benchInst(name, seconds, false);
      Trap_dis_pmwx & benchable = benchInst.get_benchable();
      run_bench<Trap_dis_pmwx_bench,Trap_dis_pmwx>(benchInst, benchable,
                                                   fullname, format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }

    // Naive point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Naive_dis_pmwx_bench benchInst(name, seconds, false);
      Naive_dis_pmwx & benchable = benchInst.get_benchable();
      run_bench<Naive_dis_pmwx_bench,Naive_dis_pmwx>(benchInst, benchable,
                                                     fullname, format,
                                                     samples, iterations,
                                                     verbose, postscript);
    }

    // Walk point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Walk_dis_pmwx_bench benchInst(name, seconds, false);
      Walk_dis_pmwx & benchable = benchInst.get_benchable();
      run_bench<Walk_dis_pmwx_bench,Walk_dis_pmwx>(benchInst, benchable,
                                                   fullname, format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }

    // Dummy point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_DUMMY;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Dummy_dis_pmwx_bench benchInst(name, seconds, false);
      Dummy_dis_pmwx & benchable = benchInst.get_benchable();
      run_bench<Dummy_dis_pmwx_bench,Dummy_dis_pmwx>(benchInst, benchable,
                                                     fullname, format,
                                                     samples, iterations,
                                                     verbose, postscript);
    }
  }
  
  return 0;
}
