#include <CGAL/config.h>
#include "short_names.h"
#include <CGAL/basic.h>
#include <CGAL/IO/Color.h>
#include "bench_config.h"
#include "numberType.h"

enum MaxFilesNumber {
  MAX_FILES_NUMBER = 20
};

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

#if BENCH_PM == PLANAR_MAP
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_default_dcel.h>
#elif BENCH_PM == PLANAR_MAP_WITH_INTERSECTIONS
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/Pm_default_dcel.h>
#elif BENCH_PM == ARRANGEMENT
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_2_default_dcel.h>
#endif

#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>
#if defined(POSTSCRIPT_SUPPORTED)
#include <CGAL/IO/Pm_Postscript_file_stream.h>
#endif

#include <CGAL/Pm_trapezoid_dag_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_dummy_point_location.h>
#include <CGAL/Pm_triangle_point_location.h>
#include <CGAL/Pm_simple_point_location.h>

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
#if BENCH_TRAITS == SEGMENT_TRAITS || BENCH_TRAITS == SEGMENT_CACHED_TRAITS
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
#if BENCH_PM == PLANAR_MAP
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Planar_map;
typedef Planar_map                                      Pm;
#define PLANAR_MAP_TYPE "Pm"
#elif BENCH_PM == PLANAR_MAP_WITH_INTERSECTIONS
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Planar_map;
#define PLANAR_MAP_TYPE "Pmwx"
#elif BENCH_PM == ARRANGEMENT
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Planar_map;
typedef Planar_map::Planar_map                          Pm;
#define PLANAR_MAP_TYPE "Arr"
#endif

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef std::list<Curve_2>                              Curve_list;
typedef std::list<Point_2>                              Point_list;
typedef Planar_map::Locate_type                         Locate_type;
typedef Planar_map::Halfedge_iterator                   Pm_Halfedge_iterator;

// Point location strategies:
typedef CGAL::Pm_trapezoid_dag_point_location<Pm>       Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm>               Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm>     Walk_point_location;
typedef CGAL::Pm_dummy_point_location<Pm>               Dummy_point_location;
typedef CGAL::Pm_triangle_point_location<Pm>            Triangle_point_location;
typedef CGAL::Pm_simple_point_location<Pm>              Simple_point_location;

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

inline std::ostream & operator<<(std::ostream & os,
                                 const Planar_map::Vertex & vertex)
{
  os << vertex.point();
  return os;
}

CGAL_END_NAMESPACE

/*! */
inline Window_stream & operator<<(Window_stream & os, Planar_map & pm)
{
  Planar_map::Edge_iterator ei;
  os << CGAL::BLUE;
  for (ei = pm.edges_begin(); ei != pm.edges_end(); ++ei)
    os << (*ei).curve();
  Planar_map::Vertex_iterator vi;
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
  //virtual void op() = 0;

  /*! */
  int init()
  {
#if BENCH_TRAITS == SEGMENT_TRAITS || BENCH_TRAITS == SEGMENT_CACHED_TRAITS
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

/*! 
This class is supplied to the bench as a templatwe parameter.
it defines the operations we want to measure.
in this case we measure only the incremental construction 
of the planar map.
*/
template <class Strategy>
class Increment_pm : public Basic_Pm {
public:
  // op is the main function that must be implemented in each benchable class
  // that is supplied to the Bench class. 
  // the other functions that should be implemented are inherited from Basic_Pm
  // This (op) is the main method, that is different from all other classes 
  // inherite from Basic_Pm
  void op()
  {
    Strategy strategy;
    Planar_map pm(&strategy);
    
    Curve_list::const_iterator i;
    //iterator on all the curves got from the input, 
    //insert each one - incrementally - to the pm.
    for (i = m_curve_list.begin(); i != m_curve_list.end(); i++) 
    {
      pm.insert(*i);
    }
    
    if (m_verbose) {      //print to cout
      if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
      std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
    }

#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
    // Clear the exacus cache:
    CST::Algebraic_curve_2::clear_cache();
    CST::Point_2::clear_cache();
#endif
  }
};

typedef Increment_pm<Trap_point_location>       Trap_inc_pm;
typedef CGAL::Bench<Trap_inc_pm>                Trap_inc_pm_bench;

typedef Increment_pm<Naive_point_location>      Naive_inc_pm;
typedef CGAL::Bench<Naive_inc_pm>               Naive_inc_pm_bench;

typedef Increment_pm<Walk_point_location>       Walk_inc_pm;
typedef CGAL::Bench<Walk_inc_pm>                Walk_inc_pm_bench;

typedef Increment_pm<Triangle_point_location>   Triangle_inc_pm;
typedef CGAL::Bench<Triangle_inc_pm>            Triangle_inc_pm_bench;

typedef Increment_pm<Simple_point_location>     Simple_inc_pm;
typedef CGAL::Bench<Simple_inc_pm>              Simple_inc_pm_bench;

/*! 
*/
template <class Strategy>
class Aggregate_pm : public Basic_Pm {
public:
  void op()
  {
    if (m_verbose) {
      std::cout << "Inserting Aggregate" << std::endl;
    }
    Strategy strategy;
    Planar_map pm(&strategy);
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
#if BENCH_TRAITS == EXACUS_CONIC_TRAITS
    // Clear the exacus cache:
    CST::Algebraic_curve_2::clear_cache();
    CST::Point_2::clear_cache();
#endif
  }
};

/*! */
typedef Aggregate_pm<Dummy_point_location>      Dummy_agg_pm;
typedef CGAL::Bench<Dummy_agg_pm>               Dummy_agg_pm_bench;

typedef Aggregate_pm<Naive_point_location>      Naive_agg_pm;
typedef CGAL::Bench<Naive_agg_pm>               Naive_agg_pm_bench;

typedef Aggregate_pm<Simple_point_location>     Simple_agg_pm;
typedef CGAL::Bench<Simple_agg_pm>              Simple_agg_pm_bench;

typedef Aggregate_pm<Triangle_point_location>   Triangle_agg_pm;
typedef CGAL::Bench<Triangle_agg_pm>            Triangle_agg_pm_bench;

typedef Aggregate_pm<Trap_point_location>       Trap_agg_pm;
typedef CGAL::Bench<Trap_agg_pm>                Trap_agg_pm_bench;

typedef Aggregate_pm<Walk_point_location>       Walk_agg_pm;
typedef CGAL::Bench<Walk_agg_pm>                Walk_agg_pm_bench;


/*! */
template <class Strategy>
class Display_pm : public Basic_Pm {
public:
  /*! */
  void op()
  {
    Strategy strategy;
    Planar_map pm(&strategy);
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

      std::copy(pm.vertices_begin(), pm.vertices_end(),
                std::ostream_iterator<Planar_map::Vertex>(std::cout, "\n"));
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
      CGAL::Pm_drawer<Planar_map, CGAL::Postscript_file_stream> drawer(ps_stream);
      // drawer.draw_faces(pm.faces_begin(), pm.faces_end());
      ps_stream << CGAL::BLUE;
      drawer.draw_halfedges(pm.halfedges_begin(), pm.halfedges_end());
      ps_stream << CGAL::RED;
      drawer.draw_vertices(pm.vertices_begin(), pm.vertices_end());

      // draw_pm(pm, drawer, ps_stream);
      // ps_stream << pm;
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

    m_window->set_redraw(&Display_pm::redraw);
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

typedef Display_pm<Trap_point_location>         Trap_dis_pm;
typedef CGAL::Bench<Trap_dis_pm>                Trap_dis_pm_bench;

typedef Display_pm<Naive_point_location>        Naive_dis_pm;
typedef CGAL::Bench<Naive_dis_pm>               Naive_dis_pm_bench;

typedef Display_pm<Walk_point_location>         Walk_dis_pm;
typedef CGAL::Bench<Walk_dis_pm>                Walk_dis_pm_bench;

typedef Display_pm<Dummy_point_location>        Dummy_dis_pm;
typedef CGAL::Bench<Dummy_dis_pm>               Dummy_dis_pm_bench;

typedef Display_pm<Triangle_point_location>     Triangle_dis_pm;
typedef CGAL::Bench<Triangle_dis_pm>            Triangle_dis_pm_bench;

typedef Display_pm<Simple_point_location>       Simple_dis_pm;
typedef CGAL::Bench<Simple_dis_pm>              Simple_dis_pm_bench;

/*************** ixx ********************/
/*! */
template <class Strategy>
class Locate_Pm {
public:
  //  enum MaxFilesNumber {
  //MAX_FILES_NUMBER = 20
  //};

  /*! */
  Locate_Pm() :
    m_verbose(false), m_postscript(false), 
    m_strategy(), m_pm(&m_strategy)
    //m_bbox(0.0,0.0,0.0,0.0),
    //m_width(1024), m_height(1024)
  { 
    for (int j=0; j<MAX_FILES_NUMBER; j++) {
      m_filename[j] = 0;
    }
  }

  virtual ~Locate_Pm() {}
  
  /*! */
  virtual void op() 
  {
    if (m_verbose) std::cout << "op Locate_Pm" << std::endl;
    //    Strategy strategy;
    Locate_type lt;
    Pm_Halfedge_iterator e;
    Point_list::const_iterator i;

    //iterator on all the points got from the input, 
    //go over each one and locate it in the pm .
    for (i = m_point_list.begin(); i != m_point_list.end(); i++) 
    {
      e = m_pm.locate(*i,lt);
      if (m_verbose) {
	std::cout << "locate point"<< (*i) <<" at " ;
	//print output
	if (lt==Planar_map::UNBOUNDED_FACE) std::cout << "Unbounded face" << std::endl;
	else if (lt==Planar_map::FACE) std::cout << "Face that is left of " 
	 << e->source()->point()  <<" towards "<< e->target()->point() << std::endl;
	else if (lt==Planar_map::EDGE) std::cout << "Edge :" 
	 << e->source()->point()  <<" towards "<< e->target()->point() << std::endl;
	else if (lt==Planar_map::VERTEX) std::cout << "vertex: "
				  << e->target()->point() << std::endl;
	else std::cout << "Unknown locate type" << std::endl;
      }
    }
  }

  /*! */
  int init()
  {
    if (m_verbose) {
      std::cout << "init Locate_Pm " << std::endl;
      std::cout << "file[0] = "<< m_filename[0] 
		<< "file[1] = " << m_filename[1]  << std::endl;
    }
    //read the planar map from the first input file
    std::ifstream inp(m_filename[0]);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << m_filename[0] << "!" << std::endl;
      return -1;
    }
    //inp >> pm;
    m_pm.read(inp);

    //ixx: create points reader ?
    Point_reader<Traits> reader;
    int rc = reader.read_data(m_filename[1], std::back_inserter(m_point_list),
                              m_format);
    if (rc < 0) return rc;
    if (m_verbose) std::cout << m_point_list.size() << " points" << std::endl;

    return 0;
  }
    
  /* */
  void clean() 
  {
    if (m_verbose) std::cout << "clean Locate_Pm" << std::endl;
    m_point_list.clear(); 
    m_pm.clear();
  }
  void sync(){}

  void set_format(Format_id fmt) { m_format = fmt; }
  
  void set_file_name(const char * filename, int file_index = 0) 
  { m_filename[file_index] = filename; }

  void set_verbose(const bool verbose) { m_verbose = verbose; }
  void set_postscript(const bool postscript) { m_postscript = postscript; }

protected:
  const char * m_filename[MAX_FILES_NUMBER];
  Point_list m_point_list;
  bool m_verbose;
  bool m_postscript;
  Format_id m_format;
  //  CGAL::Bbox_2 m_bbox;
  //  Pmwx m_pm;
  Strategy m_strategy;
  //  Pmwx pm(&strategy);
  Planar_map m_pm;

  //  int m_width, m_height;
  //  float m_x0, m_x1, m_y0, m_y1;
};

typedef Locate_Pm<Trap_point_location>     Trap_loc_pm;
typedef CGAL::Bench<Trap_loc_pm>           Trap_loc_pm_bench;

typedef Locate_Pm<Naive_point_location>    Naive_loc_pm;
typedef CGAL::Bench<Naive_loc_pm>          Naive_loc_pm_bench;

typedef Locate_Pm<Walk_point_location>     Walk_loc_pm;
typedef CGAL::Bench<Walk_loc_pm>           Walk_loc_pm_bench;

typedef Locate_Pm<Triangle_point_location> Triangle_loc_pm;
typedef CGAL::Bench<Triangle_loc_pm>       Triangle_loc_pm_bench;

typedef Locate_Pm<Simple_point_location>   Simple_loc_pm;
typedef CGAL::Bench<Simple_loc_pm>         Simple_loc_pm_bench;

/***************************************/


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
  //  std::cout << "main " << std::endl;
  //get arguments
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
  int files_number = parse_args.get_files_number();
  //  std::cout << "files_number: " << files_number << std::endl;
  const char * filename[MAX_FILES_NUMBER];
  const char * fullname[MAX_FILES_NUMBER];
  for (int i=0; i< files_number; i++)
  {
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
    // Trapezoidal point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_id)) {
      //set the name of this run
      //get type name (increment or aggregate or display)
      //get strategy (trapezoidal / naive ...)
      //get traits type, kernel type, number type, input filename, etc ...
      std::string name =
	std::string(parse_args.get_type_name(type_id)) + " " +
	std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      //create an instance of the Bench class that will be performed
      Trap_inc_pm_bench benchInst(name, seconds, false);
      //create an instance of the benchable class that is passed as a template
      //parameter to the Bench class, and that its op function is performed
      Trap_inc_pm & benchable = benchInst.get_benchable();
      //run the bench with the template classes passed as above
      //and other parameters taken from the input
      run_bench<Trap_inc_pm_bench,Trap_inc_pm>(benchInst, benchable,
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
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Naive_inc_pm_bench benchInst(name, seconds, false);
      Naive_inc_pm & benchable = benchInst.get_benchable();
      run_bench<Naive_inc_pm_bench,Naive_inc_pm>(benchInst, benchable,
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
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Walk_inc_pm_bench benchInst(name, seconds, false);
      Walk_inc_pm & benchable = benchInst.get_benchable();
      run_bench<Walk_inc_pm_bench,Walk_inc_pm>(benchInst, benchable,
                                               fullname[0], format,
                                               samples, iterations,
                                               verbose, postscript);
    }
     
    // Triangle point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRIANGLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
	std::string(parse_args.get_type_name(type_id)) + " " +
	std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	NUMBER_TYPE +" " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Triangle_inc_pm_bench benchInst(name, seconds, false);
      Triangle_inc_pm & benchable = benchInst.get_benchable();
      run_bench<Triangle_inc_pm_bench,Triangle_inc_pm>(benchInst, benchable,
						       fullname[0], format,
						       samples, iterations,
						       verbose, postscript);
    }

    // Simple point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_SIMPLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
	std::string(parse_args.get_type_name(type_id)) + " " +
	std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	NUMBER_TYPE +" " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Simple_inc_pm_bench benchInst(name, seconds, false);
      Simple_inc_pm & benchable = benchInst.get_benchable();
      run_bench<Simple_inc_pm_bench,Simple_inc_pm>(benchInst, benchable,
						   fullname[0], format,
						   samples, iterations,
						   verbose, postscript);
    }
  }
  
  // Construct Aggregately
  type_id = CGAL::Bench_parse_args::TYPE_AGGREGATE;
  if (type_mask & (0x1 << type_id)) {
    if (verbose) std::cout << "TYPE_AGGREGATE " << std::endl;
    // Dummy point location:
    Strategy_id strategy_id = CGAL::Bench_parse_args::STRATEGY_DUMMY;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Dummy_agg_pm_bench benchInst(name, seconds, false);
      Dummy_agg_pm & benchable = benchInst.get_benchable();
      run_bench<Dummy_agg_pm_bench,Dummy_agg_pm>(benchInst, benchable,
                                                 fullname[0], format,
                                                 samples, iterations,
                                                 verbose, postscript);
    }
    // Trapezoidal point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_id)) {
      //set the name of this run
      //get type name (increment or aggregate or display)
      //get strategy (trapezoidal / naive ...)
      //get traits type, kernel type, number type, input filename, etc ...
      std::string name =
	std::string(parse_args.get_type_name(type_id)) + " " +
	std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      //create an instance of the Bench class that will be performed
      Trap_agg_pm_bench benchInst(name, seconds, false);
      //create an instance of the benchable class that is passed as a template
      //parameter to the Bench class, and that its op function is performed
      Trap_agg_pm & benchable = benchInst.get_benchable();
      //run the bench with the template classes passed as above
      //and other parameters taken from the input
      run_bench<Trap_agg_pm_bench,Trap_agg_pm>(benchInst, benchable,
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
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Naive_agg_pm_bench benchInst(name, seconds, false);
      Naive_agg_pm & benchable = benchInst.get_benchable();
      run_bench<Naive_agg_pm_bench,Naive_agg_pm>(benchInst, benchable,
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
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Walk_agg_pm_bench benchInst(name, seconds, false);
      Walk_agg_pm & benchable = benchInst.get_benchable();
      run_bench<Walk_agg_pm_bench,Walk_agg_pm>(benchInst, benchable,
                                               fullname[0], format,
                                               samples, iterations,
                                               verbose, postscript);
    }
     
    // Triangle point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRIANGLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
	std::string(parse_args.get_type_name(type_id)) + " " +
	std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	NUMBER_TYPE +" " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Triangle_agg_pm_bench benchInst(name, seconds, false);
      Triangle_agg_pm & benchable = benchInst.get_benchable();
      run_bench<Triangle_agg_pm_bench,Triangle_agg_pm>(benchInst, benchable,
						       fullname[0], format,
						       samples, iterations,
						       verbose, postscript);
    }

    // Simple point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_SIMPLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
	std::string(parse_args.get_type_name(type_id)) + " " +
	std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	NUMBER_TYPE +" " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Simple_agg_pm_bench benchInst(name, seconds, false);
      Simple_agg_pm & benchable = benchInst.get_benchable();
      run_bench<Simple_agg_pm_bench,Simple_agg_pm>(benchInst, benchable,
						   fullname[0], format,
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
    if (verbose) std::cout << "TYPE_DISPLAY " << std::endl;
    // Trapezoidal point location:
    Strategy_id strategy_id = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Trap_dis_pm_bench benchInst(name, seconds, false);
      Trap_dis_pm & benchable = benchInst.get_benchable();
      run_bench<Trap_dis_pm_bench,Trap_dis_pm>(benchInst, benchable,
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
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Naive_dis_pm_bench benchInst(name, seconds, false);
      Naive_dis_pm & benchable = benchInst.get_benchable();
      run_bench<Naive_dis_pm_bench,Naive_dis_pm>(benchInst, benchable,
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
          PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	  NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Walk_dis_pm_bench benchInst(name, seconds, false);
      Walk_dis_pm & benchable = benchInst.get_benchable();
      run_bench<Walk_dis_pm_bench,Walk_dis_pm>(benchInst, benchable,
                                                   fullname[0], format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }

    // Triangle point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRIANGLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	  NUMBER_TYPE +" " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Triangle_dis_pm_bench benchInst(name, seconds, false);
      Triangle_dis_pm & benchable = benchInst.get_benchable();
      run_bench<Triangle_dis_pm_bench,Triangle_dis_pm>(benchInst, benchable,
                                                   fullname[0], format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }

    // Simple point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_SIMPLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          PLANAR_MAP_TYPE +" " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	  NUMBER_TYPE +" " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Simple_dis_pm_bench benchInst(name, seconds, false);
      Simple_dis_pm & benchable = benchInst.get_benchable();
      run_bench<Simple_dis_pm_bench,Simple_dis_pm>(benchInst, benchable,
                                                   fullname[0], format,
                                                   samples, iterations,
                                                   verbose, postscript);
    }

    // Dummy point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_DUMMY;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
        std::string(parse_args.get_type_name(type_id)) + " " +
        std::string(parse_args.get_strategy_name(strategy_id)) + " " +
        PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " +
        NUMBER_TYPE + " " + INSERT_TYPE + 
	" (" + std::string(filename[0]) + ")";
      Dummy_dis_pm_bench benchInst(name, seconds, false);
      Dummy_dis_pm & benchable = benchInst.get_benchable();
      run_bench<Dummy_dis_pm_bench,Dummy_dis_pm>(benchInst, benchable,
                                                 fullname[0], format,
                                                 samples, iterations,
                                                 verbose, postscript);
    }

  }
  
  /************************ ixx **************************************/

  // Locate - using Point Location's strategies 
  //(only if type_id == point_location)
  type_id = CGAL::Bench_parse_args::TYPE_POINT_LOCATION;
  if (type_mask & (0x1 << type_id)) {
    if (verbose) std::cout << "TYPE_POINT_LOCATION" << std::endl;
    // Trapezoidal point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategy_mask & (0x1 << strategy_id)) {
      //set the name of this run
      //get type name (increment or aggregate or display)
      //get strategy (trapezoidal / naive ...)
      //get traits type, kernel type, number type, input filename, etc ...
      std::string name =
	  std::string(parse_args.get_type_name(type_id)) + " " +
	  std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	  PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	  NUMBER_TYPE + " " + INSERT_TYPE + 
	  " (" + std::string(filename[0]) + ")" +
	  " (" + std::string(filename[1]) + ")";
      //create an instance of the Bench class that will be performed
      Trap_loc_pm_bench benchInst(name, seconds, false);
      //create an instance of the benchable class that is passed as a template
      //parameter to the Bench class, and that its op function is performed
      Trap_loc_pm & benchable = benchInst.get_benchable();
      //run the bench with the template classes passed as above
      //and other parameters taken from the input
      run_bench<Trap_loc_pm_bench,Trap_loc_pm>(benchInst, benchable,
					       fullname[0], format,
					       samples, iterations,
					       verbose, postscript,
					       fullname[1]);
    }
    
    // Naive point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
          std::string(parse_args.get_type_name(type_id)) + " " +
          std::string(parse_args.get_strategy_name(strategy_id)) + " " +
          PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	  NUMBER_TYPE + " " + INSERT_TYPE + 
	  " (" + std::string(filename[0]) + ")" +
	  " (" + std::string(filename[1]) + ")";
      Naive_loc_pm_bench benchInst(name, seconds, false);
      Naive_loc_pm & benchable = benchInst.get_benchable();
      run_bench<Naive_loc_pm_bench,Naive_loc_pm>(benchInst, benchable,
						 fullname[0], format,
						 samples, iterations,
						 verbose, postscript,
						 fullname[1]);
    }

    // Walk point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
	  std::string(parse_args.get_type_name(type_id)) + " " +
	  std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	  PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
          NUMBER_TYPE + " " + INSERT_TYPE + 
          " (" + std::string(filename[0]) + ")" +
	  " (" + std::string(filename[1]) + ")";
      Walk_loc_pm_bench benchInst(name, seconds, false);
      Walk_loc_pm & benchable = benchInst.get_benchable();
      run_bench<Walk_loc_pm_bench,Walk_loc_pm>(benchInst, benchable,
					       fullname[0], format,
					       samples, iterations,
					       verbose, postscript,
					       fullname[1]);
    }
     
    // Triangle point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_TRIANGLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
	  std::string(parse_args.get_type_name(type_id)) + " " +
	  std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	  PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
	  NUMBER_TYPE + " " + INSERT_TYPE + 
          " (" + std::string(filename[0]) + ")" +
	  " (" + std::string(filename[1]) + ")";
      Triangle_loc_pm_bench benchInst(name, seconds, false);
      Triangle_loc_pm & benchable = benchInst.get_benchable();
      run_bench<Triangle_loc_pm_bench,Triangle_loc_pm>(benchInst, benchable,
						       fullname[0], format,
						       samples, iterations,
						       verbose, postscript,
						       fullname[1]);
    }

    // Simple point location:
    strategy_id = CGAL::Bench_parse_args::STRATEGY_SIMPLE;
    if (strategy_mask & (0x1 << strategy_id)) {
      std::string name =
	  std::string(parse_args.get_type_name(type_id)) + " " +
	  std::string(parse_args.get_strategy_name(strategy_id)) + " " +
	  PLANAR_MAP_TYPE + " " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + 
          NUMBER_TYPE + " " + INSERT_TYPE + 
          " (" + std::string(filename[0]) + ")" +
	  " (" + std::string(filename[1]) + ")";
      Simple_loc_pm_bench benchInst(name, seconds, false);
      Simple_loc_pm & benchable = benchInst.get_benchable();
      run_bench<Simple_loc_pm_bench,Simple_loc_pm>(benchInst, benchable,
						   fullname[0], format,
						   samples, iterations,
						   verbose, postscript,
						   fullname[1]);
    }
  }

  return 0;
}
