#include <CGAL/config.h>

#include "short_names.h"

#include <CGAL/basic.h>

// Kernel:
#if defined(USE_LEDA_KERNEL)
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#elif defined(USE_MY_KERNEL)
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>
#if !defined(USE_LEDA_SEGMENT_TRAITS)
#error Must define USE_LEDA_SEGMENT_TRAITS!
#endif
#elif defined(USE_SIMPLE_CARTESIAN_KERNEL)
#include <CGAL/Simple_cartesian.h>
#else
#include <CGAL/Cartesian.h>
#endif

// Traits:
#if defined(USE_CONIC_TRAITS)
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/IO/Conic_arc_2_Window_stream.h>
#elif defined(USE_LEDA_SEGMENT_TRAITS)
#include <CGAL/Arr_leda_segment_traits_2.h>
#elif defined(USE_SEGMENT_CACHED_TRAITS)
#include <CGAL/Arr_segment_cached_traits_2.h>
#elif defined(USE_POLYLINE_TRAITS)
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#else
#include <CGAL/Arr_segment_traits_2.h>
#endif

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>

#include <CGAL/Pm_default_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_dummy_point_location.h>

#include <CGAL/Bench.h>
#include <CGAL/Bench.C>
#include <CGAL/Bench_parse_args.h>
#include <CGAL/Bench_parse_args.C>

#include <stdlib.h>
#include <iostream>
#include <list>

#include "numberType.h"

#if defined(USE_CONIC_TRAITS)
#include "Conic_reader.h"
#elif defined(USE_POLYLINE_TRAITS)
#include "Polyline_reader.h"
#else
#include "Segment_reader.h"
#endif

#if defined(USE_LEDA_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Kernel;
#define KERNEL_TYPE "Leda"
#elif defined(USE_MY_KERNEL)
typedef CGAL::Pm_segment_traits_leda_kernel_2           Kernel;
#define KERNEL_TYPE "Partial Leda"
#elif defined(USE_SIMPLE_CARTESIAN_KERNEL)
typedef CGAL::Simple_cartesian<WNT>                     Kernel;
#define KERNEL_TYPE "Simple Cartesian"
#else
typedef CGAL::Cartesian<WNT>                            Kernel;
#define KERNEL_TYPE "Cartesian"
#endif

#if defined (USE_INSERT_OLD)
#define INSERT_TYPE "Old"
#else
#define INSERT_TYPE ""
#endif

#if defined(USE_CONIC_TRAITS)
typedef CGAL::Arr_conic_traits_2<Kernel>                Traits;
#define TRAITS_TYPE "Conics"
#elif defined(USE_LEDA_SEGMENT_TRAITS)
typedef CGAL::Arr_leda_segment_traits_2<Kernel>         Traits;
#define TRAITS_TYPE "Leda Segments"
#elif defined(USE_SEGMENT_CACHED_TRAITS)
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Traits;
#define TRAITS_TYPE "Cached Segments"
#elif defined(USE_POLYLINE_TRAITS)
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       SegmentTraits;
typedef CGAL::Arr_polyline_traits_2<SegmentTraits>      Traits;
#define TRAITS_TYPE "Polylines"
#else
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#define TRAITS_TYPE "Segments"
#endif

typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;
typedef Traits::Point_2                                 Point;
typedef Traits::X_monotone_curve_2                      Curve;
typedef std::list<Curve>                                CurveList;

typedef CGAL::Pm_default_point_location<Pm>             Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm>               Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm>     Walk_point_location;
typedef CGAL::Pm_dummy_point_location<Pm>               Dummy_point_location;

typedef CGAL::Bench_parse_args::TypeId                  TypeId;
typedef CGAL::Bench_parse_args::StrategyId              StrategyId;
typedef CGAL::Bench_parse_args::FormatId                FormatId;

#if defined(USE_LEDA_KERNEL) || defined(USE_MY_KERNEL)
inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os,
                                        const Point & p)
{
  os << leda_point(p.xcoordD(), p.ycoordD());
  return os;
}

#if defined(USE_SEGMENT_CACHED_TRAITS)
inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os,
                                        const Curve & curve)
{
  Kernel::Segment_2 seg = static_cast<Kernel::Segment_2>(curve);
  os << leda_segment(seg.xcoord1D(), seg.ycoord1D(),
                     seg.xcoord2D(), seg.ycoord2D());
  return os;
}
#else
inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os,
                                        const Kernel::Segment_2 & seg)
{
  os << leda_segment(seg.xcoord1D(), seg.ycoord1D(),
                     seg.xcoord2D(), seg.ycoord2D());
  return os;
}
#endif

#else

#if defined(USE_SEGMENT_CACHED_TRAITS)
inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os,
                                        const Curve & curve)
{
  os << static_cast<Kernel::Segment_2>(curve);  
  return os;
}
#endif

#endif

#if defined(USE_POLYLINE_TRAITS)
inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os, 
                                        const Curve & cv)
{
  Curve::const_iterator ps = cv.begin();
  Curve::const_iterator pt = ps; pt++;

  while (pt != cv.end())
  {
    os << Kernel::Segment_2(*ps, *pt);
    ps++; pt++;
  }
  return (os);
}
#endif

inline CGAL::Window_stream & operator<<(CGAL::Window_stream & os, Pmwx & pm)
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

/*
 */
class Basic_Pm {
public:
  /*
   */
  Basic_Pm() : m_filename(0), m_verbose(false), m_bbox(0.0,0.0,0.0,0.0) { }

  virtual ~Basic_Pm() {}
  
  /*
   */
  virtual void op() = 0;

  /*
   */
  int init()
  {
#if defined(USE_CONIC_TRAITS)
    Conic_reader<Traits> reader;
#elif defined(USE_POLYLINE_TRAITS)
    Polyline_reader<Traits> reader;
#else
    Segment_reader<Traits> reader;
#endif
    int rc = reader.ReadData(m_filename, m_curveList, m_format, m_bbox);
    if (rc < 0) return rc;
    if (m_verbose) std::cout << m_curveList.size() << " curves" << std::endl;

    return 0;
  }
    
  /*
   */
  void clean() { m_curveList.clear(); }
  void sync(){}

  void set_format(CGAL::Bench_parse_args::FormatId fmt) { m_format = fmt; }
  void set_file_name(const char * filename) { m_filename = filename; }
  void set_verbose(const bool verbose) { m_verbose = verbose; }

protected:
  const char * m_filename;
  CurveList m_curveList;
  bool m_verbose;
  CGAL::Bench_parse_args::FormatId m_format;
  CGAL::Bbox_2 m_bbox;
};

/*!
 */
template <class Strategy>
class Increment_pmwx : public Basic_Pm {
public:
  virtual void op()
  {
    Strategy strategy;
    Pmwx pm(&strategy);
    CurveList::const_iterator i;
    for (i = m_curveList.begin(); i != m_curveList.end(); i++)
        pm.insert(*i);
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

/*!
 */
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
    pm.insert_old(m_curveList.begin(), m_curveList.end());
#else
    pm.insert(m_curveList.begin(), m_curveList.end());
#endif
    if (m_verbose) {
      if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
      std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
    }
  }
};

typedef Aggregate_pmwx<Dummy_point_location>    Dummy_agg_pmwx;
typedef CGAL::Bench<Dummy_agg_pmwx>             Dummy_agg_pmwx_bench;

/*!
 */
template <class Strategy>
class Display_pmwx : public Basic_Pm {
private:
  typedef CGAL::Window_stream Window_stream;

public:
  /*!
   */
  virtual void op()
  {
    Strategy strategy;
    Pmwx pm(&strategy);
#if defined(USE_INSERT_OLD)
    pm.insert_old(m_curveList.begin(), m_curveList.end());
#else
    pm.insert(m_curveList.begin(), m_curveList.end());
#endif
    if (m_verbose) {
      if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
      std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
      std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
      std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
    }

    m_window->set_flush(0);
    (*m_window) << pm;
    m_window->set_flush(1);
    m_window->flush();
  }
  
  /*!
   */
  int init()
  {
    int rc = Basic_Pm::init();
    if (rc < 0) return rc;

    float x_range = m_bbox.xmax() - m_bbox.xmin();
    float y_range = m_bbox.ymax() - m_bbox.ymin();
    float width = 640;
    float height = (y_range * width) / x_range;
    
    m_window = new Window_stream(static_cast<int>(width),
                                 static_cast<int>(height));
    if (!m_window) return -1;

    float min_range = (x_range < y_range) ? x_range : y_range;
    float x_margin = min_range / 4;
    float y_margin = (height * x_margin) / width;
        
    float x0 = m_bbox.xmin() - x_margin;
    float x1 = m_bbox.xmax() + x_margin;
    float y0 = m_bbox.ymin() - y_margin;
    m_window->init(x0, x1, y0);   // logical window size 

    m_window->set_redraw(&Display_pmwx::redraw);
    m_window->set_mode(leda_src_mode);
    m_window->set_node_width(3);
    m_window->set_point_style(leda_cross_point);
    m_window->set_line_width(1);
    m_window->display(leda_window::center, leda_window::center);
    return 0;
  }

  /*!
   */
  void clean()
  {
    Basic_Pm::clean();
    delete m_window;
  }
  
private:
  /*!
   */
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

/*
 */
template <class Bench_inst, class Benchable>
void runBench(Bench_inst & benchInst, Benchable & benchable,
              const char * fullname, FormatId format,
              int samples, int iterations, bool verbose)
{
    // Bench_inst benchInst(name, seconds, false);
    // Benchable & benchable = benchInst.get_benchable();
  benchable.set_format(format);
  benchable.set_file_name(fullname);
  benchable.set_verbose(verbose);

  if (samples > 0) benchInst.set_samples(samples);
  else if (iterations > 0) benchInst.set_iterations(iterations);

  benchInst();

  if (verbose) std::cout << "(" << benchInst.get_samples() << ") "
                         << std::endl;
}

/*
 */
int main(int argc, char * argv[])
{
  CGAL::Bench_parse_args parseArgs(argc, argv);
  int rc = parseArgs.parse();
  if (rc > 0) return 0;
  if (rc < 0) return rc;
  
  bool verbose = parseArgs.get_verbose();
  unsigned int typeMask = parseArgs.get_type_mask();
  unsigned int strategyMask = parseArgs.get_strategy_mask();
  FormatId format = parseArgs.get_input_format();
  int samples = parseArgs.get_samples();
  int iterations = parseArgs.get_iterations();
  int seconds = parseArgs.get_seconds();
  bool printHeader = parseArgs.get_print_header();
  int nameLength = parseArgs.get_name_length();
  const char * filename = parseArgs.get_file_name();
  const char * fullname = parseArgs.get_full_name();

  if (!filename || !fullname) return -1;
      
  CGAL::Bench_base::set_name_length(nameLength);
  if (printHeader) CGAL::Bench_base::print_header();
  
  // Construct Incrementaly
  TypeId typeId = CGAL::Bench_parse_args::TYPE_INCREMENT;
  if (typeMask & (0x1 << typeId)) {
    // Trapezoidal point location:
    StrategyId strategyId = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Trap_inc_pmwx_bench benchInst(name, seconds, false);
      Trap_inc_pmwx & benchable = benchInst.get_benchable();
      runBench<Trap_inc_pmwx_bench,Trap_inc_pmwx>(benchInst, benchable,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }
    
    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Naive_inc_pmwx_bench benchInst(name, seconds, false);
      Naive_inc_pmwx & benchable = benchInst.get_benchable();
      runBench<Naive_inc_pmwx_bench,Naive_inc_pmwx>(benchInst, benchable,
                                                    fullname, format,
                                                    samples, iterations,
                                                    verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Walk_inc_pmwx_bench benchInst(name, seconds, false);
      Walk_inc_pmwx & benchable = benchInst.get_benchable();
      runBench<Walk_inc_pmwx_bench,Walk_inc_pmwx>(benchInst, benchable,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }
  }

  // Construct Aggregately
  typeId = CGAL::Bench_parse_args::TYPE_AGGREGATE;
  if (typeMask & (0x1 << typeId)) {
    // Dummy point location:
    StrategyId strategyId = CGAL::Bench_parse_args::STRATEGY_DUMMY;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Dummy_agg_pmwx_bench benchInst(name, seconds, false);
      Dummy_agg_pmwx & benchable = benchInst.get_benchable();
      runBench<Dummy_agg_pmwx_bench,Dummy_agg_pmwx>(benchInst, benchable,
                                                    fullname, format,
                                                    samples, iterations,
                                                    verbose);
    }
  }
  
  // Construct and Display
  typeId = CGAL::Bench_parse_args::TYPE_DISPLAY;
  if (typeMask & (0x1 << typeId)) {
    // Trapezoidal point location:
    StrategyId strategyId = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Trap_dis_pmwx_bench benchInst(name, seconds, false);
      Trap_dis_pmwx & benchable = benchInst.get_benchable();
      runBench<Trap_dis_pmwx_bench,Trap_dis_pmwx>(benchInst, benchable,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }

    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Naive_dis_pmwx_bench benchInst(name, seconds, false);
      Naive_dis_pmwx & benchable = benchInst.get_benchable();
      runBench<Naive_dis_pmwx_bench,Naive_dis_pmwx>(benchInst, benchable,
                                                    fullname, format,
                                                    samples, iterations,
                                                    verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Walk_dis_pmwx_bench benchInst(name, seconds, false);
      Walk_dis_pmwx & benchable = benchInst.get_benchable();
      runBench<Walk_dis_pmwx_bench,Walk_dis_pmwx>(benchInst, benchable,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }

    // Dummy point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_DUMMY;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.get_type_name(typeId)) + " " +
          std::string(parseArgs.get_strategy_name(strategyId)) + " " +
          "PMWX " + TRAITS_TYPE + " " + KERNEL_TYPE + " " + NUMBER_TYPE +
          " " + INSERT_TYPE + " (" + std::string(filename) + ")";
      Dummy_dis_pmwx_bench benchInst(name, seconds, false);
      Dummy_dis_pmwx & benchable = benchInst.get_benchable();
      runBench<Dummy_dis_pmwx_bench,Dummy_dis_pmwx>(benchInst, benchable,
                                                    fullname, format,
                                                    samples, iterations,
                                                    verbose);
    }
  }
  
  return 0;
}

