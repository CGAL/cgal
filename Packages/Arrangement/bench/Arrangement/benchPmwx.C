#include <CGAL/config.h>

#include "short_names.h"

#include <CGAL/basic.h>
#include <CGAL/leda_rational.h>

#if defined(USE_LEDA_KERNEL)
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#else
#if defined(USE_MY_KERNEL)
#include <CGAL/Arr_leda_segment_exact_traits.h>
#else
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_exact_traits.h>
#endif
#endif

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>

#include <CGAL/Pm_default_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>

#include <CGAL/Bench.h>
#include <CGAL/Bench.C>
#include <CGAL/Bench_parse_args.h>
#include <CGAL/Bench_parse_args.C>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>

typedef leda_rational                                   NT;

#if defined(USE_LEDA_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Traits;
#define PM_TYPE "Leda Kernel"
#else
#if defined(USE_MY_KERNEL)
typedef CGAL::Arr_leda_segment_exact_traits<NT>         Traits;
#define PM_TYPE "My Leda Kernel"
#else
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>          Traits;
#define PM_TYPE "CGAL Kernel"
#endif
#endif

typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;
typedef Traits::Point_2                                 Point;
typedef Traits::X_curve_2                               Curve;
typedef std::list<Curve>                                CurveList;

typedef CGAL::Pm_default_point_location<Pm>             Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm>               Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm>     Walk_point_location;

typedef CGAL::Bench_parse_args::TypeId                  TypeId;
typedef CGAL::Bench_parse_args::StrategyId              StrategyId;
typedef CGAL::Bench_parse_args::FormatId                FormatId;

/*
 */
class Basic_Pm {
public:
  /*
   */
  Basic_Pm() : m_filename(0), m_verbose(false),
    m_x0(0), m_x1(0), m_y0(0), m_y1(0) {}
  virtual ~Basic_Pm() {}
  
  /*
   */
  virtual void op() = 0;

  /*
   */
  void setExtreme(double x, double y)
  {
    m_x0 = x;
    m_y0 = y;
    m_x1 = x;
    m_y1 = y;
  }

  /*
   */
  void compareExtreme(double x, double y)
  {
    if (x < m_x0) m_x0 = x;
    if (y < m_y0) m_y0 = y;
    if (m_x1 < x) m_x1 = x;
    if (m_y1 < y) m_y1 = y;
  }
    
  /*
   */
  int init()
  {
    std::ifstream inp(m_filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << m_filename << "!" << std::endl;
      return -1;
    }
    int count;
    inp >> count;
    
    int i;
    for (i = 0; i < count; i++) {
      leda_rational x0, y0, x1, y1;
      if (m_format == CGAL::Bench_parse_args::FORMAT_RAT) {
        inp >> x0 >> y0 >> x1 >> y1;
      } else if (m_format == CGAL::Bench_parse_args::FORMAT_INT) {
        int ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else if (m_format == CGAL::Bench_parse_args::FORMAT_FLT) {
        float ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else {
        std::cerr << "Illegal format!" << std::endl;
        return -1;
      }

      Point p1(x0, y0);
      Point p2(x1, y1);
      // if (p1 == p2) continue;
      Curve curve(p1, p2);
      m_curveList.push_back(curve);
      if (i == 0) setExtreme(CGAL::to_double(x0), CGAL::to_double(y0));
      else compareExtreme(CGAL::to_double(x0), CGAL::to_double(y0));
      compareExtreme(CGAL::to_double(x1), CGAL::to_double(y1));
    }
    inp.close();
    if (m_verbose) std::cout << m_curveList.size() << " curves" << std::endl;

    return 0;
  }
    
  /*
   */
  void clean() { m_curveList.clear(); }
  void sync(){}

  void setFormat(CGAL::Bench_parse_args::FormatId fmt) { m_format = fmt; }
  void setFilename(const char * filename) { m_filename = filename; }
  void setVerbose(const bool verbose) { m_verbose = verbose; }

protected:
  const char * m_filename;
  CurveList m_curveList;
  bool m_verbose;
  CGAL::Bench_parse_args::FormatId m_format;

  double m_x0;
  double m_x1;
  double m_y0;
  double m_y1;
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
    // pm.insert(m_curveList.begin(), m_curveList.end());
    // if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
    //std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
    //std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
    //std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
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
    Strategy strategy;
    Pmwx pm(&strategy);
    pm.insert(m_curveList.begin(), m_curveList.end());
    // if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
  }
};

typedef Aggregate_pmwx<Trap_point_location>     Trap_agg_pmwx;
typedef CGAL::Bench<Trap_agg_pmwx>              Trap_agg_pmwx_bench;

typedef Aggregate_pmwx<Naive_point_location>    Naive_agg_pmwx;
typedef CGAL::Bench<Naive_agg_pmwx>             Naive_agg_pmwx_bench;

typedef Aggregate_pmwx<Walk_point_location>     Walk_agg_pmwx;
typedef CGAL::Bench<Walk_agg_pmwx>              Walk_agg_pmwx_bench;

#if defined(USE_LEDA_KERNEL) || defined(USE_MY_KERNEL)
CGAL::Window_stream & operator<<(CGAL::Window_stream & os, const Point & p)
{ return os << leda_point(p.xcoordD(), p.ycoordD()); } 
CGAL::Window_stream & operator<<(CGAL::Window_stream & os, const Curve & c)
{ return os << leda_segment(c.xcoord1D(),c.ycoord1D(),c.xcoord2D(),c.ycoord2D()); }
#endif

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
    Pm pm(&strategy);
    Traits traits;
    CGAL::sweep_to_construct_planar_map_2(m_curveList.begin(),
                                          m_curveList.end(),
                                          traits, pm);
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

    float x_range = m_x1 - m_x0;
    float y_range = m_y1 - m_y0;
    float width = 640;
    float height = (y_range * width) / x_range;
    
    m_window = new Window_stream(static_cast<int>(width),
                                 static_cast<int>(height));
    if (!m_window) return -1;

    float min_range = (x_range < y_range) ? x_range : y_range;
    float x_margin = min_range / 4;
    float y_margin = (height * x_margin) / width;
        
    float x0 = m_x0 - x_margin;
    float x1 = m_x1 + x_margin;
    float y0 = m_y0 - y_margin;
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
    Basic_Pm::init();
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

/*
 */
template <class Bench_inst, class Bench_user>
void runBench(Bench_inst & benchInst, Bench_user & benchUser,
              const char * fullname, FormatId format,
              int samples, int iterations, bool verbose)
{
    // Bench_inst benchInst(name, seconds, false);
    // Bench_user & benchUser = benchInst.getBenchUser();
  benchUser.setFormat(format);
  benchUser.setFilename(fullname);
  benchUser.setVerbose(verbose);

  if (samples > 0) benchInst.setSamples(samples);
  else if (iterations > 0) benchInst.setIterations(iterations);

  benchInst();

  if (verbose) std::cout << "(" << benchInst.getSamples() << ") " << std::endl;
}

/*
 */
int main(int argc, char * argv[])
{
  CGAL::Bench_parse_args parseArgs(argc, argv);
  int rc = parseArgs.parse();
  if (rc > 0) return 0;
  if (rc < 0) return rc;
  
  bool verbose = parseArgs.getVerbose();
  unsigned int typeMask = parseArgs.getTypeMask();
  unsigned int strategyMask = parseArgs.getStrategyMask();
  FormatId format = parseArgs.getInputFormat();
  int samples = parseArgs.getSamples();
  int iterations = parseArgs.getIterations();
  int seconds = parseArgs.getSeconds();
  bool printHeader = parseArgs.getPrintHeader();
  int nameLength = parseArgs.getNameLength();
  const char * filename = parseArgs.getFilename();
  const char * fullname = parseArgs.getFullname();
      
  CGAL::Bench_base::setNameLength(nameLength);
  if (printHeader) CGAL::Bench_base::printHeader();
  
  // Construct Incrementaly
  TypeId typeId = CGAL::Bench_parse_args::TYPE_INCREMENT;
  if (typeMask & (0x1 << typeId)) {
    // Trapezoidal point location:
    StrategyId strategyId = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Trap_inc_pmwx_bench benchInst(name, seconds, false);
      Trap_inc_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Trap_inc_pmwx_bench,Trap_inc_pmwx>(benchInst, benchUser,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }
    
    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Naive_inc_pmwx_bench benchInst(name, seconds, false);
      Naive_inc_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Naive_inc_pmwx_bench,Naive_inc_pmwx>(benchInst, benchUser,
                                                    fullname, format,
                                                    samples, iterations,
                                                    verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Walk_inc_pmwx_bench benchInst(name, seconds, false);
      Walk_inc_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Walk_inc_pmwx_bench,Walk_inc_pmwx>(benchInst, benchUser,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }
  }

  // Construct Aggregately
  typeId = CGAL::Bench_parse_args::TYPE_AGGREGATE;
  if (typeMask & (0x1 << typeId)) {
    // Trapezoidal point location:
    StrategyId strategyId = CGAL::Bench_parse_args::STRATEGY_TRAPEZOIDAL;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Trap_agg_pmwx_bench benchInst(name, seconds, false);
      Trap_agg_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Trap_agg_pmwx_bench,Trap_agg_pmwx>(benchInst, benchUser,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }

    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Naive_agg_pmwx_bench benchInst(name, seconds, false);
      Naive_agg_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Naive_agg_pmwx_bench,Naive_agg_pmwx>(benchInst, benchUser,
                                                    fullname, format,
                                                    samples, iterations,
                                                    verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Walk_agg_pmwx_bench benchInst(name, seconds, false);
      Walk_agg_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Walk_agg_pmwx_bench,Walk_agg_pmwx>(benchInst, benchUser,
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
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Trap_dis_pmwx_bench benchInst(name, seconds, false);
      Trap_dis_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Trap_dis_pmwx_bench,Trap_dis_pmwx>(benchInst, benchUser,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }

    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Naive_dis_pmwx_bench benchInst(name, seconds, false);
      Naive_dis_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Naive_dis_pmwx_bench,Naive_dis_pmwx>(benchInst, benchUser,
                                                    fullname, format,
                                                    samples, iterations,
                                                    verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PMWX " + PM_TYPE + " (" + std::string(filename) + ")";
      Walk_dis_pmwx_bench benchInst(name, seconds, false);
      Walk_dis_pmwx & benchUser = benchInst.getBenchUser();
      runBench<Walk_dis_pmwx_bench,Walk_dis_pmwx>(benchInst, benchUser,
                                                  fullname, format,
                                                  samples, iterations,
                                                  verbose);
    }
  }
    
  return 0;
}

