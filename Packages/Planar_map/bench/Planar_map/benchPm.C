#include <CGAL/config.h>

#include "short_names.h"

#include <CGAL/basic.h>
#include <CGAL/leda_rational.h>

#if defined(USE_LEDA_KERNEL)
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#else
#if defined(USE_MY_KERNEL)
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>
#else
#include <CGAL/Cartesian.h>
#endif
#endif

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Pm_Window_stream.h>

#include <CGAL/Pm_default_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_naive_point_location.h>

#include <CGAL/Bench.C>
#include <CGAL/Bench_parse_args.h>
#include <CGAL/Bench_parse_args.C>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>

typedef leda_rational                                   NT;

#if defined(USE_LEDA_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Kernel;
#define PM_TYPE "Leda Kernel"
#else
#if defined(USE_MY_KERNEL)
typedef CGAL::Pm_segment_traits_leda_kernel_2<NT>       Kernel;
#define PM_TYPE "My Leda Kernel"
#else
typedef CGAL::Cartesian<NT>                             Kernel;
#define PM_TYPE "CGAL Kernel"
#endif
#endif

typedef CGAL::Pm_segment_traits_2<Kernel>               Traits;
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Planar_map;
typedef Traits::Point_2                                 Point;
typedef Traits::X_curve_2                               Curve;
typedef std::list<Curve>                                CurveList;

typedef CGAL::Pm_default_point_location<Planar_map>     Trap_point_location;
typedef CGAL::Pm_naive_point_location<Planar_map>       Naive_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Planar_map>
                                                        Walk_point_location;

typedef CGAL::Bench_parse_args::TypeId                  TypeId;
typedef CGAL::Bench_parse_args::StrategyId              StrategyId;
typedef CGAL::Bench_parse_args::FormatId                FormatId;

/*
 */
class Basic_pm {
public:
  /*
   */
  Basic_pm() : m_filename(0), m_verbose(false),
    m_x0(0), m_x1(0), m_y0(0), m_y1(0) {}
  virtual ~Basic_pm() {}
  
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
    if (m_verbose) std::cout << m_curveList.size() << std::endl;

    return 0;
  }
    
  /*
   */
  void clean() { m_curveList.clear(); }
  void sync(){}

  void setFormat(FormatId fmt) { m_format = fmt; }
  void setFilename(const char * filename) { m_filename = filename; }
  void setVerbose(const bool verbose) { m_verbose = verbose; }

protected:
  const char * m_filename;
  CurveList m_curveList;
  bool m_verbose;
  FormatId m_format;

  double m_x0;
  double m_x1;
  double m_y0;
  double m_y1;
};

/*! Construct Incrementaly
 */
template <class Strategy>
class Increment_pm : public Basic_pm {
public:
  virtual void op()
  {
    Strategy strategy;
    Planar_map pm(&strategy);
    CurveList::const_iterator i;
    for (i = m_curveList.begin(); i != m_curveList.end(); i++) pm.insert(*i);
    // if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
    //std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
    //std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
    //std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
  }
};

typedef Increment_pm<Trap_point_location>       Trap_inc_pm;
typedef CGAL::Bench<Trap_inc_pm>                Trap_inc_pm_bench;

typedef Increment_pm<Naive_point_location>      Naive_inc_pm;
typedef CGAL::Bench<Naive_inc_pm>               Naive_inc_pm_bench;

typedef Increment_pm<Walk_point_location>       Walk_inc_pm;
typedef CGAL::Bench<Walk_inc_pm>                Walk_inc_pm_bench;

/*! Construct Aggregately
 */
template <class Strategy>
class Aggregate_pm : public Basic_pm {
public:
  virtual void op()
  {
    Strategy strategy;
    Planar_map pm(&strategy);
    pm.insert(m_curveList.begin(), m_curveList.end());
  }
};

typedef Aggregate_pm<Trap_point_location>       Trap_agg_pm;
typedef CGAL::Bench<Trap_agg_pm>                Trap_agg_pm_bench;

typedef Aggregate_pm<Naive_point_location>      Naive_agg_pm;
typedef CGAL::Bench<Naive_agg_pm>               Naive_agg_pm_bench;

typedef Aggregate_pm<Walk_point_location>       Walk_agg_pm;
typedef CGAL::Bench<Walk_agg_pm>                Walk_agg_pm_bench;

#if defined(USE_LEDA_KERNEL) || defined(USE_MY_KERNEL)
CGAL::Window_stream & operator<<(CGAL::Window_stream & os, const Point & p)
{ return os << leda_point(p.xcoordD(), p.ycoordD()); } 
CGAL::Window_stream & operator<<(CGAL::Window_stream & os, const Curve & c)
{ return os << leda_segment(c.xcoord1D(),c.ycoord1D(),c.xcoord2D(),c.ycoord2D()); }
#endif

/*!
 */
template <class Strategy>
class Display_pm : public Basic_pm {
private:
  typedef CGAL::Window_stream Window_stream;

public:
  /*!
   */
  virtual void op()
  {
    Strategy strategy;
    Planar_map pm(&strategy);
    pm.insert(m_curveList.begin(), m_curveList.end());
    m_window->set_flush(0);
    (*m_window) << pm;
    m_window->set_flush(1);
    m_window->flush();
  }
  
  /*!
   */
  int init()
  {
    int rc = Basic_pm::init();
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

    m_window->set_redraw(&Display_pm::redraw);
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
    Basic_pm::clean();
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

typedef Display_pm<Trap_point_location>         Trap_dis_pm;
typedef CGAL::Bench<Trap_dis_pm>                Trap_dis_pm_bench;

typedef Display_pm<Naive_point_location>        Naive_dis_pm;
typedef CGAL::Bench<Naive_dis_pm>               Naive_dis_pm_bench;

typedef Display_pm<Walk_point_location>         Walk_dis_pm;
typedef CGAL::Bench<Walk_dis_pm>                Walk_dis_pm_bench;

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
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Trap_inc_pm_bench benchInst(name, seconds, false);
      Trap_inc_pm & benchUser = benchInst.getBenchUser();
      runBench<Trap_inc_pm_bench,Trap_inc_pm>(benchInst, benchUser,
                                              fullname, format,
                                              samples, iterations, verbose);
    }
    
    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Naive_inc_pm_bench benchInst(name, seconds, false);
      Naive_inc_pm & benchUser = benchInst.getBenchUser();
      runBench<Naive_inc_pm_bench,Naive_inc_pm>(benchInst, benchUser,
                                                fullname, format,
                                                samples, iterations, verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Walk_inc_pm_bench benchInst(name, seconds, false);
      Walk_inc_pm & benchUser = benchInst.getBenchUser();
      runBench<Walk_inc_pm_bench,Walk_inc_pm>(benchInst, benchUser,
                                              fullname, format,
                                              samples, iterations, verbose);
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
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Trap_agg_pm_bench benchInst(name, seconds, false);
      Trap_agg_pm & benchUser = benchInst.getBenchUser();
      runBench<Trap_agg_pm_bench,Trap_agg_pm>(benchInst, benchUser,
                                              fullname, format,
                                              samples, iterations, verbose);
    }

    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Naive_agg_pm_bench benchInst(name, seconds, false);
      Naive_agg_pm & benchUser = benchInst.getBenchUser();
      runBench<Naive_agg_pm_bench,Naive_agg_pm>(benchInst, benchUser,
                                                fullname, format,
                                                samples, iterations, verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Walk_agg_pm_bench benchInst(name, seconds, false);
      Walk_agg_pm & benchUser = benchInst.getBenchUser();
      runBench<Walk_agg_pm_bench,Walk_agg_pm>(benchInst, benchUser,
                                              fullname, format,
                                              samples, iterations, verbose);
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
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Trap_dis_pm_bench benchInst(name, seconds, false);
      Trap_dis_pm & benchUser = benchInst.getBenchUser();
      runBench<Trap_dis_pm_bench,Trap_dis_pm>(benchInst, benchUser,
                                              fullname, format,
                                              samples, iterations, verbose);
    }

    // Naive point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_NAIVE;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Naive_dis_pm_bench benchInst(name, seconds, false);
      Naive_dis_pm & benchUser = benchInst.getBenchUser();
      runBench<Naive_dis_pm_bench,Naive_dis_pm>(benchInst, benchUser,
                                                fullname, format,
                                                samples, iterations, verbose);
    }

    // Walk point location:
    strategyId = CGAL::Bench_parse_args::STRATEGY_WALK;
    if (strategyMask & (0x1 << strategyId)) {
      std::string name =
          std::string(parseArgs.getTypeName(typeId)) + " " +
          std::string(parseArgs.getStrategyName(strategyId)) + " " +
          "PM " + PM_TYPE + " (" + std::string(filename) + ")";
      Walk_dis_pm_bench benchInst(name, seconds, false);
      Walk_dis_pm & benchUser = benchInst.getBenchUser();
      runBench<Walk_dis_pm_bench,Walk_dis_pm>(benchInst, benchUser,
                                              fullname, format,
                                              samples, iterations, verbose);
    }
  }
  
  return 0;
}
