#include <CGAL/config.h>

#include "CGAL/Short_names.h"

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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>

#include "CGAL/Bench.h"
#include "CGAL/Bench.C"
#include "CGAL/Bench_parse_args.h"
#include "CGAL/Bench_parse_args.C"

typedef leda_rational                                   NT;

#if defined(USE_LEDA_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Kernel;
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
class Increment_pmwx : public Basic_Pm {
public:
  virtual void op()
  {
    Pmwx pm;
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

typedef CGAL::Bench<Increment_pmwx> Increemnt_pmwx_bench;

/*!
 */
class Aggregate_pm : public Basic_Pm {
public:
  virtual void op()
  {
    Pm pm;
    Traits traits;
    CGAL::sweep_to_construct_planar_map_2(m_curveList.begin(),
                                          m_curveList.end(),
                                          traits, pm);
    // if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
  }
};

typedef CGAL::Bench<Aggregate_pm> Aggregate_pm_bench;

#if defined(USE_LEDA_KERNEL) || defined(USE_MY_KERNEL)
CGAL::Window_stream & operator<<(CGAL::Window_stream & os, const Point & p)
{ return os << leda_point(p.xcoordD(), p.ycoordD()); } 
CGAL::Window_stream & operator<<(CGAL::Window_stream & os, const Curve & c)
{ return os << leda_segment(c.xcoord1D(),c.ycoord1D(),c.xcoord2D(),c.ycoord2D()); }
#endif

/*!
 */
class Display_Pm : public Basic_Pm {
private:
  typedef CGAL::Window_stream Window_stream;

public:
  /*!
   */
  virtual void op()
  {
    Pm pm;
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

    m_window->set_redraw(&Display_Pm::redraw);
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

typedef CGAL::Bench<Display_Pm> DisplayPmBench;

/*
 */
int main(int argc, char * argv[])
{
  CGAL::Bench_parse_args parseArgs(argc, argv);
  int rc = parseArgs.parse();
  if (rc > 0) return 0;
  if (rc < 0) return rc;
  
  bool verbose = parseArgs.getVerbose();
  unsigned int benchMask = parseArgs.getBenchMask();
  CGAL::Bench_parse_args::FormatId inputFormat = parseArgs.getInputFormat();
  int samples = parseArgs.getSamples();
  int iterations = parseArgs.getIterations();
  int seconds = parseArgs.getSeconds();
  bool printHeader = parseArgs.getPrintHeader();
  int nameLength = parseArgs.getNameLength();
  const char * filename = parseArgs.getFilename();
  const std::string * fullname = parseArgs.getFullname();
      
  // Construct
  const char * bname =
    parseArgs.getBenchName(CGAL::Bench_parse_args::BENCH_INCREMENT);
  Increemnt_pmwx_bench benchIncrement((std::string(bname) + " PMWX " +PM_TYPE +
                                       "(" + std::string(filename) + ")"),
                                      seconds, false);
  Increment_pmwx & incerementPmwx = benchIncrement.getBenchUser();
  incerementPmwx.setFormat(inputFormat);
  incerementPmwx.setFilename(fullname->c_str());
  incerementPmwx.setVerbose(verbose);

  // Sweep
  bname = parseArgs.getBenchName(CGAL::Bench_parse_args::BENCH_AGGREGATE);
  Aggregate_pm_bench benchAggregate((std::string(bname) + " PMWX " PM_TYPE +
                                     "(" + std::string(filename) + ")"),
                                    seconds, false);
  Aggregate_pm & aggregatePm = benchAggregate.getBenchUser();
  aggregatePm.setFormat(inputFormat);
  aggregatePm.setFilename(fullname->c_str());
  aggregatePm.setVerbose(verbose);
  
  // Construct and Display
  bname = parseArgs.getBenchName(CGAL::Bench_parse_args::BENCH_DISPLAY);
  DisplayPmBench benchDisplay((std::string(bname) + " PMWX " + PM_TYPE +
                               "(" + std::string(filename) + ")"),
                              seconds, false);
  Display_Pm & displayPm = benchDisplay.getBenchUser();
  displayPm.setFormat(inputFormat);
  displayPm.setFilename(fullname->c_str());
  displayPm.setVerbose(verbose);

  if (samples > 0) {
    benchIncrement.setSamples(samples);
    benchAggregate.setSamples(samples);
    benchDisplay.setSamples(samples);
  } else {
    if (iterations > 0) {
      benchIncrement.setIterations(iterations);
      benchAggregate.setIterations(iterations);
      benchDisplay.setIterations(iterations);
    }
  }

  CGAL::Bench_base::setNameLength(nameLength);
  if (printHeader) CGAL::Bench_base::printHeader();
  if (benchMask & (0x1 << CGAL::Bench_parse_args::BENCH_INCREMENT))
    benchIncrement();
  if (benchMask & (0x1 << CGAL::Bench_parse_args::BENCH_AGGREGATE))
    benchAggregate();
  if (benchMask & (0x1 << CGAL::Bench_parse_args::BENCH_DISPLAY))
    benchDisplay();
  
  // Ensure the compiler doesn't optimize the code away...
  if (verbose) {
    std::cout << "(" << benchIncrement.getSamples() << ")" << std::endl;
    std::cout << "(" << benchAggregate.getSamples() << ") " << std::endl;
    std::cout << "(" << benchDisplay.getSamples() << ") " << std::endl;
  }
  
  return 0;
}

