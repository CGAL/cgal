#include <CGAL/config.h>

#include "CGAL/Short_names.h"

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
#include "CGAL/Bench.h"
#include "CGAL/Dir_search.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>

#include "getopt.h"

typedef leda_rational                                   NT;

#if defined(USE_LEDA_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Kernel;
#else
#if defined(USE_MY_KERNEL)
typedef CGAL::Pm_segment_traits_leda_kernel_2<NT>       Kernel;
#else
typedef CGAL::Cartesian<NT>                             Kernel;
#endif
#endif

typedef CGAL::Pm_segment_traits_2<Kernel>               Traits;
typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Planar_map;
typedef Traits::Point_2                                 Point;
typedef Traits::X_curve_2                               Curve;
typedef std::list<Curve>                                CurveList;

// Global variable declarations:
static
#if (defined _MSC_VER)
const
#endif
char * IOOpts[] = {
  "format", "f", NULL
};

enum IOId  {
  IO_FORMAT = 0,
  IO_F
};

static
#if (defined _MSC_VER)
const
#endif
char * FormatOpts[] = {
  "rat", "int", "float", "r", "i", "f", NULL
};

enum FormatId  {
  FORMAT_RAT = 0,
  FORMAT_INT,
  FORMAT_FLT,
  FORMAT_R,
  FORMAT_I,
  FORMAT_F
};

static
#if (defined _MSC_VER)
const
#endif
char * BenchOpts[] = {
  "construct", "display", "co", "di", NULL
};

enum BenchId  {
  BENCH_CONSTRUCT = 0,
  BENCH_DISPLAY,
  BENCH_CO,
  BENCH_DI,
  BENCH_ALL
};

static char OptionStr[] = "b:hi:s:t:v";

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
      if (m_format == FORMAT_RAT) {
        inp >> x0 >> y0 >> x1 >> y1;
      } else if (m_format == FORMAT_INT) {
        int ix0, iy0, ix1, iy1;
        inp >> ix0 >> iy0 >> ix1 >> iy1;
        x0 = ix0; y0 = iy0; x1 = ix1; y1 = iy1;
      } else if (m_format == FORMAT_FLT) {
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

  void setFormat(FormatId format) { m_format = format; }
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

/*!
 */
class Construct_Pm : public Basic_Pm {
public:
  virtual void op()
  {
    Planar_map pm;
    pm.insert(m_curveList.begin(), m_curveList.end());
    // if (!pm.is_valid()) std::cerr << "map invalid!" << std::endl;
    //std::cout << "# of vertices: " << pm.number_of_vertices() << std::endl;
    //std::cout << "# of halfedges: " << pm.number_of_halfedges() << std::endl;
    //std::cout << "# of faces: " << pm.number_of_faces() << std::endl;
  }
};

typedef CGAL::Bench<Construct_Pm> ConstructPmBench;

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
    Planar_map pm;
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

static char * progName = 0;

/*!
 */
static void printHelp(void)
{
  printf("Usage: %s [options]\n\
  -b <name>\tset bench name to <name> (default all)\n\
  \t\t\t<name> is one of {co[nstruct], di[splay]}\n\
  -h\t\tprint this help message\n\
  -s <samples>\tset number of samples to <samples> (default 10)\n\
  -t <seconds>\tset number of seconds to <seconds> (default 1)\n\
  -v\t\ttoggle verbosity (default no)\n", progName);
}

/*!
 */
int getFormat(FormatId & format, char * value)
{
  char * subvalue;
  switch (getsubopt(&value, FormatOpts, &subvalue)) {
    case FORMAT_INT:
    case FORMAT_I: format = FORMAT_INT; break;
    case FORMAT_RAT:
    case FORMAT_R: format = FORMAT_RAT; break;
    case FORMAT_FLT:
    case FORMAT_F: format = FORMAT_FLT; break;
    default:
      std::cerr << "Unrecognized Format option '" << optarg << "'!"
                << std::endl;
      std::cerr << "Try `" << progName << " -h' for more information."
                << std::endl;
      return -1;
  }
  return 0;
}

/*!
 */
int getIOParm(FormatId & format, char * optarg)
{
  char * options = optarg;
  char * value;
  if (*options == '\0') return 0;
  while (*options != '\0') {
    switch (getsubopt(&options, IOOpts, &value)) {
      case IO_FORMAT:
      case IO_F: if (getFormat(format, value) < 0) return -1; break;
      default:
        std::cerr << "Unrecognized IO option '" << optarg << "'!" << std::endl;
        std::cerr << "Try `" << progName << " -h' for more information."
                  << std::endl;
        return -1;
    }
  }
  return 0;
}

/*!
 */
int getBenchParam(BenchId & benchId, char * optarg)
{
  char * options = optarg;
  char * value;
  if (*options == '\0') return 0;
  while (*options != '\0') {
    switch(getsubopt(&options, BenchOpts, &value)) {
      case BENCH_CONSTRUCT:
      case BENCH_CO:
        benchId = BENCH_CONSTRUCT; break;
      case BENCH_DISPLAY:
      case BENCH_DI:
        benchId = BENCH_DISPLAY; break;
      default:
        std::cerr << "Unrecognized Bench option '" << optarg << "'!"
                  << std::endl;
        std::cerr << "Try `" << progName << " -h' for more information."
                  << std::endl;
        return -1;
    }
  }
  return 0;
}

/*
 */
int main(int argc, char * argv[])
{
  progName = strrchr(argv[0], '\\');
  progName = (progName) ? progName+1 : argv[0];

  bool verbose = false;
  BenchId benchId = BENCH_ALL;
  FormatId inputFormat = FORMAT_RAT;

  Dir_search dirs(".");
  const char * root = getenv("ROOT");
  if (root) dirs.add(std::string(root) + "/data/Segments_2");
      
#if (defined _MSC_VER)
  int samples = 10;
#else
  int samples = 0;
#endif
  
  int seconds = 0;
  int c;
  while ((c = getopt(argc, argv, OptionStr)) != EOF) {
    switch (c) {
      case 'b': if (getBenchParam(benchId, optarg) < 0) return -1; break;
      case 'h': printHelp(); return 0;
      case 'i': if (getIOParm(inputFormat, optarg)) return -1; break;
      case 's': samples = atoi(optarg); break;
      case 't': seconds = atoi(optarg); break;
      case 'v': verbose = !verbose; break;
      default:
        std::cerr << progName << ": invalid option -- "
                  << static_cast<char>(c) << std::endl;
        std::cerr << "Try `" << progName << " -h' for more information."
                  << std::endl;
        return -1;
    }
  }
  
  if (argc < 2) {
    std::cerr << "Data file missing!" << std::endl;
    return -1;
  }

  const char * filename = argv[optind];
  std::string fullname;
  if (!dirs.find(filename, fullname)) {
    std::cerr << "Cannot find file " << filename << "!" << std::endl;
    return -1;
  }
      
  // Construct
  ConstructPmBench benchConstruct((std::string(BenchOpts[BENCH_CONSTRUCT]) +
                                  " PM (" + std::string(filename) + ")"),
                                  seconds, true);
  Construct_Pm & construct_pm = benchConstruct.getBenchUser();
  construct_pm.setFormat(inputFormat);
  construct_pm.setFilename(fullname.c_str());
  construct_pm.setVerbose(verbose);

  // Construct and Display
  DisplayPmBench benchDisplay((std::string(BenchOpts[BENCH_DISPLAY]) +
                               " PM (" + std::string(filename) + ")"),
                              seconds, false);
  Display_Pm & display_pm = benchDisplay.getBenchUser();
  display_pm.setFormat(inputFormat);
  display_pm.setFilename(fullname.c_str());
  display_pm.setVerbose(verbose);

  if (samples > 0) {
    benchConstruct.setSamples(samples);
    benchDisplay.setSamples(samples);
  }

  if (benchId == BENCH_ALL) {
    benchConstruct();
    benchDisplay();
  } else switch(benchId) {
    case BENCH_CONSTRUCT: benchConstruct(); break;
    case BENCH_DISPLAY: benchDisplay(); break;
  }
  
  // Ensure the compiler doesn't optimize the code away...
  if (verbose) {
    std::cout << "(" << benchConstruct.getIterations() << ") " << std::endl;
    std::cout << "(" << benchDisplay.getIterations() << ") " << std::endl;
  }
  
  return 0;
}
