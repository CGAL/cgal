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
#include <CGAL/IO/Planar_map_iostream.h>
#include "CGAL/bench.h"
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

static
#if (defined _MSC_VER)
const
#endif
char * BenchOpts[] = {
  "construct", "display", "co", "di", NULL
};

#define BENCH_CONSTRUCT 0
#define BENCH_DISPLAY   1
#define BENCH_CO        2
#define BENCH_DI        3
#define NUM_BENCHS      2

static char OptionStr[] = "b:hs:t:v";

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
      NT x, y;
      inp >> x >> y;
      double dx = CGAL::to_double(x);
      double dy = CGAL::to_double(y);
      if (i == 0) {
        m_x0 = dx;
        m_y0 = dy;
        m_x1 = dx;
        m_y1 = dy;
      } else {
        if (dx < m_x0) m_x0 = dx;
        if (dy < m_y0) m_y0 = dy;
        if (m_x1 < dx) m_x1 = dx;
        if (m_y1 < dy) m_y1 = dy;
      }
      Point p1(x,y);
      inp >> x >> y;
      Point p2(x,y);
      // if (p1 == p2) continue;
      Curve curve(p1, p2);
      m_curveList.push_back(curve);
    }
    inp.close();
    if (m_verbose) std::cout << m_curveList.size() << std::endl;

    return 0;
  }

  /*
   */
  void clean() { m_curveList.clear(); }
  void sync(){}

  void setFilename(const char * filename) { m_filename = filename; }
  void setVerbose(const bool verbose) { m_verbose = verbose; }

protected:
  const char * m_filename;
  CurveList m_curveList;
  bool m_verbose;

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
  typedef CGAL::Window_stream           Window_stream;
  typedef Planar_map::Vertex_iterator   Vertex_iterator;
  typedef Planar_map::Edge_iterator     Edge_iterator;

public:
  /*!
   */
  virtual void op()
  {
    Planar_map pm;
    pm.insert(m_curveList.begin(), m_curveList.end());
    m_window->set_flush(0);
    Vertex_iterator vi;
    for (vi = pm.vertices_begin(); vi != pm.vertices_end(); vi++) {
      (*m_window) << CGAL::GREEN;
      (*m_window) << (*vi).point();
    }

    Edge_iterator ei;
    for (ei = pm.edges_begin(); ei != pm.edges_end(); ei++) {
      (*m_window) << CGAL::BLUE;
      (*m_window) << ((*ei).curve());
    }

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

/*
 */
int main(int argc, char * argv[])
{
  progName = strrchr(argv[0], '\\');
  progName = (progName) ? progName+1 : argv[0];

  bool verbose = false;
  int benchId = -1;
  
#if (defined _MSC_VER)
  int samples = 10;
#else
  int samples = 0;
#endif
  
  int seconds = 0;
  int c;
  char * options, * value;
  while ((c = getopt(argc, argv, OptionStr)) != EOF) {
    switch (c) {
      case 'b':
	options = optarg;
	if (*options == '\0') break;
	while (*options != '\0') {
          switch(getsubopt(&options, BenchOpts, &value)) {
            case BENCH_CONSTRUCT:
            case BENCH_CO:
              benchId = BENCH_CONSTRUCT; break;
            case BENCH_DISPLAY:
            case BENCH_DI:
              benchId = BENCH_DISPLAY; break;
          }
        }
        break;
      case 'h': printHelp(); return 0;
      case 's': samples = atoi(optarg); break;
      case 't': seconds = atoi(optarg); break;
      case 'v': verbose = !verbose; break;
      default:
        fprintf(stderr, "%s: invalid option -- %c\n", progName, c);
        fprintf(stderr, "Try `%s -h' for more information.\n", progName);
        return -1;
    }
  }
  
  if (argc < 2) {
    std::cerr << "Data file missing!" << std::endl;
    return -1;
  }
  const char * filename = argv[optind];

  // Construct
  ConstructPmBench benchConstruct((std::string(BenchOpts[BENCH_CONSTRUCT]) +
                                  std::string(" PM (") +
                                  std::string(filename) + std::string(")")),
                                  seconds, true);
  Construct_Pm & construct_pm = benchConstruct.getBenchUser();
  construct_pm.setFilename(filename);
  construct_pm.setVerbose(verbose);

  // Construct and Display
  DisplayPmBench benchDisplay((std::string(BenchOpts[BENCH_DISPLAY]) +
                              std::string(" PM (") +
                              std::string(filename) + std::string(")")),
                              seconds, false);
  Display_Pm & display_pm = benchDisplay.getBenchUser();
  display_pm.setFilename(filename);
  display_pm.setVerbose(verbose);

  if (samples > 0) {
    benchConstruct.setSamples(samples);
    benchDisplay.setSamples(samples);
  }

  if (benchId == -1) {
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
