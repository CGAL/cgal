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
#include <CGAL/IO/Pm_iostream.h>
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


static char OptionStr[] = "hs:t:v";

/*
 */
class Basic_Pm {
public:
  /*
   */
  Basic_Pm() : m_filename(0), m_verbose(false) {}
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

/*!
 */
static void printHelp(void)
{
  printf("Usage: bench [options]\n\
  -h\t\tprint this help message\n\
  -s <samples>\tset number of samples to <samples>\n\
  -t <seconds>\tset number of seconds to <seconds>\n");
}

/*
 */
int main(int argc, char * argv[])
{
  char * progName = strrchr(argv[0], '\\');
  progName = (progName) ? progName+1 : argv[0];

  bool verbose = false;
#if (defined _MSC_VER)
  int samples = 10;
#else
  int samples = 0;
#endif
  
  int seconds = 0;
  int c;
  while ((c = getopt(argc, argv, OptionStr)) != EOF) {
    switch (c) {
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

  ConstructPmBench bench(std::string("Construct PM (") +
                         std::string(filename) + std::string(")"));
  Construct_Pm & construct_pm = bench.getBenchUser();
  construct_pm.setFilename(filename);
  construct_pm.setVerbose(verbose);

  if (samples > 0) bench.setSamples(samples);
  if (seconds > 0) bench.setSeconds(seconds);

  bench();

  // Ensure the compiler doesn't optimize the code away...
  if (verbose) std::cout << "(" << bench.getIterations() << ") " << std::endl;
  
  return 0;
}
