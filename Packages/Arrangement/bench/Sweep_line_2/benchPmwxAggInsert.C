#include <CGAL/config.h>

#include "short_names.h"

#include <CGAL/basic.h>

#if defined(USE_CONIC_TRAITS)
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/leda_real.h>
//#include <CGAL/IO/Conic_arc_2_Window_stream.h>
#else
#include <CGAL/leda_rational.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#if defined(USE_INSERT_FAST)
#include <CGAL/Arr_segment_traits_tight_2.h>
#else
#include <CGAL/Arr_segment_traits_2.h>
#endif
#endif

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
//#include <CGAL/IO/Window_stream.h>
//#include <CGAL/IO/Pm_iostream.h>
//#include <CGAL/IO/Pm_Window_stream.h>

#include <CGAL/Pm_naive_point_location.h>

#include <CGAL/Bench.h>
#include <CGAL/Bench.C>
#include <CGAL/Bench_parse_args.h>
#include <CGAL/Bench_parse_args.C>

#include <stdlib.h>
#include <iostream>
#include <list>

#if defined(USE_CONIC_TRAITS)
#include "Conic_reader.h"
#else
#include "Segment_reader.h"
#endif

#if defined(USE_CONIC_TRAITS)
typedef leda_real                                       NT;
#else
typedef leda_rational                                   NT;
#endif

#if defined(USE_CONIC_TRAITS)
typedef CGAL::Cartesian<NT>               Kernel;
typedef CGAL::Arr_conic_traits_2<Kernel>  Traits;

#if defined(USE_INSERT_FAST)
#define PM_TYPE "Tight Sweep (conics)"
#else
#define PM_TYPE "Sweep       (conics)"
#endif

#else

typedef CGAL::leda_rat_kernel_traits                    Kernel;
#if defined(USE_INSERT_FAST)
typedef CGAL::Arr_segment_traits_tight_2<Kernel>        Traits;
#define PM_TYPE "Tight Sweep (seg)"
#else
typedef CGAL::Arr_segment_traits_2<Kernel>          Traits;
#define PM_TYPE "Sweep       (segs)"
#endif

#endif

typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;
typedef Traits::Point_2                                 Point;
typedef Traits::X_curve_2                               Curve;
typedef std::list<Curve>                                CurveList;

typedef CGAL::Pm_naive_point_location<Pm>               NPL;

typedef CGAL::Bench_parse_args::FormatId                FormatId;


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

  void setFormat(CGAL::Bench_parse_args::FormatId fmt) { m_format = fmt; }
  void setFilename(const char * filename) { m_filename = filename; }
  void setVerbose(const bool verbose) { m_verbose = verbose; }

protected:
  const char * m_filename;
  CurveList m_curveList;
  bool m_verbose;
  CGAL::Bench_parse_args::FormatId m_format;
  CGAL::Bbox_2 m_bbox;
};


/*!
 */
class Aggregate_pmwx : public Basic_Pm {
public:
  virtual void op()
  {
    if (m_verbose) {
      std::cout << "Inserting Aggregate" << std::endl;
    }

    NPL pl;
    Pmwx pm(&pl); 
#if defined(USE_INSERT_FAST)
    pm.insert_tight(m_curveList.begin(), m_curveList.end());
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

typedef CGAL::Bench<Aggregate_pmwx>             Agg_pmwx_bench;


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
  //unsigned int typeMask = parseArgs.getTypeMask();
  //unsigned int strategyMask = parseArgs.getStrategyMask();
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
  
  
  // Naive point location:
  std::string name = std::string(filename) + " - " + std::string("Agg. insert, ") + 
    PM_TYPE;
  Agg_pmwx_bench benchInst(name, seconds, false);
  Aggregate_pmwx & benchUser = benchInst.getBenchUser();
  runBench<Agg_pmwx_bench,Aggregate_pmwx>(benchInst, benchUser,
				    fullname, format,
				    samples, iterations,
				    verbose);
  
  return 0;
}

