#include <CGAL/config.h>

#include "short_names.h"

#include <CGAL/basic.h>

#if defined(USE_CONIC_TRAITS)
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/leda_real.h>
#include <CGAL/IO/Conic_arc_2_Window_stream.h>
#include <CGAL/Cartesian.h>
#else
#include <CGAL/leda_rational.h>
#if defined(USE_LEDA_KERNEL)
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#if defined(USE_INSERT_FAST)
#include <CGAL/Arr_segment_traits_tight_2.h>
#else
#include <CGAL/Arr_segment_traits_2.h>
#endif
#else
#if defined(USE_MY_KERNEL)
#include <CGAL/Arr_leda_segment_traits_2.h>
#else
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#endif
#endif
#endif

#if defined(USE_INSERT_FAST)
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_tight_2.h>
#else
#include <CGAL/Sweep_line_2.h>
#endif

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

#if defined(USE_LEDA_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Kernel;
#if defined(USE_INSERT_FAST)
#define PM_TYPE "Tight Sweep (seg)"
typedef CGAL::Arr_segment_traits_tight_2<Kernel>         Traits;
#else
#define PM_TYPE "Sweep       (segs)"
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#endif

#else

#if defined(USE_MY_KERNEL)
typedef CGAL::Arr_leda_segment_traits_2                 Traits;
#define PM_TYPE "My Leda Kernel"
#else
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#define PM_TYPE "CGAL Kernel"
#endif
#endif
#endif

typedef Traits::Point_2                                 Point;
typedef Traits::X_curve_2                               Curve;
typedef std::list<Curve>                                CurveList;
typedef CurveList::iterator                             CurveListIter;

typedef CGAL::Bench_parse_args::TypeId                  TypeId;
typedef CGAL::Bench_parse_args::StrategyId              StrategyId;
typedef CGAL::Bench_parse_args::FormatId                FormatId;

#if defined(USE_INSERT_FAST)
typedef CGAL::Sweep_line_subcurve<Traits> CurveWrap;
typedef CGAL::Sweep_line_event<Traits, CurveWrap> SweepEvent;

typedef CGAL::Sweep_line_tight_2<CurveListIter, Traits, 
                                SweepEvent, CurveWrap> SweepLine;
#else
typedef CGAL::Sweep_line_2<CurveListIter, Traits> SweepLine;
#endif

/*
 */
class Basic_Sweep {
public:
  /*
   */
  Basic_Sweep() : m_filename(0), m_verbose(false), m_bbox(0.0,0.0,0.0,0.0) { }

  virtual ~Basic_Sweep() {}
  
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

class Subcurves_sweep : public Basic_Sweep {
public:
  virtual void op()
  {
    SweepLine sl;
    std::list<Curve> subcurves;
    sl.get_subcurves(m_curveList.begin(),
		     m_curveList.end(), 
		     std::back_inserter(subcurves), true);
    if (m_verbose) {
      std::cout << "# of subcurves: " << subcurves.size() << std::endl;
    }
  }
};
typedef CGAL::Bench<Subcurves_sweep> Subcurves_sweep_bench;

/*!
 */

class Points_sweep : public Basic_Sweep {
public:
  virtual void op()
  {
    SweepLine sl;
    std::vector<Point> points;
    sl.get_intersection_points(m_curveList.begin(),
			       m_curveList.end(), 
			       std::back_inserter(points), true);
    if (m_verbose) {
      std::cout << "# of subcurves: " << points.size() << std::endl;
    }
  }
};
/*
typedef Display_pmwx<Trap_point_location>     Trap_dis_pmwx;
typedef CGAL::Bench<Trap_dis_pmwx>            Trap_dis_pmwx_bench;

typedef Display_pmwx<Naive_point_location>    Naive_dis_pmwx;
typedef CGAL::Bench<Naive_dis_pmwx>           Naive_dis_pmwx_bench;

typedef Display_pmwx<Walk_point_location>     Walk_dis_pmwx;
typedef CGAL::Bench<Walk_dis_pmwx>            Walk_dis_pmwx_bench;
*/
typedef CGAL::Bench<Points_sweep> Points_sweep_bench;
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
  TypeId typeId = CGAL::Bench_parse_args::TYPE_SUBCURVES;
  if (typeMask & (0x1 << typeId)) {

    //std::string name = std::string(parseArgs.getTypeName(typeId)) + " " +
    //  "SL " + PM_TYPE + " (" + std::string(filename) + ")";
    std::string name = std::string(filename) + " - " + PM_TYPE + " " + std::string(parseArgs.getTypeName(typeId));

    Subcurves_sweep_bench benchInst(name, seconds, false);
    Subcurves_sweep & benchUser = benchInst.getBenchUser();
    runBench<Subcurves_sweep_bench, Subcurves_sweep>(benchInst, benchUser,
						     fullname, format,
						     samples, iterations,
						     verbose);
  }

  // Construct Aggregately
  typeId = CGAL::Bench_parse_args::TYPE_POINTS;
  if (typeMask & (0x1 << typeId)) {

    // std::string name = std::string(parseArgs.getTypeName(typeId)) + " " +
    //  "   SL " + PM_TYPE + " (" + std::string(filename) + ")";
    std::string name = std::string(filename) + " - " + PM_TYPE + " " + std::string(parseArgs.getTypeName(typeId));
    Points_sweep_bench benchInst(name, seconds, false);
    Points_sweep & benchUser = benchInst.getBenchUser();
    runBench<Points_sweep_bench, Points_sweep>(benchInst, benchUser,
					       fullname, format,
					       samples, iterations,
					       verbose);
  }

  return 0;
}

