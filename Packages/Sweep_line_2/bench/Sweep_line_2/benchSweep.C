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

#include <CGAL/Arr_segment_traits_2.h>

#else
#if defined(USE_MY_KERNEL)
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/Arr_leda_segment_traits_2.h>
#else
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#endif
#endif
#endif

#if defined(USE_INSERT_OLD)
#include <CGAL/Sweep_line_2_old.h>
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

#if defined(USE_INSERT_OLD)
#define PM_TYPE "Old Sweep (conics)"
#else
#define PM_TYPE "Sweep     (conics)"
#endif

#else

#if defined(USE_LEDA_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#if defined(USE_INSERT_OLD)
#define PM_TYPE "Old Sweep (seg)"
#else
#define PM_TYPE "Sweep     (segs)"
#endif

#else

#if defined(USE_MY_KERNEL)
typedef CGAL::leda_rat_kernel_traits                    Kernel;
typedef CGAL::Arr_leda_segment_traits_2<Kernel>         Traits;
#define PM_TYPE "My Leda Kernel"
#else
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
#define PM_TYPE "CGAL Kernel"
#endif
#endif
#endif

typedef Traits::Point_2                                 Point;
typedef Traits::X_monotone_curve_2                      Curve;
typedef std::list<Curve>                                CurveList;
typedef CurveList::iterator                             CurveListIter;

typedef CGAL::Bench_parse_args::TypeId                  TypeId;
typedef CGAL::Bench_parse_args::StrategyId              StrategyId;
typedef CGAL::Bench_parse_args::FormatId                FormatId;

#if defined(USE_INSERT_OLD)
typedef CGAL::Sweep_line_2_old<CurveListIter, Traits> SweepLine;
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
template <class Bench_inst, class Benchable>
void runBench(Bench_inst & benchInst, Benchable & benchable,
              const char * fullname, FormatId format,
              int samples, int iterations, bool verbose)
{
    // Bench_inst benchInst(name, seconds, false);
    // Benchable & benchable = benchInst.getBenchUser();
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
  FormatId format = parseArgs.get_input_format();
  int samples = parseArgs.get_samples();
  int iterations = parseArgs.get_iterations();
  int seconds = parseArgs.get_seconds();
  bool printHeader = parseArgs.get_print_header();
  int nameLength = parseArgs.get_name_length();
  const char * filename = parseArgs.get_file_name();
  const char * fullname = parseArgs.get_full_name();
      
  CGAL::Bench_base::set_name_length(nameLength);
  if (printHeader) CGAL::Bench_base::print_header();
  

  // Construct Incrementaly
  TypeId typeId = CGAL::Bench_parse_args::TYPE_SUBCURVES;
  if (typeMask & (0x1 << typeId)) {

    //std::string name = std::string(parseArgs.get_type_name(typeId)) + " " +
    //  "SL " + PM_TYPE + " (" + std::string(filename) + ")";
    std::string name = std::string(filename) + " - " + PM_TYPE + " " +
      std::string(parseArgs.get_type_name(typeId));

    Subcurves_sweep_bench benchInst(name, seconds, false);
    Subcurves_sweep & benchable = benchInst.get_benchable();
    runBench<Subcurves_sweep_bench, Subcurves_sweep>(benchInst, benchable,
						     fullname, format,
						     samples, iterations,
						     verbose);
  }

  // Construct Aggregately
  typeId = CGAL::Bench_parse_args::TYPE_POINTS;
  if (typeMask & (0x1 << typeId)) {

    // std::string name = std::string(parseArgs.get_type_name(typeId)) + " " +
    //  "   SL " + PM_TYPE + " (" + std::string(filename) + ")";
    std::string name = std::string(filename) + " - " + PM_TYPE + " " +
      std::string(parseArgs.get_type_name(typeId));
    Points_sweep_bench benchInst(name, seconds, false);
    Points_sweep & benchable = benchInst.get_benchable();
    runBench<Points_sweep_bench, Points_sweep>(benchInst, benchable,
					       fullname, format,
					       samples, iterations,
					       verbose);
  }

  return 0;
}

