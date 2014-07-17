#include <iostream>
#include <sstream>
#include <cassert>

//#define CGAL_SDG_VERBOSE

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

// add this to profile
//#define CGAL_PROFILE 1

#define CGAL_SDG_NO_FACE_MAP 1
#define USE_INPLACE_LIST 1
#define CGAL_SDG_USE_SIMPLIFIED_ARRANGEMENT_TYPE_PREDICATE 1
//#define CGAL_SDG_SORT_POINTS_IN_SITE2 1

// choose the kernel
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#ifdef CGAL_USE_LEDA
#  include <CGAL/Leda_real.h>
#else
//#  include <CGAL/CORE_Expr.h>
#  include <CGAL/Gmpq.h>
#endif

typedef CGAL::Simple_cartesian<double> K;
typedef  K::Point_2 Point_2;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>

typedef CGAL::Field_with_sqrt_tag MTag;
#ifdef CGAL_USE_LEDA
typedef CGAL::Field_with_sqrt_tag EMTag;
typedef CGAL::Simple_cartesian<leda::real> EK;
#else
//typedef CGAL::Field_with_sqrt_tag EMTag;
//typedef CGAL::Simple_cartesian<CORE::Expr> EK;
typedef CGAL::Integral_domain_without_division_tag EMTag;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> EK;
#endif
typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2
        <K, MTag, EK, EMTag>  GtLinf;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2
        <K, MTag, EK, EMTag>  GtL2;

#include <tclap/CmdLine.h>


template <typename Gt>
void
run_benchmark(const size_t repetitions, const typename Gt::Site_2 & p,
    const typename Gt::Site_2 & q, const typename Gt::Site_2 & r,
    const typename Gt::Site_2 & t)
{
  typedef typename Gt::Vertex_conflict_2 Vertex_conflict_2;
  Gt gt;
  Vertex_conflict_2 incircle = gt.vertex_conflict_2_object();
  CGAL::Timer timer;
  timer.start();
  for (size_t i = 0; i < repetitions; ++i) {
    (void) incircle(p, q, r, t);
  }
  timer.stop();
  std::cerr << "Test time = " << timer.time() << "s\n";
}

int main(int argc, const char *argv[])
{
  try {

  TCLAP::CmdLine cmd("Run an incircle test multiple times",
                     ' ', "0.1");
  TCLAP::SwitchArg linfSwitch("", "linf", "Linf incircle test (default)",
      cmd, false);
  TCLAP::SwitchArg l2Switch("", "l2", "L2 incircle test", cmd, false);
  TCLAP::ValueArg<size_t> repsArg("r", "reps", "Number of repetitions",
      false, 1000000, "positive integer", cmd);
  TCLAP::UnlabeledMultiArg<std::string> sitesArg("sites", "Sites",
      true, "sites", cmd);
  cmd.parse(argc, argv);

  const bool is_l2 = l2Switch.getValue();
  if (is_l2 and (linfSwitch.getValue() == true)) {
    std::cerr << "error: both L2 and Linf specified" << std::endl;
    return -2;
  }
  const bool is_linf = not is_l2;
  const size_t reps = repsArg.getValue();
  std::vector<std::string> svec = sitesArg.getValue();

  if (svec.size() != 4) {
    std::cerr << "error: 4 sites must be specified" << std::endl;
    return -3;
  }
  std::istringstream strp(svec[0]);
  std::istringstream strq(svec[1]);
  std::istringstream strr(svec[2]);
  std::istringstream strt(svec[3]);

  if (is_linf) {
    GtLinf::Site_2 p, q, r, t;
    strp >> p;
    strq >> q;
    strr >> r;
    strt >> t;
    run_benchmark<GtLinf>(reps, p, q, r, t);
  } else {
    GtL2::Site_2 p, q, r, t;
    strp >> p;
    strq >> q;
    strr >> r;
    strt >> t;
    run_benchmark<GtL2>(reps, p, q, r, t);
  }

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error()
              << " for arg " << e.argId() << std::endl;
    return -1;
  }
  return 0;
}
