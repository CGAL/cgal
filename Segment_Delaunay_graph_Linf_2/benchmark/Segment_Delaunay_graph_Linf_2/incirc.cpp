#include <iostream>
#include <sstream>
#include <cassert>

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

template <typename Gt>
void
run_benchmark(const size_t repetitions, const typename Gt::Site_2 & p,
    const typename Gt::Site_2 & q,
    const typename Gt::Site_2 & t)
{
  typedef typename Gt::Vertex_conflict_2 Vertex_conflict_2;
  Gt gt;
  Vertex_conflict_2 incircle = gt.vertex_conflict_2_object();
  CGAL::Timer timer;
  timer.start();
  for (size_t i = 0; i < repetitions; ++i) {
    (void) incircle(p, q, t);
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
    std::cerr << "error: " << svec.size() << " specified" << std::endl;
    return -3;
  }
  size_t count_inf = 0;
  size_t inf_i = 4; // out of bounds value
  for (size_t i = 0; i < 4; ++i) {
    if (svec[i] == "inf") {
      ++count_inf;
      inf_i = i;
    } else {
      if (svec[i].length() < 5) {
        std::cerr << "error: argument " << i << " too short" << std::endl;
        return -4;
      }
      if (not ((svec[i][0] == 'p') or (svec[i][0] == 's'))) {
        std::cerr << "error: argument " << i << " non site" << std::endl;
        return -5;
      }
      GtLinf::Site_2 read_test;
      std::istringstream stream_test(svec[i]);
      stream_test >> read_test;
      if (!stream_test) {
        std::cerr << "error: reading site " << i << std::endl;
        return -11;
      }
    }
  }
  if ((count_inf > 0) and (inf_i == 3)) {
    std::cerr << "error: test site cannot be at infinity" << std::endl;
    return -6;
  }
  if (count_inf > 1) {
    std::cerr << "error: " << count_inf << " infinite sites" << std::endl;
    return -11;
  }

  std::cout << "Running " << (is_linf ? "Linf" : "L2") << " test:  "
    << svec[0] << "  "  << svec[1] << "  " << svec[2] << "   " << svec[3]
    << std::endl;

  std::istringstream strt(svec[3]);
  if (count_inf == 0) {
    std::istringstream strp(svec[0]);
    std::istringstream strq(svec[1]);
    std::istringstream strr(svec[2]);
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
  } else {
    CGAL_assertion(count_inf == 1);
    std::istringstream strp(
        (inf_i == 0) ? svec[1] : (inf_i == 1) ? svec[2] : svec[0]);
    std::istringstream strq(
        (inf_i == 0) ? svec[2] : (inf_i == 1) ? svec[0] : svec[1]);
    if (is_linf) {
      GtLinf::Site_2 p, q, t;
      strp >> p;
      strq >> q;
      strt >> t;
      run_benchmark<GtLinf>(reps, p, q, t);
    } else {
      GtL2::Site_2 p, q, t;
      strp >> p;
      strq >> q;
      strq >> t;
      run_benchmark<GtL2>(reps, p, q, t);
    }
  }

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error()
              << " for arg " << e.argId() << std::endl;
    return -1;
  }
  return 0;
}
