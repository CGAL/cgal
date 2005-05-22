#ifndef VDA_TEST_H
#define VDA_TEST_H 1

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include "vda_print_report.h"
#include "vda_test_concept.h"

template<class VD, class Projector, class Dual_primal_projector>
class VDA_Tester
{
 private:
  typedef typename VD::Dual_graph   DG;

  void compute_dg(char* fname, DG& dg) const
  {
    std::ifstream ifs(fname);
    assert( fname );

    typename Projector::Site_2  s;

    while ( ifs >> s ) { dg.insert(s); }

    ifs.close();
  }

  void test_vd(const DG& dg) const
  {
    
    VD vd(dg);

    test_dual_graph_concept( vd.dual() );
    test_voronoi_traits_concept( vd.voronoi_traits() );

    std::ofstream nos("");

    print_report(vd, project_, dp_project_, nos);
  }


  void print_separators() const
  {
    std::cout << std::endl << std::endl;
    std::cout << "================================="
	      << "=================================" << std::endl;
    std::cout << "================================="
	      << "=================================" << std::endl;
    std::cout << std::endl << std::endl;
  }

 public:
  VDA_Tester(const Projector& project,
	     const Dual_primal_projector& dp_project)
    : project_(project), dp_project_(dp_project) {}

  void operator()(char* fname) const
  {
    std::cout << "*** Testing data file: " << fname << std::endl
	      << std::endl;
    dg_timer_.start();
    DG dg;
    compute_dg(fname, dg);
    dg_timer_.stop();

    vda_timer_.start();
    test_vd(dg);
    vda_timer_.stop();
    print_separators();
  }

  template<class P>
  void operator()(char* fname, const P& p) const
  {
    std::cout << "*** Testing data file: " << fname << std::endl
	      << std::endl;

    // if the predicate p is true on the Delaunay graph we do not
    // test further...
    dg_timer_.start();
    DG dg;
    compute_dg(fname, dg);
    dg_timer_.stop();

    vda_timer_.start();
    if ( !p(dg) ) {
      test_vd(dg);
    }

    vda_timer_.stop();
    print_separators();
  }

  double dg_time() const { return dg_timer_.time(); }
  double vda_time() const { return vda_timer_.time(); }
  double total_time() const { return dg_time() + vda_time(); }
  void reset_timers() const {
    dg_timer_.reset();
    vda_timer_.reset();
  }

  void print_times() const {
    std::cerr << "Elapsed time for the Delaunay graph (sec): "
	      << dg_time() << std::endl << std::endl;
    std::cerr << "Elapsed time for the Voronoi diagram adaptor (sec): "
	      << vda_time() << std::endl << std::endl;
    std::cerr << "Total elapsed time (sec): " << total_time()
	      << std::endl << std::endl;
  }

 private:
  Projector              project_;
  Dual_primal_projector  dp_project_;
  mutable CGAL::Timer    dg_timer_;
  mutable CGAL::Timer    vda_timer_;
};

#endif // VDA_TEST_H
