#ifndef VDA_TEST_H
#define VDA_TEST_H 1

#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include "vda_print_report.h"

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

    DG dg;
    compute_dg(fname, dg);
    test_vd(dg);
    print_separators();
  }

  template<class P>
  void operator()(char* fname, const P& p) const
  {
    std::cout << "*** Testing data file: " << fname << std::endl
	      << std::endl;

    // if the predicate p is true on the Delaunay graph we do not
    // test further...
    DG dg;
    compute_dg(fname, dg);
    
    if ( !p(dg) ) {
      test_vd(dg);
    }

    print_separators();
  }

 private:
  Projector              project_;
  Dual_primal_projector  dp_project_;
};

#endif // VDA_TEST_H
