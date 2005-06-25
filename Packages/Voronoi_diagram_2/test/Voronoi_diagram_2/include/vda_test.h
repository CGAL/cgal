// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef VDA_TEST_H
#define VDA_TEST_H 1

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include "vda_print_report.h"
#include "vda_test_vda.h"
#include "vda_test_concept.h"
#include "vda_test_locate.h"

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
    test_voronoi_traits_concept( vd.dual(), vd.voronoi_traits() );

    std::ofstream nos("");

    test_vda(vd);
    print_report(vd, project_, dp_project_, nos);
  }

  void test_loc(char* fname, char* qfname, const CGAL::Tag_false&) const
  {
    static int i = 0;
    if ( i == 0 ) {
      i++;
      std::cout << std::endl;
      std::cout << "  ***********************************" << std::endl;
      std::cout << "  * Point location is not supported *" << std::endl;
      std::cout << "  ***********************************" << std::endl;
      std::cout << std::endl;
    }

  }

  void test_loc(char* fname, char* qfname, const CGAL::Tag_true&) const
  {
    std::cout << "*** Testing data file (for point location): "
	      << fname << std::endl << std::endl;

    dg_timer_.start();
    DG dg;
    compute_dg(fname, dg);
    dg_timer_.stop();

    VD vd(dg);

    loc_timer_.start();
    std::ifstream qfs(qfname);
    test_locate(vd, project_, qfs, std::cout);
    loc_timer_.stop();
    print_separators();
  }

 public:
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

  void operator()(char* fname, char* qfname) const
  {
    typename VD::Voronoi_traits::Has_point_locator has_pl;
    test_loc(fname, qfname, has_pl);
  }

#if 0
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
#endif

  double dg_time() const { return dg_timer_.time(); }
  double loc_time() const { return loc_timer_.time(); }
  double vda_time() const { return vda_timer_.time(); }
  double total_time() const { return dg_time() + vda_time(); }
  void reset_timers() const {
    dg_timer_.reset();
    vda_timer_.reset();
    loc_timer_.reset();
  }

  void print_times() const {
    std::cerr << "Elapsed time for the Delaunay graph (sec): "
	      << dg_time() << std::endl << std::endl;
    std::cerr << "Elapsed time for the Voronoi diagram adaptor (sec): "
	      << vda_time() << std::endl << std::endl;
    std::cerr << "Total elapsed time (sec): " << total_time()
	      << std::endl << std::endl;
  }

  void print_loc_times() const {
    std::cerr << "Elapsed time for the Delaunay graph (sec): "
	      << dg_time() << std::endl << std::endl;
    std::cerr << "Elapsed time for point location (sec): "
	      << loc_time() << std::endl << std::endl;
    std::cerr << "Total elapsed time (sec): " << dg_time() + loc_time()
	      << std::endl << std::endl;
  }

 private:
  Projector              project_;
  Dual_primal_projector  dp_project_;
  mutable CGAL::Timer    dg_timer_;
  mutable CGAL::Timer    vda_timer_;
  mutable CGAL::Timer    loc_timer_;
};

#endif // VDA_TEST_H
