// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef VDA_TEST_H
#define VDA_TEST_H 1

#include <CGAL/tags.h>
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
  typedef typename VD::Delaunay_graph   DG;

  template<class OutputIt>
  OutputIt read_from_file(const char* fname, OutputIt it) const {
    std::ifstream ifs(fname);
    assert( fname );

    typename Projector::Site_2  s;

    while ( ifs >> s ) {
      *it++ = s;
    }

    ifs.close();
    return it;
  }

  template<class Iterator>
  void compute_dg(DG& dg, Iterator first, Iterator beyond) const
  {
#ifdef BENCHMARKING
    typename DG::All_edges_iterator eit;
    typename DG::All_vertices_iterator vit;
    for (Iterator it = first; it != beyond; ++it) {
      dg.insert(*it);
      for (eit = dg.all_edges_begin(); eit != dg.all_edges_end(); ++eit) {}
      for (vit = dg.all_vertices_begin();
           vit != dg.all_vertices_end(); ++vit) {}
    }

    typename DG::size_type counter = 0;
    for (eit = dg.all_edges_begin(); eit != dg.all_edges_end(); ++eit) {
      counter++;
    }
    std::cout << "# of vertices: " << dg.number_of_vertices() << std::endl;
    std::cout << "# of faces   : " << dg.number_of_faces() << std::endl;
    std::cout << "# of edges   : " << counter << std::endl;
#else
    for (Iterator it = first; it != beyond; ++it) {
      dg.insert(*it);
    }
#endif
  }

  template<class Iterator>
  void compute_vd(VD& vd, Iterator first, Iterator beyond) const
  {
#ifdef BENCHMARKING
    typename VD::Halfedge_iterator hit;
    typename VD::Face_iterator fit;
    for (Iterator it = first; it != beyond; ++it) {
      vd.insert(*it);
      for (hit = vd.halfedges_begin(); hit != vd.halfedges_end(); ++hit) {}
      for (fit = vd.faces_begin(); fit != vd.faces_end(); ++fit) {}
    }

    assert( vd.is_valid() );

    std::cout << "# of vertices: " << vd.number_of_vertices() << std::endl;
    std::cout << "# of faces   : " << vd.number_of_faces() << std::endl;
    std::cout << "# of h/edges : " << vd.number_of_halfedges() << std::endl;
#else
    for (Iterator it = first; it != beyond; ++it) {
      vd.insert(*it);
    }
#endif
  }

  template<class Iterator>
  VD* compute_vd(const DG& , Iterator first, Iterator beyond,
                 CGAL::Tag_true) const
  {
    // insertion is supported
    VD* vd = new VD();
    compute_vd(*vd, first, beyond);
    return vd;
  }

  template<class Iterator>
  VD* compute_vd(const DG& dg, Iterator, Iterator, CGAL::Tag_false) const
  {
    // insertion is not supported
    VD* vd = new VD(dg);
    return vd;
  }

  void test_vd(const VD& vd) const
  {
    test_dual_graph_concept( vd.dual(), vd.adaptation_traits() );
    test_adaptation_traits_concept( vd.dual(), vd.adaptation_traits() );
    test_adaptation_policy_concept( vd.dual(), vd.adaptation_traits(),
                                    vd.adaptation_policy() );

    std::ofstream nos("");

    test_vda(vd);
    print_report(vd, project_, dp_project_, nos);
  }

  void test_loc(const char* /*fname*/, const char* /*qfname*/, const CGAL::Tag_false&, bool) const
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

  void test_loc(const char* fname, const char* qfname, const CGAL::Tag_true&,
                bool print_sites) const
  {
    std::cout << "*** Testing data file (for point location): "
              << fname << std::endl << std::endl;

    std::vector<typename VD::Adaptation_traits::Site_2> vec_s;
    read_from_file(fname, std::back_inserter(vec_s));

    dg_timer_.start();
    DG dg;
    compute_dg(dg, vec_s.begin(), vec_s.end());
    dg_timer_.stop();

    VD vd(dg);

    loc_timer_.start();
    std::ifstream qfs(qfname);
    assert( qfname );
    test_locate(vd, project_, qfs, std::cout, print_sites);
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

  void operator()(const char* fname) const
  {
    std::cout << "*** Testing data file: " << fname << std::endl
              << std::endl;

    std::vector<typename VD::Adaptation_traits::Site_2> vec_s;
    read_from_file(fname, std::back_inserter(vec_s));

    dg_timer_.start();
    DG dg;
    compute_dg(dg, vec_s.begin(), vec_s.end());
    dg_timer_.stop();

#ifdef BENCHMARKING
    std::cout << std::endl;

    vda_timer_.start();
    VD* vd = compute_vd(dg, vec_s.begin(), vec_s.end(),
                        typename VD::Adaptation_traits::Has_insert());
    vda_timer_.stop();

    std::cout << std::endl << std::endl;

    bool b_dg(false), b_vd(false);
    std::cout << "Is Delaunay graph valid? " << std::flush;
    std::cout << ((b_dg = dg.is_valid()) ? "yes" : "no") << std::endl;

    std::cout << "Is Voronoi diagram valid? " << std::flush;
    std::cout << ((b_vd = vd->is_valid()) ? "yes" : "no") << std::endl;

    assert( b_dg );
    assert( b_vd );
#else
    vda_timer_.start();

    VD* vd = compute_vd(dg, vec_s.begin(), vec_s.end(),
                        typename VD::Adaptation_policy::Has_site_inserter());
    vda_timer_.stop();

    test_vd(*vd);
#endif

    delete vd;

    print_separators();
  }

  void operator()(const char* fname, const char* qfname, bool print_sites = false) const
  {
    typename VD::Adaptation_traits::Has_nearest_site_2 has_ns;
    test_loc(fname, qfname, has_ns, print_sites);
  }

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
