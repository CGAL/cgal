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

#define VDA_TEST_SDG

#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"

#include <CGAL/MP_Float.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>

typedef CGAL::MP_Float  NT;
typedef CGAL::Integral_domain_without_division_tag  MTag;

typedef CGAL::Simple_cartesian<NT>      K;
typedef CGAL::Simple_cartesian<double>  DK;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<DK> Gt;

//CGAL::Segment_Delaunay_graph_traits_without_intersections_2<K,MTag>  Gt;

#if 1 // definitions for hierarchy

#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>

typedef CGAL::Segment_Delaunay_graph_hierarchy_2<Gt>                 SDG;
#else
typedef CGAL::Segment_Delaunay_graph_2<Gt>                           SDG;
#endif
typedef CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG>         AT;
typedef CGAL::Identity_policy_2<SDG,AT>                               IP;

typedef CGAL::Segment_Delaunay_graph_caching_degeneracy_removal_policy_2<SDG>
CDRP;

typedef CGAL::Voronoi_diagram_2<SDG,AT,IP>                            IVDA;
typedef CGAL::Voronoi_diagram_2<SDG,AT,CDRP>                          CDRVDA;

template<class VD>
void run_tests()
{
  typedef typename VD::Delaunay_graph                       DG;
  typedef typename VD::Adaptation_traits::Site_2            Site_2;
  typedef typename VD::Adaptation_traits::Point_2           Point_2;
  typedef Project_site<typename DG::Vertex_handle,Site_2>   Project_site;
  typedef Project_primal<VD,Point_2>                        Project_primal;
  typedef VDA_Tester<VD,Project_site,Project_primal>        Tester;

  Project_site    project_site;
  Project_primal  project_primal;

  Tester test(project_site, project_primal);

  test.reset_timers();

  test("data/empty.cin");
  test("data/singleton.sdg.cin");
  test("data/1D.sdg.cin");
  test("data/complicated.sdg.cin");
  test("data/non-degenerate.sdg.cin");
  test("data/data0.sdg.cin");
  test("data/data1.sdg.cin");
  test("data/data2.sdg.cin");
  test("data/data3.sdg.cin");
  test("data/data4.sdg.cin");
  test("data/data5.sdg.cin");
  test("data/data6.sdg.cin");
  test("data/data7.sdg.cin");
  test("data/data8.sdg.cin");
  test("data/data9.sdg.cin");
  test("data/data10.sdg.cin");
  test("data/data11.sdg.cin");
  test("data/degenerate1.sdg.cin");
  test("data/degenerate2.sdg.cin");

  test.print_times();

  test.print_separators();
  test.print_separators();

  test.reset_timers();

  test("data/singleton.sdg.cin", "data/queries2.cin");
  test("data/1D.sdg.cin", "data/queries3.cin");
  test("data/data10.sdg.cin", "data/queries2.cin");
  test("data/data11.sdg.cin", "data/queries3.cin");
  test("data/multiple-edges.sdg.cin", "data/queries9.cin", true);
  test("data/multiple-edges2.sdg.cin", "data/queries10.cin", true);

  test.print_loc_times();
}

int main()
{
  run_tests<IVDA>();
  print_separator(std::cout);
  run_tests<CDRVDA>();

  return 0;
}
