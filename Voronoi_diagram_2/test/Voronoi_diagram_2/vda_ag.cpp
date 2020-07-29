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

#define VDA_TEST_AG

#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>
#include <CGAL/Apollonius_graph_adaptation_traits_2.h>
#include <CGAL/Apollonius_graph_adaptation_policies_2.h>

typedef CGAL::Simple_cartesian<double> Rep;
typedef CGAL::Apollonius_graph_filtered_traits_2<Rep> K;

#if 1 // definitions for hierarchy

#include <CGAL/Apollonius_graph_hierarchy_2.h>

typedef CGAL::Apollonius_graph_hierarchy_2<K>                  AG;
#else
typedef CGAL::Apollonius_graph_2<K>                            AG;
#endif

typedef CGAL::Apollonius_graph_adaptation_traits_2<AG>         AT;
typedef CGAL::Identity_policy_2<AG,AT>                         IP;

typedef CGAL::Apollonius_graph_caching_degeneracy_removal_policy_2<AG> CDRP;

typedef CGAL::Voronoi_diagram_2<AG,AT,IP>                      IVDA;
typedef CGAL::Voronoi_diagram_2<AG,AT,CDRP>                    CDRVDA;

template<class VD>
void run_tests()
{
  typedef typename VD::Delaunay_graph                        DG;
  typedef typename DG::Site_2                                Site_2;
  typedef Project_site<typename DG::Vertex_handle,Site_2>    Project_site;
  typedef Project_ag_dual<VD,Site_2>                         Project_ag_dual;
  typedef VDA_Tester<VD,Project_site,Project_ag_dual>        Tester;

  Project_site      project_site;
  Project_ag_dual   project_ag_dual;

  Tester test(project_site, project_ag_dual);

  test.reset_timers();

  test("data/empty.cin");
  test("data/singleton.ag.cin");
  test("data/1D.ag.cin");
  test("data/data1.ag.cin");
  test("data/data2.ag.cin");
  test("data/data3.ag.cin");
  test("data/data4.ag.cin");
  test("data/data5.ag.cin");
  test("data/degenerate.ag.cin");
  test("data/multiple-edges.ag.cin");

  test.print_times();

  test.print_separators();
  test.print_separators();

  test.reset_timers();

  test("data/singleton.ag.cin", "data/queries1.cin");
  test("data/1D.ag.cin", "data/queries1.cin");
  test("data/data5.ag.cin", "data/queries1.cin");
  test("data/multiple-edges.ag.cin", "data/queries7.cin", true);
  test("data/data6.ag.cin", "data/queries8.cin", true);
  test("data/data7.ag.cin", "data/queries8.cin", true);
  test("data/data8.ag.cin", "data/queries8.cin", true);

  test.print_loc_times();
}

int main()
{
  run_tests<IVDA>();
  print_separator(std::cout);
  run_tests<CDRVDA>();

  return 0;
}
