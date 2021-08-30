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

#define VDA_TEST_RT


#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"
#include "vda_print_report.h"


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_policies_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      Gt;

#if 1 // definitions for hierarchy

#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>

typedef CGAL::Regular_triangulation_vertex_base_2<Gt>              VBB;
typedef CGAL::Regular_triangulation_face_base_2<Gt>                FB;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<VBB>           VB;
typedef CGAL::Triangulation_data_structure_2<VB,FB>                TDS;
typedef CGAL::Regular_triangulation_2<Gt,TDS>                      RTB;
typedef CGAL::Triangulation_hierarchy_2<RTB>                       RT;
#else
typedef CGAL::Regular_triangulation_2<Gt>                          RT;
#endif

typedef CGAL::Regular_triangulation_adaptation_traits_2<RT>        AT;
typedef CGAL::Identity_policy_2<RT,AT>                             IP;
typedef CGAL::Regular_triangulation_caching_degeneracy_removal_policy_2<RT>
CDRP;

typedef CGAL::Voronoi_diagram_2<RT,AT,IP>                          IVDA;
typedef CGAL::Voronoi_diagram_2<RT,AT,CDRP>                        CDRVDA;


template<class VD>
void run_tests()
{
  typedef typename VD::Delaunay_graph                          DG;
  typedef typename DG::Geom_traits::Weighted_point_2           Site_2;
  typedef Project_point<typename RT::Vertex_handle,Site_2>     Project_point;
  typedef Project_dual<VD,typename RT::Geom_traits::Point_2>   Project_dual;
  typedef VDA_Tester<VD,Project_point,Project_dual>            Tester;

  Project_point   project_point;
  Project_dual    project_dual;

  Tester test(project_point, project_dual);

  test.reset_timers();

  test("data/empty.cin");
  test("data/singleton.rt.cin");
  test("data/1D.rt.cin");
  test("data/data1.rt.cin");
  test("data/data2.rt.cin");
  test("data/data3.rt.cin");
  test("data/degenerate.rt.cin");

  test.print_times();

  test.print_separators();
  test.print_separators();

  test.reset_timers();

  test("data/1D.rt.cin", "data/queries6.cin");
  test("data/data3.rt.cin", "data/queries5.cin");

  test.print_loc_times();
}

int main()
{
  run_tests<IVDA>();
  print_separator(std::cout);
  run_tests<CDRVDA>();

  return 0;
}
