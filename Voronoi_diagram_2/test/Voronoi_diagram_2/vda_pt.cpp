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

#define VDA_TEST_PT

#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"

#include <CGAL/Simple_cartesian.h>
#include "Delaunay_graph_concept.h"
#include "Adaptation_traits_concept.h"
#include "Adaptation_policy_concept.h"


typedef CGAL::Simple_cartesian<double>              K;
typedef CGAL::Delaunay_graph_concept<K>             DG;
typedef CGAL::Adaptation_traits_concept<DG>         AT;
typedef CGAL::Adaptation_policy_concept<DG,AT>      AP;
typedef CGAL::Voronoi_diagram_2<DG,AT,AP>           VDA;

int main()
{
  typedef Project_site<DG::Vertex_handle,DG::Site_2>   Project_site;
  typedef Project_dual<VDA,DG::Site_2>                 Project_dual;
  typedef VDA_Tester<VDA,Project_site,Project_dual>    Tester;

  Project_site   project_site;
  Project_dual   project_dual;

  Tester test(project_site, project_dual);

  test.reset_timers();

  test("data/empty.cin");
  test("data/data1.pt.cin");

  test.print_times();

  test.print_separators();
  test.print_separators();

  test.reset_timers();

  test("data/data1.pt.cin", "data/queries1.cin");
  test("data/data1.pt.cin", "data/queries1.cin");

  test.print_loc_times();

  return 0;
}
