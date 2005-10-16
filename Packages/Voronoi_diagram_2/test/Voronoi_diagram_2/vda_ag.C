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

#define VDA_TEST_AG

#include <CGAL/basic.h>

#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>
#include <CGAL/Apollonius_graph_Voronoi_traits_2.h>

typedef CGAL::Simple_cartesian<double>                        Rep;
typedef CGAL::Apollonius_graph_filtered_traits_2<Rep>         K;

#if 1 // definitions for hierarchy

#include <CGAL/Apollonius_graph_hierarchy_2.h>

typedef CGAL::Apollonius_graph_hierarchy_2<K>                 AG;
#else
typedef CGAL::Apollonius_graph_2<K>                           AG;
#endif

typedef CGAL::Apollonius_graph_Voronoi_traits_2<AG>           UVT;
typedef CGAL::Apollonius_graph_caching_Voronoi_traits_2<AG>   CVT;
typedef CGAL::Apollonius_graph_identity_Voronoi_traits_2<AG>  IVT;
typedef CGAL::Voronoi_diagram_2<AG,CVT>                       VDA;

int main()
{
  typedef Project_site<AG::Vertex_handle,AG::Site_2>    Project_site;
  typedef Project_ag_dual<VDA,AG::Site_2>               Project_ag_dual;
  typedef VDA_Tester<VDA,Project_site,Project_ag_dual>  Tester;

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

  test.print_times();

  test.print_separators();
  test.print_separators();

  test.reset_timers();

  test("data/singleton.ag.cin", "data/queries1.cin");
  test("data/1D.ag.cin", "data/queries1.cin");
  test("data/data5.ag.cin", "data/queries1.cin");

  test.print_loc_times();

  return 0;
}
