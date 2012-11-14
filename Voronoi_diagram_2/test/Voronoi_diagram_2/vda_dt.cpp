// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#define VDA_TEST_DT

#include <CGAL/basic.h>

#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
#if 1 // definitions for hierarchy

#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>

typedef CGAL::Triangulation_vertex_base_2<K>                        VBB;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<VBB>            VB;
typedef CGAL::Triangulation_data_structure_2<VB>                    TDS;
typedef CGAL::Delaunay_triangulation_2<K,TDS>                       DTB;
typedef CGAL::Triangulation_hierarchy_2<DTB>                        DT;
#else
typedef CGAL::Delaunay_triangulation_2<K>                           DT;
#endif

typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>        AT;
typedef CGAL::Identity_policy_2<DT,AT>                              IP;

typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT>
CDRP;

typedef CGAL::Voronoi_diagram_2<DT,AT,IP>                           IVDA;
typedef CGAL::Voronoi_diagram_2<DT,AT,CDRP>                         CDRVDA;


template<class VD>
void run_tests()
{
  typedef typename VD::Delaunay_graph                         DG;
  typedef typename DG::Geom_traits::Point_2                   Point_2;
  typedef Project_point<typename DG::Vertex_handle,Point_2>   Project_point;
  typedef Project_dual<VD,Point_2>                            Project_dual;
  typedef VDA_Tester<VD,Project_point,Project_dual>           Tester;

  Project_point   project_point;
  Project_dual    project_dual;

  Tester test(project_point, project_dual);

  test.reset_timers();

  test("data/empty.cin");
  test("data/singleton.dt.cin");
  test("data/1D.dt.cin");
  test("data/data1.dt.cin");
  test("data/data2.dt.cin");
  test("data/degenerate1.dt.cin");
  test("data/degenerate2.dt.cin");

  test.print_times();

  test.print_separators();
  test.print_separators();

  test.reset_timers();

  test("data/singleton.dt.cin", "data/queries4.cin");
  test("data/1D.dt.cin", "data/queries3.cin");
  test("data/degenerate1.dt.cin", "data/queries4.cin");

  test.print_loc_times();
}


int main()
{
  run_tests<IVDA>();
  print_separator(std::cout);
  run_tests<CDRVDA>();

  return 0;
}
