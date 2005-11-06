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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#define VDA_TEST_SVD

#include <CGAL/basic.h>

#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"

#include <CGAL/MP_Float.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_traits_2.h>
#include <CGAL/Segment_Voronoi_diagram_Voronoi_traits_2.h>
#include <CGAL/Segment_Voronoi_diagram_adaptation_policies_2.h>

typedef CGAL::MP_Float  NT;
typedef CGAL::Ring_tag  MTag;

typedef CGAL::Simple_cartesian<NT>      K;
typedef CGAL::Simple_cartesian<double>  DK;
typedef
CGAL::Segment_Voronoi_diagram_filtered_traits_without_intersections_2<DK>
Gt;
//CGAL::Segment_Voronoi_diagram_traits_without_intersections_2<K,MTag>  Gt;

#if 1 // definitions for hierarchy

#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>

typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>                 SVD;
#else
typedef CGAL::Segment_Voronoi_diagram_2<Gt>                           SVD;
#endif
typedef CGAL::Segment_Voronoi_diagram_Voronoi_traits_2<SVD>           VT;
typedef CGAL::Identity_policy_2<SVD,VT>                               IP;

typedef CGAL::Segment_Voronoi_diagram_caching_degeneracy_removal_policy_2<SVD>
CDRP;

typedef CGAL::Voronoi_diagram_2<SVD,VT,IP>                            IVDA;
typedef CGAL::Voronoi_diagram_2<SVD,VT,CDRP>                          CDRVDA;

template<class VD>
void run_tests()
{
  typedef typename VD::Delaunay_graph                       DG;
  typedef typename VD::Voronoi_traits::Site_2               Site_2;
  typedef typename VD::Voronoi_traits::Point_2              Point_2;
  typedef Project_site<typename DG::Vertex_handle,Site_2>   Project_site;
  typedef Project_primal<VD,Point_2>                        Project_primal;
  typedef VDA_Tester<VD,Project_site,Project_primal>        Tester;

  Project_site    project_site;
  Project_primal  project_primal;

  Tester test(project_site, project_primal);

  test.reset_timers();

  test("data/empty.cin");
  test("data/singleton.svd.cin");
  test("data/1D.svd.cin");
  test("data/complicated.svd.cin");
  test("data/non-degenerate.svd.cin");
  test("data/data0.svd.cin");
  test("data/data1.svd.cin");
  test("data/data2.svd.cin");
  test("data/data3.svd.cin");
  test("data/data4.svd.cin");
  test("data/data5.svd.cin");
  test("data/data6.svd.cin");
  test("data/data7.svd.cin");
  test("data/data8.svd.cin");
  test("data/data9.svd.cin");
  test("data/data10.svd.cin");
  test("data/data11.svd.cin");
  test("data/degenerate1.svd.cin");
  test("data/degenerate2.svd.cin");

  test.print_times();

  test.print_separators();
  test.print_separators();

  test.reset_timers();

  test("data/singleton.svd.cin", "data/queries2.cin");
  test("data/1D.svd.cin", "data/queries3.cin");
  test("data/data10.svd.cin", "data/queries2.cin");
  test("data/data11.svd.cin", "data/queries3.cin");

  test.print_loc_times();
}

int main()
{
  run_tests<IVDA>();
  print_separator(std::cout);
  run_tests<CDRVDA>();

  return 0;
}
