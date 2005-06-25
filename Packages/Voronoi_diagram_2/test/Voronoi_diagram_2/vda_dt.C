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

#include <CGAL/basic.h>

#include <CGAL/Voronoi_diagram_adaptor_2.h>
#include "vda_test.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_Voronoi_traits_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
#if 1 // definitions for hierarchy

#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>

typedef CGAL::Triangulation_vertex_base_2<K>                      VBB;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<VBB>          VB;
typedef CGAL::Triangulation_data_structure_2<VB>                  TDS;
typedef CGAL::Delaunay_triangulation_2<K,TDS>                     DTB;
typedef CGAL::Triangulation_hierarchy_2<DTB>                      DT;
#else
typedef CGAL::Delaunay_triangulation_2<K>                         DT;
#endif
//typedef CGAL::Delaunay_triangulation_cached_Voronoi_traits_2<DT>  VT;
typedef CGAL::Delaunay_triangulation_Voronoi_traits_2<DT>         VT;
typedef CGAL::Voronoi_diagram_adaptor_2<DT,VT>                    VDA;


struct DT_Predicate
{
  bool operator()(const DT& dt) const {
    if ( dt.dimension() == 1 ) {
      std::cerr << "The Delaunay triangulation is 1-dimensional."
		<< std::endl;
      std::cerr << "Cannot view the Delaunay triangulation as an arrangement."
		<< std::endl;
      return true;
    }
    return false;
  }
};

int main()
{
  typedef DT::Geom_traits::Point_2 Point_2;
  typedef Project_point<DT::Vertex_handle,Point_2>    Project_point;
  typedef Project_dual<VDA,Point_2>                   Project_dual;
  typedef VDA_Tester<VDA,Project_point,Project_dual>  Tester;

  Project_point   project_point;
  Project_dual    project_dual;

  Tester test(project_point, project_dual);

  //  DT_Predicate p;

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

  return 0;
}
