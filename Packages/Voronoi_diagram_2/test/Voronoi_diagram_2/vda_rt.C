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
#include "vda_print_report.h"


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_Voronoi_traits_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef CGAL::Regular_triangulation_euclidean_traits_2<K>        Gt;

#if 1 // definitions for hierarchy

#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>

typedef CGAL::Regular_triangulation_vertex_base_2<Gt>            VBB;
typedef CGAL::Regular_triangulation_face_base_2<Gt>              FB;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<VBB>         VB;
typedef CGAL::Triangulation_data_structure_2<VB,FB>              TDS;
typedef CGAL::Regular_triangulation_2<Gt,TDS>                    RTB;
typedef CGAL::Triangulation_hierarchy_2<RTB>                     RT;
#else
typedef CGAL::Regular_triangulation_2<Gt>                        RT;
#endif
typedef CGAL::Regular_triangulation_Voronoi_traits_2<RT>  VT;
//typedef CGAL::Regular_triangulation_cached_Voronoi_traits_2<RT>  VT;
typedef CGAL::Voronoi_diagram_adaptor_2<RT,VT>                   VDA;


struct RT_Predicate
{
  bool operator()(const RT& rt) const {
    if ( rt.dimension() == 1 ) {
      std::cerr << "The regular triangulation is 1-dimensional."
		<< std::endl;
      std::cerr << "Cannot view the regular triangulation as an arrangement."
		<< std::endl;
      return true;
    }
    return false;
  }
};


int main()
{
  typedef RT::Geom_traits::Weighted_point_2           Site_2;
  typedef Project_point<RT::Vertex_handle,Site_2>     Project_point;
  typedef Project_dual<VDA,RT::Geom_traits::Point_2>  Project_dual;
  typedef VDA_Tester<VDA,Project_point,Project_dual>  Tester;

  Project_point   project_point;
  Project_dual    project_dual;

  Tester test(project_point, project_dual);

  //  RT_Predicate p;

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

  return 0;
}
