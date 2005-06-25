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

#include <CGAL/Simple_cartesian.h>
#include "Delaunay_graph_concept.h"
#include "Voronoi_traits_concept.h"


typedef CGAL::Simple_cartesian<double>              K;
typedef CGAL::Delaunay_graph_concept<K>             DG;
typedef CGAL::Voronoi_traits_concept<DG>            VT;
typedef CGAL::Voronoi_diagram_adaptor_2<DG,VT>      VDA;


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
