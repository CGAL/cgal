// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Aymeric PELLE (aymeric.pelle@inria.fr)

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_with_weighted_circumcenter_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

#include <cassert>

#include <CGAL/_test_types.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_vertex_base_3<K>                            Vb;
typedef CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3<K>   Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                            Tds;

typedef CGAL::Regular_triangulation_3<K, Tds>             Regular_triangulation_3;

typedef CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3<K>::Rebind_TDS<Tds>::Other Cell_type;

typedef Regular_triangulation_3::Vertex_handle            Vertex_handle;
typedef Regular_triangulation_3::Cell_handle              Cell_handle;

typedef Regular_triangulation_3::Bare_point               Bare_point;
typedef Regular_triangulation_3::Point                    Point;
typedef Regular_triangulation_3::Weighted_point           Weighted_point;

// Explicit instantiation of the whole class.
template class
CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3<
  K, CGAL::Regular_triangulation_cell_base_3<
      K, CGAL::Triangulation_cell_base_3<
           K, CGAL::Triangulation_ds_cell_base_3<Tds> >,
           CGAL::Keep_hidden_points, std::vector<Weighted_point> > >;

void test_default_constructor ()
{
  Cell_type cell;

  assert(cell.vertex(0) == Vertex_handle());
  assert(cell.vertex(1) == Vertex_handle());
  assert(cell.vertex(2) == Vertex_handle());
  assert(cell.vertex(3) == Vertex_handle());
  assert(cell.neighbor(0) == Cell_handle());
  assert(cell.neighbor(1) == Cell_handle());
  assert(cell.neighbor(2) == Cell_handle());
  assert(cell.neighbor(3) == Cell_handle());
}

void test_constructor_1 ()
{
  Tds::Vertex v0, v1, v2, v3;
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);

  Cell_type cell(vh0, vh1, vh2, vh3);

  assert(cell.vertex(0) == vh0 && cell.vertex(0) != Vertex_handle());
  assert(cell.vertex(1) == vh1 && cell.vertex(1) != Vertex_handle());
  assert(cell.vertex(2) == vh2 && cell.vertex(2) != Vertex_handle());
  assert(cell.vertex(3) == vh3 && cell.vertex(3) != Vertex_handle());
  assert(cell.neighbor(0) == Cell_handle());
  assert(cell.neighbor(1) == Cell_handle());
  assert(cell.neighbor(2) == Cell_handle());
  assert(cell.neighbor(3) == Cell_handle());
}

void test_constructor_2 ()
{
  Tds::Vertex v0, v1, v2, v3;
  Tds::Cell c0, c1, c2, c3;
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);
  Cell_handle ch0 = tds.create_cell(c0);
  Cell_handle ch1 = tds.create_cell(c1);
  Cell_handle ch2 = tds.create_cell(c2);
  Cell_handle ch3 = tds.create_cell(c3);

  Cell_type cell(vh0, vh1, vh2, vh3, ch0, ch1, ch2, ch3);

  assert(cell.vertex(0) == vh0 && cell.vertex(0) != Vertex_handle());
  assert(cell.vertex(1) == vh1 && cell.vertex(1) != Vertex_handle());
  assert(cell.vertex(2) == vh2 && cell.vertex(2) != Vertex_handle());
  assert(cell.vertex(3) == vh3 && cell.vertex(3) != Vertex_handle());
  assert(cell.neighbor(0) == ch0 && cell.neighbor(0) != Cell_handle());
  assert(cell.neighbor(1) == ch1 && cell.neighbor(1) != Cell_handle());
  assert(cell.neighbor(2) == ch2 && cell.neighbor(2) != Cell_handle());
  assert(cell.neighbor(3) == ch3 && cell.neighbor(3) != Cell_handle());
}

void test_weighted_circumcenter ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);

  Cell_type cell(vh0, vh1, vh2, vh3);

  const Bare_point& circumcenter = cell.weighted_circumcenter();
  assert(circumcenter == Bare_point(1,1,1));
  const Bare_point& circumcenter_2 = cell.weighted_circumcenter();
  assert(&circumcenter == &circumcenter_2);
}

void test_copy_constructor ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds::Cell c0, c1, c2, c3;
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);
  Cell_handle ch0 = tds.create_cell(c0);
  Cell_handle ch1 = tds.create_cell(c1);
  Cell_handle ch2 = tds.create_cell(c2);
  Cell_handle ch3 = tds.create_cell(c3);

  Cell_type cell(vh0, vh1, vh2, vh3, ch0, ch1, ch2, ch3);
  Cell_type ccell(cell);

  assert(ccell.vertex(0) == vh0 && ccell.vertex(0) != Vertex_handle());
  assert(ccell.vertex(1) == vh1 && ccell.vertex(1) != Vertex_handle());
  assert(ccell.vertex(2) == vh2 && ccell.vertex(2) != Vertex_handle());
  assert(ccell.vertex(3) == vh3 && ccell.vertex(3) != Vertex_handle());
  assert(ccell.neighbor(0) == ch0 && ccell.neighbor(0) != Cell_handle());
  assert(ccell.neighbor(1) == ch1 && ccell.neighbor(1) != Cell_handle());
  assert(ccell.neighbor(2) == ch2 && ccell.neighbor(2) != Cell_handle());
  assert(ccell.neighbor(3) == ch3 && ccell.neighbor(3) != Cell_handle());

  const Bare_point& circumcenter = cell.weighted_circumcenter();
  const Bare_point& ccircumcenter = ccell.weighted_circumcenter();
  assert(ccircumcenter == Bare_point(1,1,1));
  assert(&circumcenter != &ccircumcenter);
}

void test_copy_constructor_2 ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds::Cell c0, c1, c2, c3;
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);
  Cell_handle ch0 = tds.create_cell(c0);
  Cell_handle ch1 = tds.create_cell(c1);
  Cell_handle ch2 = tds.create_cell(c2);
  Cell_handle ch3 = tds.create_cell(c3);

  Cell_type cell(vh0, vh1, vh2, vh3, ch0, ch1, ch2, ch3);
  const Bare_point& circumcenter = cell.weighted_circumcenter();
  Cell_type ccell(cell);

  assert(ccell.vertex(0) == vh0 && ccell.vertex(0) != Vertex_handle());
  assert(ccell.vertex(1) == vh1 && ccell.vertex(1) != Vertex_handle());
  assert(ccell.vertex(2) == vh2 && ccell.vertex(2) != Vertex_handle());
  assert(ccell.vertex(3) == vh3 && ccell.vertex(3) != Vertex_handle());
  assert(ccell.neighbor(0) == ch0 && ccell.neighbor(0) != Cell_handle());
  assert(ccell.neighbor(1) == ch1 && ccell.neighbor(1) != Cell_handle());
  assert(ccell.neighbor(2) == ch2 && ccell.neighbor(2) != Cell_handle());
  assert(ccell.neighbor(3) == ch3 && ccell.neighbor(3) != Cell_handle());

  const Bare_point& ccircumcenter = ccell.weighted_circumcenter();
  assert(ccircumcenter == Bare_point(1,1,1));
  assert(&circumcenter != &ccircumcenter);
}

void test_assignment_operator ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds::Cell c0, c1, c2, c3;
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);
  Cell_handle ch0 = tds.create_cell(c0);
  Cell_handle ch1 = tds.create_cell(c1);
  Cell_handle ch2 = tds.create_cell(c2);
  Cell_handle ch3 = tds.create_cell(c3);

  Cell_type cell(vh0, vh1, vh2, vh3, ch0, ch1, ch2, ch3);
  Cell_type ccell;
  ccell = cell;

  assert(ccell.vertex(0) == vh0 && ccell.vertex(0) != Vertex_handle());
  assert(ccell.vertex(1) == vh1 && ccell.vertex(1) != Vertex_handle());
  assert(ccell.vertex(2) == vh2 && ccell.vertex(2) != Vertex_handle());
  assert(ccell.vertex(3) == vh3 && ccell.vertex(3) != Vertex_handle());
  assert(ccell.neighbor(0) == ch0 && ccell.neighbor(0) != Cell_handle());
  assert(ccell.neighbor(1) == ch1 && ccell.neighbor(1) != Cell_handle());
  assert(ccell.neighbor(2) == ch2 && ccell.neighbor(2) != Cell_handle());
  assert(ccell.neighbor(3) == ch3 && ccell.neighbor(3) != Cell_handle());

  const Bare_point& circumcenter = cell.weighted_circumcenter();
  const Bare_point& ccircumcenter = ccell.weighted_circumcenter();
  assert(ccircumcenter == Bare_point(1,1,1));
  assert(&circumcenter != &ccircumcenter);
}

void test_assignment_operator_2 ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds::Cell c0, c1, c2, c3;
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);
  Cell_handle ch0 = tds.create_cell(c0);
  Cell_handle ch1 = tds.create_cell(c1);
  Cell_handle ch2 = tds.create_cell(c2);
  Cell_handle ch3 = tds.create_cell(c3);

  Cell_type cell(vh0, vh1, vh2, vh3, ch0, ch1, ch2, ch3);
  const Bare_point& circumcenter = cell.weighted_circumcenter();
  Cell_type ccell;
  ccell = cell;

  assert(ccell.vertex(0) == vh0 && ccell.vertex(0) != Vertex_handle());
  assert(ccell.vertex(1) == vh1 && ccell.vertex(1) != Vertex_handle());
  assert(ccell.vertex(2) == vh2 && ccell.vertex(2) != Vertex_handle());
  assert(ccell.vertex(3) == vh3 && ccell.vertex(3) != Vertex_handle());
  assert(ccell.neighbor(0) == ch0 && ccell.neighbor(0) != Cell_handle());
  assert(ccell.neighbor(1) == ch1 && ccell.neighbor(1) != Cell_handle());
  assert(ccell.neighbor(2) == ch2 && ccell.neighbor(2) != Cell_handle());
  assert(ccell.neighbor(3) == ch3 && ccell.neighbor(3) != Cell_handle());

  const Bare_point& ccircumcenter = ccell.weighted_circumcenter();
  assert(ccircumcenter == Bare_point(1,1,1));
  assert(&circumcenter != &ccircumcenter);
}

void test_set_vertex ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);

  Cell_type cell(vh0, vh1, vh2, vh3);

  const Bare_point& circumcenter = cell.weighted_circumcenter();
  assert(circumcenter == Bare_point(1,1,1));

  Vertex_handle vh0_bis = tds.create_vertex(v0);
  cell.set_vertex(0, vh0_bis);
  assert(cell.vertex(0) == vh0_bis && cell.vertex(0) != vh0 && cell.vertex(0) != Vertex_handle());

  const Bare_point& circumcenter_2 = cell.weighted_circumcenter();
  assert(circumcenter_2 == Bare_point(1,1,1));
}

void test_set_vertices ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);

  Cell_type cell(vh0, vh1, vh2, vh3);

  cell.set_vertices();
  assert(cell.vertex(0) == Vertex_handle());
  assert(cell.vertex(1) == Vertex_handle());
  assert(cell.vertex(2) == Vertex_handle());
  assert(cell.vertex(3) == Vertex_handle());
}

void test_set_vertices_with_parameters ()
{
  Weighted_point p0(Bare_point(0,0,0),1);
  Weighted_point p1(Bare_point(2,0,0),1);
  Weighted_point p2(Bare_point(0,2,0),1);
  Weighted_point p3(Bare_point(0,0,2),1);
  Tds::Vertex v0(p0), v1(p1), v2(p2), v3(p3);
  Tds tds;
  Vertex_handle vh0 = tds.create_vertex(v0);
  Vertex_handle vh1 = tds.create_vertex(v1);
  Vertex_handle vh2 = tds.create_vertex(v2);
  Vertex_handle vh3 = tds.create_vertex(v3);

  Cell_type cell(vh0, vh1, vh2, vh3);

  const Bare_point& circumcenter = cell.weighted_circumcenter();
  assert(circumcenter == Bare_point(1,1,1));

  Vertex_handle vh0_bis = tds.create_vertex(v0);
  Vertex_handle vh1_bis = tds.create_vertex(v1);
  Vertex_handle vh2_bis = tds.create_vertex(v3);
  Vertex_handle vh3_bis = tds.create_vertex(v2);
  cell.set_vertices(vh0_bis, vh1_bis, vh2_bis, vh3_bis);
  assert(cell.vertex(0) == vh0_bis && cell.vertex(0) != vh0 && cell.vertex(0) != Vertex_handle());
  assert(cell.vertex(1) == vh1_bis && cell.vertex(1) != vh1 && cell.vertex(1) != Vertex_handle());
  assert(cell.vertex(2) == vh2_bis && cell.vertex(2) != vh2 && cell.vertex(2) != Vertex_handle());
  assert(cell.vertex(3) == vh3_bis && cell.vertex(3) != vh3 && cell.vertex(3) != Vertex_handle());

  const Bare_point& circumcenter_2 = cell.weighted_circumcenter();
  assert(circumcenter_2 == Bare_point(1,1,1));
}

int main()
{
  test_default_constructor();
  test_constructor_1();
  test_constructor_2();
  test_weighted_circumcenter();
  test_copy_constructor();
  test_copy_constructor_2();
  test_assignment_operator();
  test_assignment_operator_2();
  test_set_vertex();
  test_set_vertices();
  test_set_vertices_with_parameters();

  std::cout << "EXIT_SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
