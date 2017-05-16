// Copyright (c) 1998, 2015  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Francois Rebufat
//                 Manuel Caroli
//                 Aymeric Pelle

#include "_test_cls_periodic_3_iterator.h"
#include "_test_cls_periodic_3_circulator.h"

#include <cassert>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>

template <class PeriodicTriangulation>
void
_test_periodic_3_triangulation_3_constructors(const PeriodicTriangulation &)
{
  std::cout<<"Creation"<<std::endl;
  PeriodicTriangulation PT_def;
  assert(PT_def.is_valid());

  PeriodicTriangulation PT_dom(
        typename PeriodicTriangulation::Iso_cuboid(-1,-2,-3,3,2,1));
  assert(PT_def.is_valid());

  PeriodicTriangulation PT_cp(PT_dom);
  assert(PT_dom.is_valid());
  assert(PT_cp == PT_dom);
}

template <class PeriodicTriangulation>
void
_test_cls_periodic_3_triangulation_3(const PeriodicTriangulation &,
                                     const typename PeriodicTriangulation::Point& pointtt,
                                     const char* covering_test_HOM_filename,
                                     const char* covering_test_filename,
                                     bool ex = false,
                                     bool hom = false,
                                     bool test_input_ouput = true)
{
  typedef PeriodicTriangulation                  P3T3;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename P3T3::Geometric_traits        Geometric_traits;
  typedef typename P3T3::Triangulation_data_structure
                                                 Triangulation_data_structure;

  typedef typename P3T3::Offset                  Offset;
  typedef typename P3T3::Iso_cuboid              Iso_cuboid;
  typedef typename P3T3::Covering_sheets         Covering_sheets;
  CGAL_USE_TYPE(Covering_sheets);

  typedef typename P3T3::Point_3                 Point_3;
  typedef typename P3T3::Point                   Point;
  typedef typename P3T3::Segment                 Segment;
  typedef typename P3T3::Triangle                Triangle;
  typedef typename P3T3::Tetrahedron             Tetrahedron;
  CGAL_USE_TYPE(Segment);
  CGAL_USE_TYPE(Triangle);
  CGAL_USE_TYPE(Tetrahedron);

  typedef typename P3T3::Periodic_point          Periodic_point;
  typedef typename P3T3::Periodic_segment        Periodic_segment;
  typedef typename P3T3::Periodic_triangle       Periodic_triangle;
  typedef typename P3T3::Periodic_tetrahedron    Periodic_tetrahedron;
  CGAL_USE_TYPE(Periodic_point);
  CGAL_USE_TYPE(Periodic_segment);
  CGAL_USE_TYPE(Periodic_triangle);
  CGAL_USE_TYPE(Periodic_tetrahedron);

  typedef typename P3T3::Vertex                  Vertex;
  typedef typename P3T3::Cell                    Cell;
  typedef typename P3T3::Edge                    Edge;
  typedef typename P3T3::Facet                   Facet;
  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Cell);

  typedef typename P3T3::Vertex_handle           Vertex_handle;
  typedef typename P3T3::Cell_handle             Cell_handle;

  typedef typename P3T3::size_type               size_type;
  typedef typename P3T3::difference_type         difference_type;
  CGAL_USE_TYPE(size_type);
  CGAL_USE_TYPE(difference_type);

  typedef typename P3T3::Cell_iterator           Cell_iterator;
  typedef typename P3T3::Facet_iterator          Facet_iterator;
  typedef typename P3T3::Edge_iterator           Edge_iterator;
  typedef typename P3T3::Vertex_iterator         Vertex_iterator;
  CGAL_USE_TYPE(Cell_iterator);
  CGAL_USE_TYPE(Facet_iterator);
  CGAL_USE_TYPE(Edge_iterator);
  CGAL_USE_TYPE(Vertex_iterator);

  typedef typename P3T3::Cell_circulator         Cell_circulator;
  typedef typename P3T3::Facet_circulator        Facet_circulator;
  CGAL_USE_TYPE(Cell_circulator);
  CGAL_USE_TYPE(Facet_circulator);

  typedef typename P3T3::Periodic_tetrahedron_iterator
                                                 Periodic_tetrahedron_iterator;
  typedef typename P3T3::Periodic_triangle_iterator
                                                 Periodic_triangle_iterator;
  typedef typename P3T3::Periodic_segment_iterator
                                                 Periodic_segment_iterator;
  typedef typename P3T3::Periodic_point_iterator Periodic_point_iterator;
  CGAL_USE_TYPE(Periodic_tetrahedron_iterator);
  CGAL_USE_TYPE(Periodic_triangle_iterator);
  CGAL_USE_TYPE(Periodic_segment_iterator);
  CGAL_USE_TYPE(Periodic_point_iterator);

  typedef typename P3T3::Locate_type             Locate_type;
  typedef typename P3T3::Iterator_type           Iterator_type;
  CGAL_USE_TYPE(Iterator_type);

  std::cout<<"Creating triangulations to test on"<<std::endl;
  std::cout<<"  non-degenerate, 3-sheeted covering"<<std::endl;
  P3T3 PT3(Iso_cuboid(-1,-1,-1,1,1,1));
  PT3.insert(Point(-0.578898, -0.421525,   0.841529 ));
  PT3.insert(Point(-0.341375,  0.0654286,  0.458741 ));
  PT3.insert(Point( 0.0857005, 0.51015,    0.504106 ));
  PT3.insert(Point( 0.633515,  0.271829,   0.369029 ));
  PT3.insert(Point( 0.415803, -0.251537,  -0.51599  ));
  PT3.insert(Point( 0.367819, -0.892225,   0.396091 ));
  PT3.insert(Point( 0.614818,  0.844804,  -0.0642087));
  PT3.insert(Point( 0.846336,  0.0261109, -0.754327 ));
  PT3.insert(Point( 0.154865,  0.287508,   0.377003 ));
  PT3.insert(Point(-0.255508,  0.816668,   0.823991 ));
  assert(PT3.number_of_sheets() == CGAL::make_array(3,3,3));
  assert(PT3.number_of_vertices() == 10);
  assert(PT3.is_valid());

  std::cout<<"  non-degenerate, 1-sheeted covering"<<std::endl;
  std::ifstream fin;
  if (hom)
    fin.open(covering_test_HOM_filename);
  else
    fin.open(covering_test_filename);
  P3T3 PT1;
  fin >> PT1;
//  assert(PT1.number_of_vertices() == 70);
  assert(PT1.is_valid());
//  assert(!PT1.is_extensible_triangulation_in_1_sheet_h1());
//  assert(!PT1.is_extensible_triangulation_in_1_sheet_h2());
  PT1.insert(pointtt);
  assert(PT1.is_extensible_triangulation_in_1_sheet_h1());
  assert(PT1.is_extensible_triangulation_in_1_sheet_h2());
  assert(PT1.number_of_sheets() == CGAL::make_array(1,1,1));
//  assert(PT1.number_of_vertices() == 71);
  assert(PT1.is_valid());

  std::cout<<"  degenerate, 3-sheeted covering"<<std::endl;
  P3T3 PT3_deg(Iso_cuboid(0,0,0,6,6,6));
  for (int i=0 ; i<6 ; i+=2)
    for (int j=0 ; j<6 ; j+=2)
      for (int k=0 ; k<6 ; k+=2)
        PT3_deg.insert(Point(i,j,k));

  assert(PT3_deg.number_of_sheets() == CGAL::make_array(3,3,3));
  assert(PT3_deg.number_of_vertices() == 27);
  assert(PT3_deg.is_valid());

  std::cout<<"  degenerate, 1-sheeted covering"<<std::endl;
  P3T3 PT1_deg(Iso_cuboid(0,0,0,4,4,4));;
  for (unsigned i = 0; i < 8; ++i)
    for (unsigned j = 0; j < 8; ++j)
      for (unsigned k = 0; k < 8; ++k)
        PT1_deg.insert(Point(static_cast<float>(i)*4./8.,static_cast<float>(j)*4./8.,static_cast<float>(k)*4./8.));

  assert(PT1_deg.number_of_sheets() == CGAL::make_array(1,1,1));
//  assert(PT1_deg.number_of_vertices() == 54);
  assert(PT1_deg.is_valid());

  std::cout<<"Constructor"<<std::endl;

  // RT used to make it work with homogeneous kernels
  typename Geometric_traits::RT ft1(1);
  typename Geometric_traits::RT ft2(2);
  Iso_cuboid domain(ft1, ft1, ft1, ft2, ft2, ft2, 10);
  std::cout << "x-length: "<<domain.xmax()-domain.xmin()<<'\t'
            << "y-length: "<<domain.ymax()-domain.ymin()<<'\t'
            << "z-length: "<<domain.zmax()-domain.zmin()<<std::endl;

  std::cout<<"Comparisons: "
           << (domain.xmax()-domain.xmin() == domain.ymax()-domain.ymin())
           << (domain.xmax()-domain.xmin() == domain.zmax()-domain.zmin())
           << (domain.ymax()-domain.ymin() == domain.zmax()-domain.zmin())
           << ((domain.xmax()-domain.xmin()) == (domain.ymax()-domain.ymin()))
           << ((domain.xmax()-domain.xmin()) == (domain.zmax()-domain.zmin()))
           << ((domain.ymax()-domain.ymin()) == (domain.zmax()-domain.zmin()))
           << std::endl;

  P3T3 PT_constr(domain);

  std::cout<<"Assignment"<<std::endl;

  P3T3 PT;
  PT = PT3;
  assert(PT == PT3);
  PT = PT1;
  assert(PT == PT1);
  PT = PT3_deg;
  assert(PT == PT3_deg);
  PT = PT1_deg;
  assert(PT == PT1_deg);
  assert(PT.is_valid());

  PT.swap(PT3);
  assert(PT3 == PT1_deg);
  assert(PT != PT3);
  assert(PT.is_valid());
  assert(PT3.is_valid());

  PT3.swap(PT);
  assert(PT3 != PT1_deg);
  assert(PT == PT1_deg);
  assert(PT.is_valid());
  assert(PT3.is_valid());

  PT.clear();
  assert(PT.number_of_vertices() == 0);
  assert(PT.number_of_sheets() == CGAL::make_array(3,3,3));
  assert(PT.is_valid());

  // operator== and operator!= are tested in the asserts above

  std::cout<<"Access"<<std::endl;

  Geometric_traits gt = PT.geom_traits();
  Triangulation_data_structure tds = PT.tds();
  assert(P3T3().domain() == Iso_cuboid(0,0,0,1,1,1));
  assert(PT1.domain() == Iso_cuboid(-1,-1,-1,1,1,1));

  P3T3 PT_change_domain(PT1);
  PT_change_domain.set_domain(Iso_cuboid(1,2,3,4,5,6));
  assert(PT_change_domain.number_of_vertices() == 0);
  assert(PT_change_domain.domain() == Iso_cuboid(1,2,3,4,5,6));

  assert(PT3.number_of_sheets() == CGAL::make_array(3,3,3));
  assert(PT1.number_of_sheets() == CGAL::make_array(1,1,1));
  assert(PT3_deg.number_of_sheets() == CGAL::make_array(3,3,3));
  assert(PT1_deg.number_of_sheets() == CGAL::make_array(1,1,1));

  assert(PT1.is_extensible_triangulation_in_1_sheet_h1());
  assert(!PT3.is_extensible_triangulation_in_1_sheet_h1());
  assert(PT1_deg.is_extensible_triangulation_in_1_sheet_h1());
  assert(!PT3_deg.is_extensible_triangulation_in_1_sheet_h1());

  assert(PT1.is_extensible_triangulation_in_1_sheet_h2());
  assert(!PT3.is_extensible_triangulation_in_1_sheet_h2());
  assert(PT1_deg.is_extensible_triangulation_in_1_sheet_h2());
  assert(!PT3_deg.is_extensible_triangulation_in_1_sheet_h2());

  assert(PT1.is_triangulation_in_1_sheet());
  assert(!PT3.is_triangulation_in_1_sheet());
  assert(PT1_deg.is_triangulation_in_1_sheet());
  assert(PT3_deg.is_triangulation_in_1_sheet());

  PT3_deg.convert_to_1_sheeted_covering();
  assert(PT3_deg.number_of_sheets() == CGAL::make_array(1,1,1));
  assert(PT3_deg.is_valid());
  PT3_deg.convert_to_1_sheeted_covering();
  assert(PT3_deg.number_of_sheets() == CGAL::make_array(1,1,1));
  assert(PT3_deg.is_valid());
  PT3_deg.convert_to_27_sheeted_covering();
  assert(PT3_deg.number_of_sheets() == CGAL::make_array(3,3,3));
  assert(PT3_deg.is_valid());
  PT3_deg.convert_to_27_sheeted_covering();
  assert(PT3_deg.number_of_sheets() == CGAL::make_array(3,3,3));
  assert(PT3_deg.is_valid());

  // Check the Euler relation
  assert(PT1.number_of_vertices() - PT1.number_of_edges()
         + PT1.number_of_facets() - PT1.number_of_cells() == 0);
  assert(PT1.number_of_stored_vertices() - PT1.number_of_stored_edges()
         + PT1.number_of_stored_facets() - PT1.number_of_stored_cells() == 0);
  assert(PT3.number_of_vertices() - PT3.number_of_edges()
         + PT3.number_of_facets() - PT3.number_of_cells() == 0);
  assert(PT3.number_of_stored_vertices() - PT3.number_of_stored_edges()
         + PT3.number_of_stored_facets() - PT3.number_of_stored_cells() == 0);
  assert(PT1_deg.number_of_vertices() - PT1_deg.number_of_edges()
         + PT1_deg.number_of_facets() - PT1_deg.number_of_cells() == 0);
  assert(PT1_deg.number_of_stored_vertices() - PT1_deg.number_of_stored_edges()
         + PT1_deg.number_of_stored_facets()-PT1_deg.number_of_stored_cells()==0);
  assert(PT3_deg.number_of_vertices() - PT3_deg.number_of_edges()
         + PT3_deg.number_of_facets() - PT3_deg.number_of_cells() == 0);
  assert(PT3_deg.number_of_stored_vertices() - PT3_deg.number_of_stored_edges()
         + PT3_deg.number_of_stored_facets()-PT3_deg.number_of_stored_cells()==0);

  // Checking the number copies of each item in 27-sheeted covering space
  assert(PT1.number_of_cells()    == PT1.number_of_stored_cells());
  assert(PT1.number_of_facets()   == PT1.number_of_stored_facets());
  assert(PT1.number_of_edges()    == PT1.number_of_stored_edges());
  assert(PT1.number_of_vertices() == PT1.number_of_stored_vertices());
  assert(PT3.number_of_cells()    * 27 == PT3.number_of_stored_cells());
  assert(PT3.number_of_facets()   * 27 == PT3.number_of_stored_facets());
  assert(PT3.number_of_edges()    * 27 == PT3.number_of_stored_edges());
  assert(PT3.number_of_vertices() * 27 == PT3.number_of_stored_vertices());
  assert(PT1_deg.number_of_cells()    == PT1_deg.number_of_stored_cells());
  assert(PT1_deg.number_of_facets()   == PT1_deg.number_of_stored_facets());
  assert(PT1_deg.number_of_edges()    == PT1_deg.number_of_stored_edges());
  assert(PT1_deg.number_of_vertices() == PT1_deg.number_of_stored_vertices());
  assert(PT3_deg.number_of_cells()   * 27 == PT3_deg.number_of_stored_cells());
  assert(PT3_deg.number_of_facets()  * 27 == PT3_deg.number_of_stored_facets());
  assert(PT3_deg.number_of_edges()   * 27 == PT3_deg.number_of_stored_edges());
  assert(PT3_deg.number_of_vertices()*27==PT3_deg.number_of_stored_vertices());

  std::cout<<"Geometric access functions"<<std::endl;

  Cell_handle ch = PT3.locate(Point(-1,-1,1));
  assert(PT3.periodic_point(ch->vertex(0)).second != Offset());
  assert(PT3.periodic_point(ch,2).second != Offset());
  PT3.point(PT3.periodic_point(ch,0));

  assert(PT3.periodic_segment(Edge(ch,2,3)).at(0).second != Offset());
  assert(PT3.periodic_segment(ch,2,1).at(1).second != Offset());
  PT3.construct_segment(PT3.periodic_segment(ch,0,3));

  assert(PT3.periodic_triangle(Facet(ch,0)).at(2).second != Offset());
  assert(PT3.periodic_triangle(ch,2).at(1).second != Offset());
  PT3.construct_triangle(PT3.periodic_triangle(ch,0));

  assert(PT3.periodic_tetrahedron(ch).at(3).second != Offset());
  PT3.construct_tetrahedron(PT3.periodic_tetrahedron(ch));

  std::cout<<"Queries"<<std::endl;

  Vertex_handle vh;
  assert(PT3.is_vertex(ch->vertex(0)->point(),vh));
  assert(ch->vertex(0) != vh);

  assert(!PT3.is_vertex(Point(0,0,0),vh));

  assert(PT1.is_vertex(PT1.vertices_begin()->point(),vh));
  assert(PT1.vertices_begin() == vh);

  vh = Vertex_handle();
  assert(!PT3.is_vertex(vh));
  assert(PT3.is_vertex(ch->vertex(2)));

  Cell_handle c;
  int i,j,k,l;
  Offset off0 = PT3.periodic_tetrahedron(ch).at(0).second;
  Offset off1 = PT3.periodic_tetrahedron(ch).at(1).second;
  Offset off2 = PT3.periodic_tetrahedron(ch).at(2).second;
  Offset off3 = PT3.periodic_tetrahedron(ch).at(3).second;

  assert(PT3.is_edge(ch->vertex(0),ch->vertex(1),c,i,j));
  assert(ch->vertex(0) == c->vertex(i));
  assert(ch->vertex(1) == c->vertex(j));
  assert(!PT3.is_edge(ch->vertex(0),
                      ch->neighbor(0)->vertex(ch->neighbor(0)->index(ch)),c,i,j));
  assert(PT3.is_edge(ch->vertex(0),off0,ch->vertex(1),off1,c,i,j));
  assert(ch->vertex(0) == c->vertex(i));
  assert(ch->vertex(1) == c->vertex(j));
  assert(!PT3.is_edge(ch->vertex(0),off0,ch->vertex(1),off1+Offset(0,0,1),
                      c,i,j));

  assert(PT3.is_facet(ch->vertex(0),ch->vertex(1),ch->vertex(2), c,i,j,k));
  assert(ch->vertex(0) == c->vertex(i));
  assert(ch->vertex(1) == c->vertex(j));
  assert(ch->vertex(2) == c->vertex(k));
  assert(!PT3.is_facet(ch->vertex(0),ch->vertex(1),
                       ch->neighbor(0)->vertex(ch->neighbor(0)->index(ch)),c,i,j,k));
  assert(PT3.is_facet(ch->vertex(0),off0,ch->vertex(1),off1,ch->vertex(2),off2,
                      c,i,j,k));
  assert(ch->vertex(0) == c->vertex(i));
  assert(ch->vertex(1) == c->vertex(j));
  assert(ch->vertex(2) == c->vertex(k));
  assert(!PT3.is_facet(ch->vertex(0),off0,ch->vertex(1),off1,
                       ch->vertex(2),off2+Offset(0,0,1),c,i,j,k));

  c = Cell_handle();
  assert(PT3.is_cell(ch));
  assert(!PT3.is_cell(c));
  assert(PT3.is_cell(ch->vertex(0),ch->vertex(1),ch->vertex(2),ch->vertex(3),
                     c));
  assert(!PT3.is_cell(ch->vertex(0),ch->vertex(1),ch->vertex(2),
                      ch->neighbor(0)->vertex(ch->neighbor(0)->index(ch)),c));
  assert(PT3.is_cell(ch->vertex(0),ch->vertex(1),ch->vertex(2),ch->vertex(3),
                     c,i,j,k,l));
  assert(ch->vertex(0) == c->vertex(i));
  assert(ch->vertex(1) == c->vertex(j));
  assert(ch->vertex(2) == c->vertex(k));
  assert(ch->vertex(3) == c->vertex(l));
  assert(!PT3.is_cell(ch->vertex(0),ch->vertex(1),ch->vertex(2),
                      ch->neighbor(0)->vertex(ch->neighbor(0)->index(ch)),c,i,j,k,l));
  assert(PT3.is_cell(ch->vertex(0),off0,ch->vertex(1),off1,
                     ch->vertex(2),off2,ch->vertex(3),off3,c,i,j,k,l));
  assert(PT3.is_cell(ch->vertex(0),off0,ch->vertex(1),off1,
                     ch->vertex(2),off2,ch->vertex(3),off3,c));
  assert(ch->vertex(0) == c->vertex(i));
  assert(ch->vertex(1) == c->vertex(j));
  assert(ch->vertex(2) == c->vertex(k));
  assert(ch->vertex(3) == c->vertex(l));
  assert(!PT3.is_cell(ch->vertex(0),off0,ch->vertex(1),off1+Offset(0,0,1),
                      ch->vertex(2),off2,ch->vertex(3),off3,c,i,j,k,l));
  assert(!PT3.is_cell(ch->vertex(0),off0,ch->vertex(1),off1+Offset(0,0,1),
                      ch->vertex(2),off2,ch->vertex(3),off3,c));

  assert(PT3.has_vertex(Facet(ch,0),ch->vertex(1),i));
  assert(PT3.has_vertex(ch,0,ch->vertex(1),i));
  assert(i==1);
  assert(!PT3.has_vertex(Facet(ch,0),ch->vertex(0),i));
  assert(PT3.has_vertex(Facet(ch,0),ch->vertex(1)));
  assert(PT3.has_vertex(ch,0,ch->vertex(1)));
  assert(!PT3.has_vertex(Facet(ch,0),ch->vertex(0)));

  Cell_handle nb = ch->neighbor(0);
  i = nb->index(ch);
  assert(PT3.are_equal(ch,0,nb,i));
  assert(PT3.are_equal(Facet(ch,0),Facet(nb,i)));
  assert(PT3.are_equal(Facet(ch,0),nb,i));
  assert(!PT3.are_equal(ch,1,nb,i));
  assert(!PT3.are_equal(Facet(ch,1),Facet(nb,i)));
  assert(!PT3.are_equal(Facet(ch,1),nb,i));

  std::cout<<"Point location"<<std::endl;
  ch = PT3_deg.inexact_locate(Point(0.5,0.5,0.5));
  assert( PT3_deg.construct_tetrahedron(PT3_deg.periodic_tetrahedron(ch))
          == PT3_deg.geom_traits().construct_tetrahedron_3_object()(
            Point_3(0,0,0), Point_3(0,2,0), Point_3(0,0,2), Point_3(2,0,0)) );

  ch = PT3_deg.locate(Point(0.5,0.5,0.5));
  assert( PT3_deg.construct_tetrahedron(PT3_deg.periodic_tetrahedron(ch))
          == PT3_deg.geom_traits().construct_tetrahedron_3_object()(
            Point_3(0,0,0), Point_3(0,2,0), Point_3(0,0,2), Point_3(2,0,0)) );

  Locate_type lt;
  int li,lj;
  c = PT3_deg.locate(Point(0.5,0.5,0.5),lt,li,lj);
  assert(c == ch);
  assert(lt == P3T3::CELL);
  assert(PT3_deg.side_of_cell(Point(0.5,0.5,0.5),c,lt,li,lj)
         == CGAL::ON_BOUNDED_SIDE);
  assert(lt == P3T3::CELL);
  assert(PT3_deg.side_of_cell(Point(0.5,0.5,0.5),c->neighbor(0),lt,li,lj)
         == CGAL::ON_UNBOUNDED_SIDE);

  c = PT3_deg.locate(Point(2,0.5,0.5),lt,li,lj);
  assert(lt == P3T3::FACET);
  assert( PT3_deg.construct_triangle(PT3_deg.periodic_triangle(c,li))
          == PT3_deg.geom_traits().construct_triangle_3_object()(
            Point_3(2,0,2), Point_3(2,0,0), Point_3(2,2,0)) );
  assert(PT3_deg.side_of_cell(Point(2,0.5,0.5),c,lt,li,lj)
         == CGAL::ON_BOUNDARY);
  assert(lt == P3T3::FACET);
  assert( PT3_deg.construct_triangle(PT3_deg.periodic_triangle(c,li))
          == PT3_deg.geom_traits().construct_triangle_3_object()(
            Point_3(2,0,2), Point_3(2,0,0), Point_3(2,2,0)) );

  c = PT3_deg.locate(Point(2,2,1),lt,li,lj);
  assert(lt == P3T3::EDGE);
  assert( PT3_deg.construct_segment(PT3_deg.periodic_segment(c,li,lj))
          == PT3_deg.geom_traits().construct_segment_3_object()(
            Point_3(2,2,0), Point_3(2,2,2)) );
  assert(PT3_deg.side_of_cell(Point(2,2,1),c,lt,li,lj) == CGAL::ON_BOUNDARY);
  assert(lt == P3T3::EDGE);
  assert( PT3_deg.construct_segment(PT3_deg.periodic_segment(c,li,lj))
          == PT3_deg.geom_traits().construct_segment_3_object()(
            Point_3(2,2,0), Point_3(2,2,2)) );

  c = PT3_deg.locate(Point(2,2,2),lt,li,lj);
  assert(lt == P3T3::VERTEX);
  assert(c->vertex(li)->point() == Point(2,2,2));
  assert(PT3_deg.side_of_cell(Point(2,2,2),c,lt,li,lj) == CGAL::ON_BOUNDARY);
  assert(lt == P3T3::VERTEX);
  assert(c->vertex(li)->point() == Point(2,2,2));

  c = P3T3().locate(Point(1,2,3),lt,li,lj);
  assert(c == Cell_handle());
  assert(lt == P3T3::EMPTY);
  assert(P3T3().side_of_cell(Point(1,2,3),c,lt,li,lj)
         == CGAL::ON_UNBOUNDED_SIDE);
  assert(lt == P3T3::EMPTY);

  std::cout << "Testing Iterators   "<< std::endl;

  _test_vertex_iterator(PT3);
  _test_vertex_iterator(PT1);
  _test_vertex_iterator(PT3_deg);
  _test_vertex_iterator(PT1_deg);

  _test_unique_vertex_iterator(PT3);
  _test_unique_vertex_iterator(PT1);
  _test_unique_vertex_iterator(PT3_deg);
  _test_unique_vertex_iterator(PT1_deg);

  _test_triangulation_iterator(PT3);
  _test_triangulation_iterator(PT1);
  _test_triangulation_iterator(PT3_deg);
  _test_triangulation_iterator(PT1_deg);

  std::cout << "Testing Circulator  "<< std::endl;

  _test_circulator(PT3);
  _test_circulator(PT1);
  _test_circulator(PT3_deg);
  _test_circulator(PT1_deg);

  std::cout << "Incidence and adjacency" << std::endl;
  std::cout << "  in 3-sheeted covering space" << std::endl;

  std::vector<Cell_handle> cellv;
  std::vector<Facet> facetv;
  std::vector<Edge> edgev;
  std::vector<Vertex_handle> vertexv;
  PT3.incident_cells(PT3.vertices_begin(),std::back_inserter(cellv));
  for (unsigned int n=0 ; n<cellv.size() ; n++) {
    assert( (PT3.vertices_begin() == cellv[n]->vertex(0))
            || (PT3.vertices_begin() == cellv[n]->vertex(1))
            || (PT3.vertices_begin() == cellv[n]->vertex(2))
            || (PT3.vertices_begin() == cellv[n]->vertex(3)) );
    assert(PT3.is_cell(cellv[n]));
    assert(PT3.is_cell(cellv[n]->vertex(0),cellv[n]->vertex(1),
                       cellv[n]->vertex(2),cellv[n]->vertex(3),c));
    assert(PT3.is_cell(
             cellv[n]->vertex(0),PT3.periodic_point(cellv[n],0).second,
             cellv[n]->vertex(1),PT3.periodic_point(cellv[n],1).second,
             cellv[n]->vertex(2),PT3.periodic_point(cellv[n],2).second,
             cellv[n]->vertex(3),PT3.periodic_point(cellv[n],3).second,c));
  }

  PT3.incident_facets(PT3.vertices_begin(),std::back_inserter(facetv));
  for (unsigned int n=0 ; n<facetv.size() ; n++) {
    assert( (PT3.vertices_begin()
             == facetv[n].first->vertex((facetv[n].second+1)&3))
            || (PT3.vertices_begin()
                == facetv[n].first->vertex((facetv[n].second+2)&3))
            || (PT3.vertices_begin()
                == facetv[n].first->vertex((facetv[n].second+3)&3)) );
    assert(PT3.is_facet(facetv[n].first->vertex((facetv[n].second+1)&3),
                        facetv[n].first->vertex((facetv[n].second+2)&3),
                        facetv[n].first->vertex((facetv[n].second+3)&3),c,i,j,k));
    assert(PT3.is_facet(facetv[n].first->vertex((facetv[n].second+1)&3),
                        PT3.periodic_point(facetv[n].first,(facetv[n].second+1)&3).second,
                        facetv[n].first->vertex((facetv[n].second+2)&3),
                        PT3.periodic_point(facetv[n].first,(facetv[n].second+2)&3).second,
                        facetv[n].first->vertex((facetv[n].second+3)&3),
                        PT3.periodic_point(facetv[n].first,(facetv[n].second+3)&3).second,
                        c,i,j,k));
  }
  PT3.incident_edges(PT3.vertices_begin(),std::back_inserter(edgev));
  for (unsigned int n=0 ; n<edgev.size() ; n++) {
    assert( (PT3.vertices_begin()
             == edgev[n].first->vertex(edgev[n].second))
            || (PT3.vertices_begin()
                == edgev[n].first->vertex(edgev[n].third)) );
    assert(PT3.is_edge(edgev[n].first->vertex(edgev[n].second),
                       edgev[n].first->vertex(edgev[n].third),c,i,j));
    assert(PT3.is_edge(edgev[n].first->vertex(edgev[n].second),
                       PT3.periodic_point(edgev[n].first,edgev[n].second).second,
                       edgev[n].first->vertex(edgev[n].third),
                       PT3.periodic_point(edgev[n].first,edgev[n].third).second,c,i,j));
  }
  PT3.adjacent_vertices(PT3.vertices_begin(),std::back_inserter(vertexv));
  for (unsigned int n=0 ; n<vertexv.size() ; n++) {
    assert(PT3.is_vertex(vertexv[n]));
    assert(PT3.is_edge(PT3.vertices_begin(),vertexv[n],c,i,j));
  }

  assert(PT3.degree(PT3.vertices_begin()) == vertexv.size());

  std::cout << "  in 3-sheeted covering space" << std::endl;

  cellv.clear();
  facetv.clear();
  edgev.clear();
  vertexv.clear();
  PT1.incident_cells(PT1.vertices_begin(),std::back_inserter(cellv));
  for (unsigned int n=0 ; n<cellv.size() ; n++) {
    assert( (PT1.vertices_begin() == cellv[n]->vertex(0))
            || (PT1.vertices_begin() == cellv[n]->vertex(1))
            || (PT1.vertices_begin() == cellv[n]->vertex(2))
            || (PT1.vertices_begin() == cellv[n]->vertex(3)) );
    assert(PT1.is_cell(cellv[n]));
    assert(PT1.is_cell(cellv[n]->vertex(0),cellv[n]->vertex(1),
                       cellv[n]->vertex(2),cellv[n]->vertex(3),c));
    assert(PT1.is_cell(
             cellv[n]->vertex(0),PT1.periodic_point(cellv[n],0).second,
             cellv[n]->vertex(1),PT1.periodic_point(cellv[n],1).second,
             cellv[n]->vertex(2),PT1.periodic_point(cellv[n],2).second,
             cellv[n]->vertex(3),PT1.periodic_point(cellv[n],3).second,c));
  }

  PT1.incident_facets(PT1.vertices_begin(),std::back_inserter(facetv));
  for (unsigned int n=0 ; n<facetv.size() ; n++) {
    assert( (PT1.vertices_begin()
             == facetv[n].first->vertex((facetv[n].second+1)&3))
            || (PT1.vertices_begin()
                == facetv[n].first->vertex((facetv[n].second+2)&3))
            || (PT1.vertices_begin()
                == facetv[n].first->vertex((facetv[n].second+3)&3)) );
    assert(PT1.is_facet(facetv[n].first->vertex((facetv[n].second+1)&3),
                        facetv[n].first->vertex((facetv[n].second+2)&3),
                        facetv[n].first->vertex((facetv[n].second+3)&3),c,i,j,k));
    assert(PT1.is_facet(facetv[n].first->vertex((facetv[n].second+1)&3),
                        PT1.periodic_point(facetv[n].first,(facetv[n].second+1)&3).second,
                        facetv[n].first->vertex((facetv[n].second+2)&3),
                        PT1.periodic_point(facetv[n].first,(facetv[n].second+2)&3).second,
                        facetv[n].first->vertex((facetv[n].second+3)&3),
                        PT1.periodic_point(facetv[n].first,(facetv[n].second+3)&3).second,
                        c,i,j,k));
  }
  PT1.incident_edges(PT1.vertices_begin(),std::back_inserter(edgev));
  for (unsigned int n=0 ; n<edgev.size() ; n++) {
    assert( (PT1.vertices_begin()
             == edgev[n].first->vertex(edgev[n].second))
            || (PT1.vertices_begin()
                == edgev[n].first->vertex(edgev[n].third)) );
    assert(PT1.is_edge(edgev[n].first->vertex(edgev[n].second),
                       edgev[n].first->vertex(edgev[n].third),c,i,j));
    assert(PT1.is_edge(edgev[n].first->vertex(edgev[n].second),
                       PT1.periodic_point(edgev[n].first,edgev[n].second).second,
                       edgev[n].first->vertex(edgev[n].third),
                       PT1.periodic_point(edgev[n].first,edgev[n].third).second,c,i,j));
  }
  PT1.adjacent_vertices(PT1.vertices_begin(),std::back_inserter(vertexv));
  for (unsigned int n=0 ; n<vertexv.size() ; n++) {
    assert(PT1.is_vertex(vertexv[n]));
    assert(PT1.is_edge(PT1.vertices_begin(),vertexv[n],c,i,j));
  }

  assert(PT1.degree(PT1.vertices_begin()) == vertexv.size());

  std::cout << "Mirroring" << std::endl;

  assert(ch->neighbor(0)->neighbor(PT3.mirror_index(ch,0)) == ch);
  assert(ch->neighbor(0)->neighbor(
           ch->neighbor(0)->index(PT3.mirror_vertex(ch,0))) == ch);

  assert(ch->neighbor(0)->vertex(PT3.mirror_facet(Facet(ch,0)).second)
         == PT3.mirror_vertex(ch,0));
  assert(PT3.mirror_facet(Facet(ch,0)).first == ch->neighbor(0));


  // There are problems with the IO of exact number types in binary mode.
  if(test_input_ouput)
  {
    std::cout << "I/O" << std::endl;
    std::cout << "  ascii" << std::endl;

    std::stringstream ss1;
    std::stringstream ss3;
    ss1 << PT1;
    ss3 << PT3;

    P3T3 PT1r, PT3r;
    ss1 >> PT1r;
    ss3 >> PT3r;

    assert(CGAL::is_ascii(ss1));
    assert(CGAL::is_ascii(ss3));
    if (!ex) assert(PT1 == PT1r);
    if (!ex) assert(PT3 == PT3r);

    std::cout << "  binary" << std::endl;

    PT1r.clear();
    PT3r.clear();
    if (!ex) {
      std::stringstream ss1b;
      std::stringstream ss3b;
      CGAL::set_binary_mode(ss1b);
      CGAL::set_binary_mode(ss3b);
      ss1b << PT1;
      ss3b << PT3;

      ss1b >> PT1r;
      ss3b >> PT3r;
      assert(CGAL::is_binary(ss1b));
      assert(CGAL::is_binary(ss3b));

      assert(PT1 == PT1r);
      assert(PT3 == PT3r);
    }

    std::cout << "  pretty" << std::endl;

    PT1r.clear();
    PT3r.clear();
    std::stringstream ss1p;
    std::stringstream ss3p;
    CGAL::set_pretty_mode(ss1p);
    CGAL::set_pretty_mode(ss3p);
    ss1p << PT1;
    ss3p << PT3;

    assert(CGAL::is_pretty(ss1p));
    assert(CGAL::is_pretty(ss3p));
  }
}
