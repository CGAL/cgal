// Copyright (c) 1998-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Francois Rebufat, Monique Teillaud, Sylvain Pion
//                 Mariette Yvinec

#ifndef CGAL_TEST_CLS_DELAUNAY_C
#define CGAL_TEST_CLS_DELAUNAY_C

#include <cassert>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "_test_cls_iterator.h"
#include "_test_cls_circulator.h"
#include "_test_remove_cluster.h"

#include <CGAL/Random.h>
#include <CGAL/Testsuite/use.h>

// Accessory set of functions to differentiate between
// Delaunay::nearest_vertex[_in_cell] and
//  Regular::nearest_power_vertex[_in_cell].

template < typename T, typename Weighted_tag = typename T::Weighted_tag >
struct Test_location_policy {
	typedef typename T::Location_policy Location_policy;
};

template < typename T >
struct Test_location_policy <T, CGAL::Tag_true> {
	struct Location_policy{};
};

template < typename T, typename P >
typename T::Vertex_handle
nearest_vertex(const T&t, const P&p, CGAL::Tag_true)
{
  return t.nearest_power_vertex(p);
}

template < typename T, typename P >
typename T::Vertex_handle
nearest_vertex(const T&t, const P&p, CGAL::Tag_false)
{
  return t.nearest_vertex(p);
}

template < typename T, typename P >
typename T::Vertex_handle
nearest_vertex(const T&t, const P&p)
{
  return nearest_vertex(t, p, typename T::Weighted_tag());
}

template < typename T, typename P >
typename T::Vertex_handle
nearest_vertex_in_cell(const T&t, const P&p, const typename T::Cell_handle c, CGAL::Tag_true)
{
  return t.nearest_power_vertex_in_cell(p, c);
}

template < typename T, typename P >
typename T::Vertex_handle
nearest_vertex_in_cell(const T&t, const P&p, const typename T::Cell_handle c, CGAL::Tag_false)
{
  return t.nearest_vertex_in_cell(p, c);
}

template < typename T, typename P >
typename T::Vertex_handle
nearest_vertex_in_cell(const T&t, const P&p, const typename T::Cell_handle c)
{
  return nearest_vertex_in_cell(t, p, c, typename T::Weighted_tag());
}

// Template meta programming if.
template < typename Cond, typename Then, typename Else >
struct If;

template < typename Then, typename Else >
struct If <CGAL::Tag_true, Then, Else>
{
  typedef Then type;
};

template < typename Then, typename Else >
struct If <CGAL::Tag_false, Then, Else>
{
  typedef Else type;
};

template < typename T, typename P >
void test_conflicts(T& T3_13, const P *q)
{
  typedef T Cls;
  typedef typename Cls::Locate_type   Locate_type;
  typedef typename Cls::Vertex_handle Vertex_handle;
  typedef typename Cls::Cell_handle   Cell_handle;
  typedef typename Cls::Facet         Facet;

  for (int i=0; i<22; ++i) {
    if (T3_13.dimension() < 2)
      T3_13.insert(q[i]);
    else {
      // First locate the point.
      Locate_type lt;
      int li, lj;
      Cell_handle c = T3_13.locate(q[i], lt, li, lj);
      if (lt == Cls::VERTEX) // Already exist, skip.
        continue;
      if (lt == Cls::OUTSIDE_AFFINE_HULL) { // Increases dimension.
        T3_13.insert_outside_affine_hull(q[i]);
        continue;
      }
      // Get the stuff in conflicts.
      std::vector<Cell_handle>   C;
      std::vector<Facet>         F;
      std::vector<Vertex_handle> V;

      T3_13.vertices_on_conflict_zone_boundary(q[i], c, std::back_inserter(V));
#ifndef CGAL_NO_DEPRECATED_CODE
      // test deprecated vertices_in_conflict
      std::vector<Vertex_handle> V2;
      T3_13.vertices_in_conflict(q[i], c, std::back_inserter(V2));
      assert(V2.size() == V.size());
#endif

      T3_13.find_conflicts(q[i], c, std::back_inserter(F),
                           std::back_inserter(C));

      if (T3_13.dimension() == 3)
          assert(F.size() == 2*V.size() - 4); // Euler relation.
      if (T3_13.dimension() == 2)
          assert(F.size() == V.size());

      if (i%2 == 0)
          T3_13.insert_in_hole(q[i], C.begin(), C.end(),
                               F.begin()->first, F.begin()->second);
      else {
	  // alternately test the overload which takes a Vertex_handle.
	  Vertex_handle v = T3_13.tds().create_vertex();
          T3_13.insert_in_hole(q[i], C.begin(), C.end(),
                               F.begin()->first, F.begin()->second, v);
      }
    }
  }
}

template <class Triangulation>
void
_test_cls_delaunay_3(const Triangulation &)
{
  typedef Triangulation                      Cls;

  typedef typename Test_location_policy<Cls>::Location_policy Location_policy;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested

  // typedef typename Cls::Point          Point; // Delaunay
  // typedef typename Cls::Point::Point   Point; // Regular
  typedef typename If<typename Cls::Weighted_tag,
                      typename Cls::Point, Cls>::type::Point   Point; 

  typedef typename Cls::Segment              Segment;
  typedef typename Cls::Triangle             Triangle;
  typedef typename Cls::Tetrahedron          Tetrahedron;

  typedef typename Cls::Vertex               Vertex;
  typedef typename Cls::Cell                 Cell;
  typedef typename Cls::Facet                Facet;
  typedef typename Cls::Edge                 Edge;

  typedef typename Cls::size_type            size_type;

  typedef typename Cls::Vertex_handle        Vertex_handle;
  typedef typename Cls::Cell_handle          Cell_handle; 
  typedef typename Cls::Vertex_iterator      Vertex_iterator;
  typedef typename Cls::Cell_iterator        Cell_iterator;
  typedef typename Cls::Locate_type          Locate_type;
  typedef std::list<Point>                        list_point;

  typedef typename Cls::Finite_vertices_iterator      Finite_vertices_iterator;
  typedef typename Cls::Finite_cells_iterator        Finite_cells_iterator;

  CGAL_USE_TYPE(Location_policy);
  CGAL_USE_TYPE(Segment);
  CGAL_USE_TYPE(Triangle);
  CGAL_USE_TYPE(Tetrahedron);
  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Cell);
  CGAL_USE_TYPE(Vertex_iterator);
  CGAL_USE_TYPE(Cell_iterator);
  // +++ We define now some points for building triangulations +++++//

// list of Points for T1_0 , T1_1, T1_2 :


  Point p[5]=
  { Point(0,0,0),Point(0,0,1),Point(0,0,2),Point(0,0,3),Point(0,0,4)};
  list_point l1;
  int i;
  for (i=0;i<5;i++)
    l1.push_back(p[i]);

  Point pp[5] = { Point(0,0,0),Point(0,0,4),Point(0,0,1),Point(0,0,2),Point(0,0,3)};
  list_point l2;
  for (i=0;i<5;i++)
    l2.push_back(pp[i]);

  Point ppp[5] ={Point(0,0,1),Point(0,0,2),Point(0,0,3),Point(0,0,0),Point(0,0,4)};
  list_point l3;
  for (i=0;i<5;i++)
    l3.push_back(ppp[i]);

// Points for T2_0 :
  Point p1=Point(5,5,0); 
  Point p2=Point(4,4,0);
  Point p3=Point(6,6,0); // 1- dimensional until this point
  Point p4=Point(5,3,0); // 2-dimensional
  Point p5=Point(5,7,0); 
  Point p6=Point(5,4,0); 
  Point p7=Point(5,6,0); 
  Point p8=Point(0,0,0); 
  Point p9=Point(5,5,0); 

  // Points for T3_1 :
  Point q[22] = 
  {
    Point(0,0,0), Point(4,4,0), Point(0,4,0), Point(4,0,0),
    Point(1,3,1), Point(3,1,1), Point(3,3,1), Point(1,1,1), Point(2,2,2),
    Point(1,3,3), Point(3,1,3), Point(3,3,3), Point(1,1,3), 
    Point(0,0,4), Point(4,4,4), Point(0,4,4), Point(4,0,4),
    Point(1,3,5), Point(3,1,5), Point(3,3,5), Point(1,1,5), Point(2,2,6)};

  // Points for T3_2 :

  list_point lp;
  int a, b, d;
//   for (a=0;a!=10;a++)
//     for (b=0;b!=10;b++)
//       for (d=0;d!=10;d++)
// 	lp.push_back(Point(a*b-d*a + (a-b)*10 +a ,a-b+d +5*b,
// 			   a*a-d*d+b));

  for (a=0;a!=10;a++)
    for (b=0;b!=10;b++)
      for (d=0;d!=5;d++)
	lp.push_back(Point(a*b-d*a + (a-b)*10 +a ,a-b+d +5*b,
			   a*a-d*d+b));

  // Points for T3_2 :

  list_point lp2;
  for (a=0;a!=4;a++)
    for (b=0;b!=4;b++)
      for (d=0;d!=4;d++)
	lp2.push_back(Point((a*b-d*a)*10 +a ,(a-b+d +5*b)*100,
			    a*a-d*d-b));

  //########################################################################


  /**************CONSTRUCTORS (1)*********************/
  /************** and I/O ****************************/

  std::cout << "    Constructor " << std::endl;
  // Begining with an empty triangulation and adding point until reaching
  // 3-dimentional triangulation.
  Cls T0; 
  assert(T0.dimension() == -1);
  assert(T0.number_of_vertices() == 0);
  assert(T0.is_valid());

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT1("Test1_triangulation_IO_3",std::ios::out);
        oFileT1 << T0 << std::endl;
      }
      std::ifstream iFileT1("Test1_triangulation_IO_3",std::ios::in);
      iFileT1 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == -1);
      assert(Tfromfile.number_of_vertices() == 0);
  }

  std::cout << "    Constructor1 " << std::endl;
  Point p10(0,0,0);
  Vertex_handle v0=T0.insert(p10);
  assert(T0.dimension() == 0);
  assert(T0.number_of_vertices() == 1);
  assert(T0.is_valid());

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
	std::ofstream oFileT2("Test2_triangulation_IO_3",std::ios::out);
        oFileT2 << T0 << std::endl;
      }
      std::ifstream iFileT2("Test2_triangulation_IO_3",std::ios::in);
      iFileT2 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 0);
      assert(Tfromfile.number_of_vertices() == 1);
  }

  std::cout << "    Constructor2 " << std::endl;

  Point p11(100,100,0);
  v0=T0.insert(p11);
  assert(T0.dimension() == 1);
  assert(T0.number_of_vertices() == 2);
  assert(T0.is_valid());

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT3("Test3_triangulation_IO_3",std::ios::out);
        oFileT3 << T0 << std::endl;
      }
      std::ifstream iFileT3("Test3_triangulation_IO_3",std::ios::in);
      iFileT3 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 1);
      assert(Tfromfile.number_of_vertices() == 2);
  }

  std::cout << "    Constructor3 " << std::endl;

  Point p12(100,-100,0);
  v0=T0.insert(p12);
  assert(T0.dimension() == 2);
  assert(T0.number_of_vertices() == 3);
  assert(T0.is_valid());

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT4("Test4_triangulation_IO_3",std::ios::out);
        oFileT4 << T0;
      }
      std::ifstream iFileT4("Test4_triangulation_IO_3",std::ios::in);
      iFileT4 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 2);
      assert(Tfromfile.number_of_vertices() == 3);
  }

  std::cout << "    Constructor4 " << std::endl;

  Point p13(50,0,100);
  v0=T0.insert(p13);
  assert(T0.dimension() == 3);
  assert(T0.number_of_vertices() == 4);
  assert(T0.is_valid());

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT5("Test5_triangulation_IO_3",std::ios::out);
        oFileT5 << T0;
      }
      std::ifstream iFileT5("Test5_triangulation_IO_3",std::ios::in);
      iFileT5 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 3);
      assert(Tfromfile.number_of_vertices() == 4);
  }

  std::cout << "    Constructor5 " << std::endl;

  Point p14(50,0,100);
  v0=T0.insert(p14);
  assert(T0.dimension() == 3);
  assert(T0.number_of_vertices() == 4);
  assert(T0.is_valid());
  
  // copy constructor
  Cls T1(T0);
  assert(T1.dimension() == 3);
  assert(T1.number_of_vertices() == 4);
  assert(T1.is_valid());
   
  T1.clear();
  assert(T1.dimension() == -1);
  assert(T1.number_of_vertices() == 0);
  assert(T1.is_valid());



   // Affectation :
  T1=T0;  
  assert(T1.dimension() == 3);
  assert(T1.number_of_vertices() == 4);
  assert(T1.is_valid());
  T1.clear();


  T1.swap(T0);
  assert(T1.dimension() == 3);
  assert(T1.number_of_vertices() == 4);
  assert(T1.is_valid());
  assert(T0.dimension() == -1);
  assert(T0.number_of_vertices() == 0);
  assert(T0.is_valid());
  T0.swap(T1);
  
  assert(T0.dimension() == 3);
  assert(T0.number_of_vertices() == 4);
  assert(T0.is_valid());
  assert(T1.dimension() == -1);

     // Building some psychotics triangulations :
  std::cout << "    Constructor6 " << std::endl;
  // triangulation 1-dimensional : vertical line.
  Cls T1_0;
  size_type n = T1_0.insert(l1.begin(),l1.end());
  assert(n==5);
  assert(T1_0.dimension()==1);
  assert(T1_0.number_of_vertices()==n);
  assert(T1_0.is_valid());
  std::cout << "    Constructor7 " << std::endl;
  Cls T1_1;
  n = T1_1.insert(l2.begin(),l2.end());
  assert(n==5);
  assert(T1_1.dimension()==1);
  assert(T1_1.number_of_vertices()==n);
  assert(T1_1.is_valid());
  std::cout << "    Constructor8 " << std::endl;
  Cls T1_2;
  n = T1_2.insert(l3.begin(),l3.end());
  assert(n==5);
  assert(T1_2.dimension()==1);
  assert(T1_2.number_of_vertices()==n);
  assert(T1_2.is_valid());

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT6("Test6_triangulation_IO_3",std::ios::out);
        oFileT6 << T1_2;
      }
      std::ifstream iFileT6("Test6_triangulation_IO_3",std::ios::in);
      iFileT6 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 1);
      assert(Tfromfile.number_of_vertices() == n);
  }

  std::cout << "    Constructor9 " << std::endl;
  // 2-dimensional triangulations 

  Cls T2_0;
  v0=T2_0.insert(p1);

  v0=T2_0.insert(p2);

  v0=T2_0.insert(p3);

  assert(T2_0.is_valid());
  assert(T2_0.dimension()==1);
  assert(T2_0.number_of_vertices()==3);


  v0=T2_0.insert(p4);
  assert(T2_0.is_valid());
  assert(T2_0.dimension()==2);
  assert(T2_0.number_of_vertices()==4);


  v0=T2_0.insert(p5);
  v0=T2_0.insert(p6);
  v0=T2_0.insert(p7);
  v0=T2_0.insert(p8);
  v0=T2_0.insert(p9);

  assert(T2_0.is_valid());
  assert(T2_0.dimension()==2);
  assert(T2_0.number_of_vertices()==8);

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT7("Test7_triangulation_IO_3",std::ios::out);
        oFileT7 << T2_0;
      }
      std::ifstream iFileT7("Test7_triangulation_IO_3",std::ios::in);
      iFileT7 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 2);
      assert(Tfromfile.number_of_vertices() == 8);
  }

  std::cout << "    Constructor10 " << std::endl;
  // test grid insert
  Cls T2_1;
  int m;
  int px=1, py=1;
  int qx=-1, qy=2;
  Point qq[400];
  for (m=0; m<20; m++)
    for (n=0; n<20; n++)
      {
	qq[m+20*n] = Point(m*px+(int)n*qx, m*py+(int)n*qy, 1);
	T2_1.insert( qq[m+20*n] );
      }
  assert( T2_1.number_of_vertices() == m*n );
  assert( T2_1.dimension()==2 );
  assert( T2_1.is_valid() );

  std::cout << "    Constructor11 " << std::endl;
  // 3-dimensional triangulations
  // This is a simple grid :
  int x,y,z;
  Cls T3_0;
  Point r[225];
  for (z=0 ; z<5 ; z++)
    for (y=0 ; y<5 ; y++)
      for (x=0 ; x<5 ; x++) 
	{
	  r[x+5*y+25*z] = Point(x,y,z);
	  v0=T3_0.insert(r[x+5*y+25*z]);
	}
  assert(T3_0.is_valid());
  assert(T3_0.number_of_vertices()==125);
  assert(T3_0.dimension()==3);

  if (del) {
    std::cout << "    deletion in Delaunay - grid case - (dim 3) " <<
      std::endl; 
    Cls Tdel( T3_0 );
    
    std::vector<Vertex_handle> vertices;
    for (Finite_vertices_iterator vi = Tdel.finite_vertices_begin();
	 vi != Tdel.finite_vertices_end(); ++vi)
      vertices.push_back(vi);

    size_type n = Tdel.number_of_vertices();
    size_type m = Tdel.remove(vertices.begin(), vertices.end());
    assert(m == n - Tdel.number_of_vertices());
    assert(Tdel.is_valid(false));
    std::cout << "    successfull" << std::endl; 
  }


  std::cout << "    Constructor12 " << std::endl;
  Cls T3_1;
  for (i=0;i<22;i++)
    T3_1.insert(q[i]);
  assert(T3_1.is_valid());
  assert(T3_1.number_of_vertices()==22);
  assert(T3_1.dimension()==3);

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT8("Test8_triangulation_IO_3",std::ios::out);
        oFileT8 << T3_1;
      }
      std::ifstream iFileT8("Test8_triangulation_IO_3",std::ios::in);
      iFileT8 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 3);
      assert(Tfromfile.number_of_vertices() == 22);
  }

  // Testing find_conflicts(), vertices_on_conflict_zone_boundary(), insert_in_hole()
  // FIXME : Note that we do not test the version of find_conflicts()
  //         which returns the internal facets too...
  std::cout << "    Testing find_conflicts/vertices_on_conflict_zone_boundary/insert_in_hole"
            << std::endl;
  Cls T3_13;

  test_conflicts(T3_13, q);

  assert(T3_13.is_valid());
  assert(T3_13.number_of_vertices()==22);
  assert(T3_13.dimension()==3);

#ifndef CGAL_NO_DEPRECATED_CODE
  {
    std::cout << "    Testing move_point()" << std::endl;
    Cls T;
    std::list<Vertex_handle> L;
    for (i=0; i<22; ++i)
      L.push_back(T.insert(q[i]));
    assert(T.is_valid());
    assert(T.number_of_vertices()==22);
    assert(T.dimension()==3);

    for (i=0; i<100; ++i) {
      assert(!L.empty());
      Vertex_handle v = L.front();
      L.pop_front();
      size_type nbv = T.number_of_vertices();
      L.push_back(T.move_point(v, q[(3*i)%22]));

      if (nbv != T.number_of_vertices())
        L.pop_back(); // it means we move onto an already existing point.

      assert(T.is_valid());
      assert(T.number_of_vertices()<=22);
    }
  }
#endif

  {
    std::cout << "    Testing move()" << std::endl;
    Cls T;
    std::list<Vertex_handle> L;
    for (i=0; i<22; ++i)
      L.push_back(T.insert(q[i]));
    assert(T.is_valid());
    assert(T.number_of_vertices()==22);
    assert(T.dimension()==3);

    for (i=0; i<100; ++i) {
      assert(!L.empty());
      Vertex_handle v = L.front();
      L.pop_front();
      size_type nbv = T.number_of_vertices();
      L.push_back(T.move(v, q[(3*i)%22]));

      if (nbv != T.number_of_vertices())
        L.pop_back(); // it means we move onto an already existing point.

      assert(T.is_valid());
      assert(T.number_of_vertices()<=22);
    }
  }

  {
    std::cout << "    Testing nearest_vertex()" << std::endl;
    // We do a nearest_vertex() and two nearest_vertex_in_cell()
    // queries on all points with integer coordinate
    // in the cube [-1;6]^3. In each case we check explicitely that the
    // output is correct by comparing distance to other vertices.
    Cell_handle c1 = T3_13.finite_cells_begin();
    Cell_handle c2 = T3_13.infinite_vertex()->cell();
    for (int x = -1; x < 7; ++x)
      for (int y = -1; y < 7; ++y)
	for (int z = -1; z < 7; ++z) {
	  Point p(x, y, z);
	  Vertex_handle v = nearest_vertex(T3_13, p);
	  for (typename Cls::Finite_vertices_iterator
	         fvit = T3_13.finite_vertices_begin();
	       fvit != T3_13.finite_vertices_end(); ++fvit)
	    assert(CGAL::squared_distance(p, v->point()) <=
		   CGAL::squared_distance(p, fvit->point()));
	  Vertex_handle v1 = nearest_vertex_in_cell(T3_13, p, c1);
	  int i1 = c1->index(v1);
 	  for(int i=0; i<4; ++i) {
	    if (i != i1) 
	      assert(CGAL::squared_distance(p, v1->point()) <=
		     CGAL::squared_distance(p, c1->vertex(i)->point()));
	  }
	  Vertex_handle v2 = nearest_vertex_in_cell(T3_13, p, c2);
	  int i2 = c2->index(v2);
	  for(int i=0; i<4; ++i) { 
	    if (i != i2 && c2->vertex(i) != T3_13.infinite_vertex())
	      assert(CGAL::squared_distance(p, v2->point()) <=
		     CGAL::squared_distance(p, c2->vertex(i)->point()));
	  }
	}
  }

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      {
        std::ofstream oFileT8("Test13_triangulation_IO_3",std::ios::out);
        oFileT8 << T3_13;
      }
      std::ifstream iFileT8("Test13_triangulation_IO_3",std::ios::in);
      iFileT8 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 3);
      assert(Tfromfile.number_of_vertices() == 22);
  }


  //#######################################################################
  std::cout << "  list insertion"<< std::endl;
  Cls T3_2_1;
  T3_2_1.insert(lp2.begin(),lp2.end());
  assert(T3_2_1.is_valid());

  assert(T3_2_1.dimension()==3);
  //  assert(T3_2.number_of_vertices()==1000);
  std::cout << "   end of insertion " << std::endl;

  std::cout << "  insertion of located point"<< std::endl;
  Locate_type lt;
  int li, lj;
  Cell_handle ch = T3_2_1.locate(p7,lt,li,lj);
  T3_2_1.insert(p7,lt,ch,li,lj);

  assert(T3_2_1.is_valid());
  assert(T3_2_1.dimension()==3);

  // Same as above but using the template ctor.
  std::cout << "  template constructor"<< std::endl;
  Cls T3_2_2(lp2.begin(), lp2.end());
  assert(T3_2_2.is_valid());

  assert(T3_2_2.dimension()==3);
  std::cout << "   end of insertion " << std::endl;

  //########################################################################


  std::cout << "  500 points insertion"<< std::endl;
  Cls T3_2;
  typename list_point::iterator it;
  int count = 0 ;
  std::cout << " number of inserted points : " ;
  for (it=lp.begin(); it!=lp.end();it++){
    count++;
    T3_2.insert(*it);
    if (count <10)
      std::cout << count << '\b' ;
    else
      if (count < 100)
        std::cout << count << '\b' << '\b' ;
      else 
        if (count < 1000)
          std::cout << count << '\b' << '\b' << '\b' ;
        else
	  std::cout << count << std::endl;
    std::cout.flush();
  }
  std::cout << std::endl;
  assert(T3_2.is_valid());
  assert(T3_2.dimension()==3);
  assert(T3_2.number_of_vertices()==500);
 
 

  Point p110(-5,5,0), p111(-2,-5,2), p112(-2,-9,6), p113(4,8,9), p114(5,-6,0),
    p115(3,0,5), p116(-9,0,-10), p117(1,6,-2), p118(-3,2,-4), p119(3,-3,-1);
  Cls T3_5;
  v0=T3_5.insert(p110);
  v0=T3_5.insert(p111);
  v0=T3_5.insert(p112);
  v0=T3_5.insert(p113);
  v0=T3_5.insert(p114);
  v0=T3_5.insert(p115);
  v0=T3_5.insert(p116);
  v0=T3_5.insert(p117);
  v0=T3_5.insert(p118);
  v0=T3_5.insert(p119);

  assert(T3_5.is_valid());
  assert(T3_5.number_of_vertices()==10);

  // %%%%%%%%%% deletion in Delaunay
  if (del) {
    std::cout << "    deletion in a 10 points Delaunay triangulation";
    Vertex_handle v;
    while ( T3_5.number_of_vertices() >= 1 ) {
      if ( T3_5.dimension() == 3 )
	v = T3_5.infinite_cell()->vertex
	  ( (T3_5.infinite_cell()->index( T3_5.infinite_vertex() ) +1 )&3 );
      else if ( T3_5.dimension() == 2 )
	v = T3_5.infinite_cell()->vertex
	  ( (T3_5.infinite_cell()->index( T3_5.infinite_vertex() ) +1 )%3 );
      else if ( T3_5.dimension() == 1 )
	  v = T3_5.infinite_cell()->vertex
	    ( (T3_5.infinite_cell()->index( T3_5.infinite_vertex() ) +1 )%2 );
	else
	  v = T3_5.infinite_cell()->neighbor(0)->vertex(0);

      T3_5.remove( v );
    }
    assert(T3_5.is_valid(false));
  }
  std::cout << " done" << std::endl;

//  // Test random triangulation :

//   Cls T3_4;
//   CGAL::Random random;
//   for (n=1;n<50;n++) {
//     x=random.get_int(-500,500);
//     y=random.get_int(-500,500);
//     z=random.get_int(-500,500);
//     v0=T3_4.insert(Point(x,y,z));
//   }
//   assert(T3_4.is_valid());
//   assert(T3_4.dimension()==3);
//   assert(T3_4.is_valid());

//   // %%%%%%%%%% deletion in Delaunay
//   bool success(true);
//   if (del) {
//     std::cout << "    deletion in a Delaunay of "
// 	      << T3_4.number_of_vertices() << " random points";
//     Vertex_handle v;
//     while ( T3_4.number_of_vertices() >= 1 ) {
//       if ( T3_4.dimension() > 1 )
// 	v = T3_4.infinite_cell()->vertex
// 	  ( (T3_4.infinite_cell()->index( T3_4.infinite_vertex() ) +1 )&3 );
//       else
// 	if ( T3_4.dimension() == 1 )
// 	  v = T3_4.infinite_cell()->vertex
// 	    ( (T3_4.infinite_cell()->index( T3_4.infinite_vertex() ) +1 )%2 );
// 	else
// 	  v = T3_4.infinite_cell()->neighbor(0)->vertex(0);

//       success = T3_4.remove( v );
//     }
//     if (success) assert(T3_4.is_valid(false));
//   }
//   std::cout << " done" << std::endl;

       // Iterator and circulator test

  Cls T0_1;
  Point p28(1,3,5);
  v0=T0_1.insert(p28);

  {
      std::cout << "    Testing Iterator   "<< std::endl;
      _test_vertex_iterator(T0_1);
      _test_triangulation_iterator(T0_1);
      _test_vertex_iterator(T0);
      _test_triangulation_iterator(T0);
      _test_vertex_iterator(T2_0);
      _test_triangulation_iterator(T2_0);
      _test_vertex_iterator(T1_0);
      _test_triangulation_iterator(T1_0);
      _test_vertex_iterator(T3_1);
      _test_triangulation_iterator(T3_1);
      _test_vertex_iterator(T3_0);
      _test_triangulation_iterator(T3_0); 
      _test_vertex_iterator(T3_2);        
      _test_triangulation_iterator(T3_2); 
      

      std::cout << "    Testing Circulator  "<< std::endl;
      _test_circulator(T0);
      _test_circulator(T3_1);
      _test_circulator(T3_0);
      _test_circulator(T3_2);
  }


  std::cout << "   Test is_Gabriel " << std::endl;
  Point q0(0.,0.,0.);
  Point q1(2.,0.,0.);
  Point q2(0.,2.,0.);
  Point q3(0.,0.,2.);
  Cls T4;
  v0 = T4.insert(q0);
  Vertex_handle v1 = T4.insert(q1);
  Vertex_handle v2 = T4.insert(q2, v1);         // testing with the hint
  Vertex_handle v3 = T4.insert(q3, v2->cell()); // testing with the hint
  Cell_handle c;
  int j,k,l;
  assert(T4.is_facet(v0,v1,v2,c,j,k,l));
  i = 6 - (j+k+l);
  Facet f = std::make_pair(c,i);
  assert(T4.is_Gabriel(c,i));
  assert(T4.is_Gabriel(f));
  assert(T4.is_facet(v1,v2,v3,c,j,k,l));
  i = 6 - (j+k+l);
  assert(!T4.is_Gabriel(c,i));
  assert(T4.is_edge(v0,v1,c,i,j));
  assert(T4.is_Gabriel(c,i,j));
  Edge e = make_triple(c,i,j);
  assert(T4.is_Gabriel(e));
  assert(T4.is_edge(v2,v3,c,i,j));
  assert(T4.is_Gabriel(c,i,j));
    
  std::cout <<"   Test dual (minimal test for now)" << std::endl;

  // We only test return types and instantiation, basically.
  {
    Cell_handle c = T4.finite_cells_begin();
    Point p = T4.dual(c);
    (void)p;
    Facet f = Facet(c, 2);
    CGAL::Object o = T4.dual(f);
    o = T4.dual(f.first, f.second);
  }

  std::cout <<"   Test destructors and quit "<< std::endl;

  T0.clear();
  assert(T0.is_valid());
  T3_0.clear();
  assert(T3_0.is_valid());
  assert(T3_0.dimension()==-1);
  assert(T3_0.number_of_vertices()==0);

  // "Determinism" test :
  // Triangulations built with the same order of input points
  // must have the same order of the vertex and cell iterator.
  {
    Cls Ta (q, q+22), Tb(q, q+22);
    assert(Ta == Tb);
    for (Finite_vertices_iterator ita = Ta.finite_vertices_begin(),
		                  itb = Tb.finite_vertices_begin(),
		                  end = Ta.finite_vertices_end();
	 ita != end; ++ita, ++itb)
      assert(ita->point() == itb->point());
    for (Finite_cells_iterator ita = Ta.finite_cells_begin(),
		               itb = Tb.finite_cells_begin(),
		               end = Ta.finite_cells_end();
	 ita != end; ++ita, ++itb) {
      assert(ita->vertex(0)->point() == itb->vertex(0)->point());
      assert(ita->vertex(1)->point() == itb->vertex(1)->point());
      assert(ita->vertex(2)->point() == itb->vertex(2)->point());
      assert(ita->vertex(3)->point() == itb->vertex(3)->point());
    }
  }
  
  /**********************/
  /******* MOVE *********/
  std::cout << "    displacements" << std::endl;

  std::cout << "    degenerate cases: " << std::endl;
  
  Cls TM_0;
  Vertex_handle tmv1 = TM_0.insert(Point(0,0,0));
  Vertex_handle tmv2 = TM_0.insert(Point(0,1,0));

	TM_0.move_if_no_collision(tmv1, Point(0, 2, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv1, Point(0, 0, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  Vertex_handle tmv3 = TM_0.insert(Point(0,2,1));

  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv3, Point(0, 2, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv3, Point(0, 2, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  Vertex_handle tmv4 = TM_0.insert(Point(0,1,1));
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv3, Point(1, 1, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv3, Point(4, 2, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv3, Point(0, 2, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(0, 2, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(0, 2, -1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(0, 3, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv3, Point(0, 1, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv3, Point(0, -1, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(0, -1, 0, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(0, -1, 0, 4));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(0, -1, 0, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(0, -1, 1, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(0, 0, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(0, 0, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(0, 0, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  assert(TM_0.move_if_no_collision(tmv1, Point(0, 3, 0)) != tmv1);

  TM_0.move_if_no_collision(tmv1, Point(0, 0, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(0, 1, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);   

  TM_0.move_if_no_collision(tmv4, Point(0, 3, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(0, 2, 3));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(0, 1, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  Vertex_handle tmv5 = TM_0.insert(Point(0,2,0));
  CGAL_USE(tmv5);
  Vertex_handle tmv6 = TM_0.insert(Point(1,0,0));
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv6, Point(0, 0, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv6, Point(2, 0, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv6, Point(2, 1, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv6, Point(0, 99, 99));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv6, Point(2, 1, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv6, Point(2, 2, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv6, Point(-2, 2, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  TM_0.move_if_no_collision(tmv6, Point(0, 1, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv6, Point(-2, 2, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 3);

  std::cout << "    random 1D: " << std::endl;
  Cls TM_1;
  // non-degenerate cases
  std::list<Point> points;
  for(int count=0; count<50; count++) {
    points.push_back(Point(0, 0, rand()%30000));
  }
  TM_1.insert(points.begin(), points.end());
  Vertex_handle vTM_1;
  for(int i=0; i<2; i++) {
    for(typename Cls::Finite_vertices_iterator 
          fvi = TM_1.finite_vertices_begin();
        fvi != TM_1.finite_vertices_end(); fvi++) {
      Point p = Point(0, 0, rand()%30000);
      vTM_1 = TM_1.move_if_no_collision(fvi, p);
      assert(TM_1.is_valid());
    }
  }
  assert(TM_1.is_valid());

  std::cout << "    random 2D: " << std::endl;
  Cls TM_2;
  // non-degenerate cases
  points.clear(); TM_2.clear();
  for(int count=0; count<10; count++) {
    points.push_back(Point(0, rand()%30000, rand()%30000));
  }
  TM_2.insert(points.begin(), points.end());
	Vertex_handle vTM_2;
  for(int i=0; i<2; i++) {
    for(typename Cls::Finite_vertices_iterator 
         fvi = TM_2.finite_vertices_begin();
         fvi != TM_2.finite_vertices_end(); fvi++) {
      Point p = Point(0, rand()%30000, rand()%30000);
      vTM_2 = TM_2.move_if_no_collision(fvi, p);
      assert(TM_2.is_valid());
    }
  }
  assert(TM_2.is_valid());

  std::cout << "    random 3D: " << std::endl;
  Cls TM_3;	
  // non-degenerate cases
  points.clear(); TM_3.clear();
  for(int count=0; count<50; count++) {
    points.push_back(Point(rand()%30000, rand()%30000, rand()%30000));
  }
  TM_3.insert(points.begin(), points.end());

  assert(TM_3.is_valid());
	
  Vertex_handle vTM_3;
  for(int i=0; i<2; i++) {
    for(typename Cls::Finite_vertices_iterator 
          fvi = TM_3.finite_vertices_begin();
        fvi != TM_3.finite_vertices_end(); fvi++) {
      Point p = Point(rand()%30000, rand()%30000, rand()%30000);
      vTM_3 = TM_3.move_if_no_collision(fvi, p);
      assert(TM_3.is_valid());
    }
  }

  // A simple test to see if move return the good vertex
  // when there is a collision
  assert(TM_3.move(TM_3.finite_vertices_begin(), vTM_3->point()) == vTM_3);

  // Test remove cluster
  {
		_test_remove_cluster<Triangulation>();
  }

}

#endif // CGAL_TEST_CLS_DELAUNAY_C
