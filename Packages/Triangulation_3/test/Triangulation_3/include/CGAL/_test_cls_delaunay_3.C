// Copyright (c) 1998-2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Francois Rebufat, Monique Teillaud, Sylvain Pion

#ifndef CGAL_TEST_CLS_DELAUNAY_C
#define CGAL_TEST_CLS_DELAUNAY_C

#include <cassert>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "_test_cls_iterator.C"
#include "_test_cls_circulator.C"

#include <CGAL/Random.h>

template <class Triangulation>
void
_test_cls_delaunay_3(const Triangulation &)
{
  typedef Triangulation                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested

  typedef typename Cls::Point                Point;
  typedef typename Cls::Segment              Segment;
  typedef typename Cls::Triangle             Triangle;
  typedef typename Cls::Tetrahedron          Tetrahedron;

  typedef typename Cls::Vertex               Vertex;
  typedef typename Cls::Cell                 Cell;
  typedef typename Cls::Facet                Facet;
  typedef typename Cls::Edge                 Edge;

  typedef typename Cls::Vertex_handle        Vertex_handle;
  typedef typename Cls::Cell_handle          Cell_handle; 
  typedef typename Cls::Vertex_iterator      Vertex_iterator;
  typedef typename Cls::Cell_iterator        Cell_iterator;
  typedef typename Cls::Locate_type          Locate_type;
  typedef std::list<Point>                        list_point;

  typedef typename Cls::Finite_vertices_iterator      Finite_vertices_iterator;
  typedef typename Cls::Finite_cells_iterator        Finite_cells_iterator;


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
      std::ofstream oFileT1("Test1_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT1("Test1_triangulation_IO_3",std::ios::in);
      oFileT1 << T0;
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
      std::ofstream oFileT2("Test2_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT2("Test2_triangulation_IO_3",std::ios::in);
      oFileT2 << T0;
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
      std::ofstream oFileT3("Test3_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT3("Test3_triangulation_IO_3",std::ios::in);
      oFileT3 << T0;
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
      std::ofstream oFileT4("Test4_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT4("Test4_triangulation_IO_3",std::ios::in);
      oFileT4 << T0;
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
      std::ofstream oFileT5("Test5_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT5("Test5_triangulation_IO_3",std::ios::in);
      oFileT5 << T0;
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
  int n = T1_0.insert(l1.begin(),l1.end());
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
      std::ofstream oFileT6("Test6_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT6("Test6_triangulation_IO_3",std::ios::in);
      oFileT6 << T1_2;
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
      std::ofstream oFileT7("Test7_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT7("Test7_triangulation_IO_3",std::ios::in);
      oFileT7 << T2_0;
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
	qq[m+20*n] = Point(m*px+n*qx, m*py+n*qy, 1);
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

    int n = Tdel.number_of_vertices();
    int m = Tdel.remove(vertices.begin(), vertices.end());
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
      std::ofstream oFileT8("Test8_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT8("Test8_triangulation_IO_3",std::ios::in);
      oFileT8 << T3_1;
      iFileT8 >> Tfromfile;
      assert(Tfromfile.is_valid());
      assert(Tfromfile.dimension() == 3);
      assert(Tfromfile.number_of_vertices() == 22);
  }

  // Testing find_conflicts() + insert_in_hole()
  std::cout << "    Constructor13 (find_conflicts/insert_in_hole)" << std::endl;
  Cls T3_13;
  for (i=0; i<22; ++i) {
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
      // Get the cells in conflicts.
      std::vector<Cell_handle> V;
      Facet facet;
      T3_13.find_conflicts(q[i], c, CGAL::Oneset_iterator<Facet>(facet),
                           std::back_inserter(V), CGAL::Emptyset_iterator());
      T3_13.insert_in_hole(q[i], V.begin(), V.end(), facet.first, facet.second);
    }
  }
  assert(T3_13.is_valid());
  assert(T3_13.number_of_vertices()==22);
  assert(T3_13.dimension()==3);

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
      int nbv = T.number_of_vertices();
      L.push_back(T.move_point(v, q[(3*i)%22]));

      if (nbv != T.number_of_vertices())
        L.pop_back(); // it means we move onto an already existing point.

      assert(T.is_valid());
      assert(T.number_of_vertices()<=22);
    }
  }

  {
      std::cout << "    Testing nearest_vertex()" << std::endl;
      // We do a nearest_vertex() query on all points with integer coordinate
      // in the cube [-1;6]^3, and we check explicitely that it's the nearest
      // by comparing to all vertices.
      for (int x = -1; x < 7; ++x)
        for (int y = -1; y < 7; ++y)
          for (int z = -1; z < 7; ++z) {
	    Point p(x, y, z);
	    Vertex_handle v = T3_13.nearest_vertex(p);
	    for (typename Cls::Finite_vertices_iterator
	         fvit = T3_13.finite_vertices_begin();
		 fvit != T3_13.finite_vertices_end(); ++fvit)
	      assert(squared_distance(p, v->point()) <=
		     squared_distance(p, fvit->point()));
	  }
  }

  {
      Cls Tfromfile;
      std::cout << "    I/O" << std::endl;
      std::ofstream oFileT8("Test13_triangulation_IO_3",std::ios::out);
      std::ifstream iFileT8("Test13_triangulation_IO_3",std::ios::in);
      oFileT8 << T3_13;
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

  std::cout <<"   Test destructors and quit "<< std::endl;

  T0.clear();
  assert(T0.is_valid());
  T3_0.clear();
  assert(T3_0.is_valid());
  assert(T3_0.dimension()==-1);
  assert(T3_0.number_of_vertices()==0);

}

#endif // CGAL_TEST_CLS_DELAUNAY_C
