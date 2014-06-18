// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
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

#include <cassert>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "_test_cls_iterator.h"
#include "_test_cls_circulator.h"

#include <CGAL/Random.h>
#include <CGAL/Testsuite/use.h>
#include <CGAL/use.h>

template <class Triangulation, class Container>
bool check_all_are_finite(Triangulation* tr, const Container& cont) 
{
  for(typename Container::const_iterator it = cont.begin(), end = cont.end();
      it != end; ++it)
  {
    if(tr->is_infinite(*it)) return false;
  }
  return true;
}

template <class Triangulation>
void
_test_cls_triangulation_3_input_output(const Triangulation & T,
				       const char* filename)
{
  const int dim = T.dimension();
  const typename Triangulation::size_type n = T.number_of_vertices();
  std::cout << "    I/O" << std::endl;
  {
    std::ofstream oFile(filename, std::ios::out);
    oFile << T << std::endl;
  }
  std::ifstream iFile(filename, std::ios::in);
  Triangulation Tfromfile;
  iFile >> Tfromfile;
  assert(Tfromfile.is_valid());
  assert(Tfromfile.dimension() == dim);
  assert(Tfromfile.number_of_vertices() == n);


  std::string filename_bin = filename;
  filename_bin = filename_bin + "_binary";

  const char* filename2 = filename_bin.c_str();
  std::cout << "    I/O (binary)" << std::endl;
  {
    std::ofstream oFileBin(filename2, std::ios::out|std::ios::binary);
    CGAL::set_binary_mode(oFileBin);
    oFileBin << T;
  }
  std::ifstream iFileBin(filename2, std::ios::in|std::ios::binary);
  CGAL::set_binary_mode(iFileBin);
  Triangulation Tfromfile_binary;
  iFileBin >> Tfromfile_binary;
  assert(Tfromfile_binary.is_valid());
  assert(Tfromfile_binary.dimension() == dim);
  assert(Tfromfile_binary.number_of_vertices() == n);
}

template <class Triangulation>
void
_test_cls_triangulation_3(const Triangulation &)
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

  typedef typename Cls::size_type            size_type;
  typedef typename Cls::difference_type      difference_type;

  typedef typename Cls::Vertex_handle        Vertex_handle;
  typedef typename Cls::Cell_handle          Cell_handle; 
  typedef typename Cls::Vertex_iterator      Vertex_iterator;
  typedef typename Cls::Cell_iterator        Cell_iterator;
  //  typedef typename Cls::Point_iterator       Point_iterator;
  typedef typename Cls::Locate_type          Locate_type;
  typedef std::list<Point>                        list_point;

  typedef typename Cls::Finite_vertices_iterator    Finite_vertices_iterator;
  typedef typename Cls::Finite_edges_iterator       Finite_edges_iterator;
  typedef typename Cls::Finite_facets_iterator      Finite_facets_iterator;
  typedef typename Cls::Finite_cells_iterator       Finite_cells_iterator;

  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Cell);
  CGAL_USE_TYPE(difference_type);
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
  // Beginning with an empty triangulation and adding points until reaching
  // 3-dimensional triangulation.
  Cls T0; 
  assert(T0.dimension() == -1);
  assert(T0.number_of_vertices() == 0);
  assert(T0.is_valid());

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T0, "Test1_triangulation_IO_3");
    }

  std::cout << "    Constructor1 " << std::endl;
  Point p10(0,0,0);
  Vertex_handle v0=T0.insert(p10);
  assert(T0.dimension() == 0);
  assert(T0.number_of_vertices() == 1);
  assert(T0.is_valid());

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T0, "Test2_triangulation_IO_3");
    }

  std::cout << "    Constructor2 " << std::endl;

  Point p11(100,100,0);
  v0=T0.insert(p11);
  assert(T0.dimension() == 1);
  assert(T0.number_of_vertices() == 2);
  assert(T0.is_valid());

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T0, "Test3_triangulation_IO_3");
    }

  std::cout << "    Constructor3 " << std::endl;

  Point p12(100,-100,0);
  v0=T0.insert(p12);
  assert(T0.dimension() == 2);
  assert(T0.number_of_vertices() == 3);
  assert(T0.is_valid());

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T0, "Test4_triangulation_IO_3");
    }

  std::cout << "    Constructor4 " << std::endl;

  Point p13(50,0,100);
  v0=T0.insert(p13);
  assert(T0.dimension() == 3);
  assert(T0.number_of_vertices() == 4);
  assert(T0.is_valid());

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T0, "Test5_triangulation_IO_3");
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

  std::cout << "    Testing operator==" << std::endl;
  assert(T0 == T0);
  assert(T0 == T1);
  assert(T1 == T0);
  assert(T1 == T1);
   
  T1.clear();
  assert(T1.dimension() == -1);
  assert(T1.number_of_vertices() == 0);
  assert(T1.is_valid());



   // Assignment
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

     // Building some psychotic triangulations :
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
//   // test for point_iterator
//   Point_iterator pit;
//   Point pppp;
//   int nbpt=0;
//   for (pit = T1_2.points_begin(); pit != T1_2.points_end(); ++pit) {
//     nbpt++;
//     pppp = *pit;
//   }
//   assert(nbpt==n);
  assert(T1_2.is_valid());

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T1_2, "Test6_triangulation_IO_3");
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

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T2_0, "Test7_triangulation_IO_3");
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


   std::cout << "    Constructor12 " << std::endl;
   Cls T3_1;
  for (i=0;i<22;i++)
    T3_1.insert(q[i]);
  assert(T3_1.is_valid());
  assert(T3_1.number_of_vertices()==22);
  assert(T3_1.dimension()==3);

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
    {
      _test_cls_triangulation_3_input_output(T3_1, "Test8_triangulation_IO_3");
    }


  //#######################################################################
  std::cout << "  list insertion"<< std::endl;
  Cls T3_2_1;
  T3_2_1.insert(lp2.begin(),lp2.end());
  assert(T3_2_1.is_valid());

  assert(T3_2_1.dimension()==3);
  //  assert(T3_2.number_of_vertices()==1000);
  std::cout << "   end of insertion " << std::endl;

  // The same using the template constructor.
  std::cout << "  template ctor"<< std::endl;
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
  v0=T3_5.insert(p118, v0->cell()); // testing with the hint
  v0=T3_5.insert(p119, v0);         // testing with the hint
  
  assert(T3_5.is_valid());
  assert(T3_5.number_of_vertices()==10);

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


  // Test inserts function separatelly.

  std::cout << "    Testing insertions   " << std::endl;
  Locate_type lt;
  int li,lj,i1,i2;
  Cls Ti = T0;
  Point p20(50,0,50);
  v0=Ti.insert_in_cell(p20,Ti.locate(Point(50,0,50)));
  assert(Ti.is_valid());
  assert(Ti.number_of_vertices() == 5);
  Ti=T0;
  Point p120(50,0,0);
  v0=Ti.insert_in_facet(p120,Ti.locate(Point(50,0,1)),3);
  assert(Ti.is_valid());
  assert(Ti.number_of_vertices() == 5);
  Ti=T0;
  Point p21(50,0,0);
  v0=Ti.insert_in_facet(p21,Facet(Ti.locate(Point(50,0,1)),3));
  assert(Ti.is_valid());
  assert(Ti.number_of_vertices() == 5);
  Ti=T0;
  Cell_handle c= Ti.locate(Point(0,0,0),lt,li,lj);
  assert(lt==Cls::VERTEX);
  i1=li;
  c= Ti.locate(Point(100,100,0),lt,li,lj);
  assert(lt==Cls::VERTEX);
  i2=li;

  // testing the locate with hint.
  Cell_handle ccc = Ti.locate(Point(100,100,0),c);
  assert(c == ccc);
  ccc = Ti.locate(Point(100,100,0),c->vertex(0));

  Point p22(50,50,0);
  v0=Ti.insert_in_edge(p22,Ti.locate(Point(50,40,1)),i1,i2);
  assert(Ti.is_valid());
  assert(Ti.number_of_vertices() == 5);

  Ti=T0;
  c= Ti.locate(Point(0,0,0),lt,li,lj);
  assert(lt==Cls::VERTEX);
  i1=li;
  c= Ti.locate(Point(100,100,0),lt,li,lj);
  assert(lt==Cls::VERTEX);
  i2=li;
  Point p23(50,50,0);
  v0=Ti.insert_in_edge(p23,Edge(Ti.locate(Point(50,50,0)),i1,i2));
  assert(Ti.is_valid());
  assert(Ti.number_of_vertices() == 5);

  Ti=T0;

  assert(T0.dimension() == 3);
  assert(T0.number_of_vertices() == 4);
  assert(T0.is_valid());
       

  c= Ti.locate(Point(50,50,50),lt,li,lj);
      
  Point p24(50,50,50);
  v0= Ti.insert_outside_convex_hull(p24,c);
  assert(Ti.is_valid());

  assert(Ti.number_of_vertices() == 5);

  Cls T3_3=T1_0;
  Point p25(2,0,0);
  v0=T3_3.insert_outside_affine_hull(p25);
  assert(T3_3.is_valid());
  assert(T3_3.dimension()==2);
  c= T3_3.locate(Point(4,0,0),lt,li,lj);
  Point p26(4,0,0);
  v0=T3_3.insert_outside_convex_hull(p26,c);
  assert(T3_3.is_valid());
  assert(T3_3.dimension()==2);
  Point p27(0,5,0);
  v0=T3_3.insert_outside_affine_hull(p27);
  assert(T3_3.is_valid());
  assert(T3_3.dimension()==3);

  // ################## Operations + newly created cells ################
  // Small test for inserting and returning the newly created cells
  // (the code is mainly the usual insert + incident_{edges,facets,cells} 
  // depending on the dimension)

  std::cout << " Test insertion + newly created cells: " << std::endl;

  std::cout << " 1 dimension" << std::endl;
  // dimension 1
  Cls TAI1;
  for(int i=0; i<50; i++)
  {	
    double x = (double) (2*i);
    TAI1.insert(Point(x, x, x));
  }
  std::list<Cell_handle> lis_tai1;
  for(int i=0; i<51; i++)
  {
    lis_tai1.clear();
    double x = (double) (2*i - 1);
    Vertex_handle taiv = 
      TAI1.insert_and_give_new_cells(
                                     Point(x, x, x), 
                                     std::back_inserter(lis_tai1));
    CGAL_USE(taiv);
    assert(TAI1.is_valid());
    assert(TAI1.dimension() == 1);
    assert(lis_tai1.size() == 2);
    while(!lis_tai1.empty())
    {
      Cell_handle c = lis_tai1.front();
      lis_tai1.pop_front();
      assert(TAI1.tds().is_simplex(c));
    }	
  }
  TAI1.clear();

  std::cout << " 2 dimensions" << std::endl;	
  CGAL::Random grand;
  for(int i=0; i<50; i++)
  {	
    double x = grand.get_double();
    double y = grand.get_double();		  
    TAI1.insert(Point(x, y, 0));
  }
  for(int i=0; i<50; i++)
  {
    lis_tai1.clear();
    double x = grand.get_double();
    double y = grand.get_double();
    Vertex_handle taiv = 
      TAI1.insert_and_give_new_cells(
                                     Point(x, y, 0), 
                                     std::back_inserter(lis_tai1));
    CGAL_USE(taiv);
    assert(TAI1.is_valid());
    assert(TAI1.dimension() == 2);
    while(!lis_tai1.empty())
    {
      Cell_handle c = lis_tai1.front();
      lis_tai1.pop_front();
      assert(TAI1.tds().is_simplex(c));
    }	
  } 
  TAI1.clear();

  std::cout << " 3 dimensions" << std::endl;	
  for(int i=0; i<50; i++)
  {	
    double x = grand.get_double();
    double y = grand.get_double();	
    double z = grand.get_double();			  
    TAI1.insert(Point(x, y, z));
  }
  for(int i=0; i<50; i++)
  {
    lis_tai1.clear();
    double x = grand.get_double();
    double y = grand.get_double();
    double z = grand.get_double();	
    Vertex_handle taiv = 
      TAI1.insert_and_give_new_cells(
                                     Point(x, y, z), 
                                     std::back_inserter(lis_tai1));
    CGAL_USE(taiv);
    assert(TAI1.is_valid());
    assert(TAI1.dimension() == 3);
    while(!lis_tai1.empty())
    {
      Cell_handle c = lis_tai1.front();
      lis_tai1.pop_front();
      assert(TAI1.tds().is_simplex(c));
    }	
  } 
  TAI1.clear(); 

  // the other two insertion methods is exactly the same
  // with different version of the basic insert method
  // Vertex_handle insert_and_give_new_cells(const Point& p,
  //       OutputItCells fit,
  //       Vertex_handle hint)	
  // Vertex_handle insert_and_give_new_cells(const Point& p,
  //       Locate_type lt, Cell_handle c, int li, int lj, 
  //       OutputItCells fit
	
	
  // ##################################################################


       // testing some simple basic methods (access functions)

  std::cout << "   Boolean and query functions " <<std::endl;
  c=T0.infinite_cell();
  assert(T0.is_infinite(c));
  int ind=c->index(T0.infinite_vertex());
       
  Facet f ; 
  for (i=0;i<4;i++) 
    if (i!=ind) {
      assert(T0.is_infinite(c,i));
      f=Facet(c,i);
      assert(T0.is_infinite(f));
    }
  int j;

  for (i=0;i<4;i++) 
    for (j=0;i<4;i++) 
      if ((i!=j) && ((i==ind) || (j==ind))) {
	assert(T0.is_infinite(c,i,j));
	assert(T0.is_infinite(Edge(c,i,j)));
      }


  v0=T0.infinite_vertex();
  assert(T0.is_infinite(v0));

  Finite_vertices_iterator vit;
  Vertex_handle w;
  for (vit=T3_1.finite_vertices_begin();vit!=T3_1.finite_vertices_end();vit++)
    assert(T3_1.is_vertex(vit->point(), w));

  // test mirrors
  c = T0.infinite_cell();
  for (i=0;i<4;i++) {
    Cell_handle d = c->neighbor(i);
    int j  = T0.mirror_index(c,i);
    assert(d->vertex(j) == T0.mirror_vertex(c,i));
    assert(Facet(d,j) == T0.mirror_facet(Facet(c,i)));
  }
	   
    
         
       // geometric functions
  std::cout << "Geometric functions " << std::endl;
  c= T0.locate(Point(50,0,1),lt,li,lj);
  Tetrahedron tr1=T0.tetrahedron(c);
  c= T0.locate(Point(10,0,1),lt,li,lj);
  Tetrahedron tr2=T0.tetrahedron(c);
  assert(tr1==tr2);
  c= T0.locate(Point(50,0,5),lt,li,lj);
  Triangle tri1=T0.triangle(c,1);
  c= T0.locate(Point(10,0,1),lt,li,lj);
  Triangle tri2=T0.triangle(Facet(c,1));
  assert(tri1==tri2);
  c= T0.locate(Point(10,0,1),lt,li,lj);
  Segment s1 = T0.segment(c,0,1);
  c= T0.locate(Point(50,0,5),lt,li,lj);
  Segment s2 = T0.segment(Edge(c,0,1));
  assert(s1==s2);
  c= T0.locate(Point(50,0,5),lt,li,lj);
  Point pt1 = T0.point(c,0);
  c= T0.locate(Point(10,0,1),lt,li,lj);
  Point pt2 = T0.point(c->vertex(0));
  assert(pt1==pt2);
  c= T0.locate(Point(20,0,2),lt,li,lj);
  Point pt3 = c->vertex(0)->point();
  assert(pt2==pt3);
     
  if (! del) { // Delaunay should not be flipped
    // or it will not be Delaunay any longer --> not valid
    std::cout << "  Test flip " << std::endl;
    assert( T3_1.is_valid());
    Finite_cells_iterator cit, cdone = T3_1.finite_cells_end();
    int nbflips=0;
    bool flipped;
    cit = T3_1.finite_cells_begin();
    Finite_cells_iterator next_cell;
    while ( cit != cdone ) {
      // NOTE : cells are deleted during loop
      // the cell_iterator is modified "by hand" (not using ++)
      flipped = false; i=0; j=1;
      next_cell = ++cit; --cit;
      while ( (! flipped) && (i<4) ) {
	if ( (i!=j) ) {
	  flipped = T3_1.flip( cit, i, j ) ;
	  if (flipped) {
	    nbflips++;
	    assert(T3_1.is_valid());
	  }
	}
	if ( j==3 ) { i++; j=0; }
	else j++;
      }
      cit = next_cell;
    }
    std::cout << nbflips << " flips 3-2" << std::endl;

    nbflips=0;
    for ( cit = T3_1.finite_cells_begin(); cit != cdone; cit++ ) {
      // NOTE : the triangulation is modified during loop
      // --> the cell_iterator does not mean a lot
      for ( i=0; i<4; i++ ) {
	flipped = T3_1.flip( cit, i );
	if (flipped) {
	  nbflips++;
	  assert(T3_1.is_valid());
	}
      }
    }
    std::cout << nbflips << " flips 2-3" << std::endl;
  }
  
  // Finite incident_* in dimension 2 test
  std::cout << "    Testing finite_incident_* in dim 2  "<< std::endl;
  Cls* T2[2];
  T2[0] = &T2_0;
  T2[1] = &T2_1;
  
    for(int k = 0; k < 2; ++k) {
    std::cout << "      with triangulation " << k + 1 << ": ";

    std::vector<Vertex_handle> f_vertices_old;
    std::vector<Vertex_handle> f_vertices;
    std::vector<Edge> f_edges;
    std::vector<Facet> f_facets;
    std::vector<Cell_handle> f_cells;

    f_vertices.clear();
    f_edges.clear();
    f_facets.clear();
    f_cells.clear();
    
    for(Finite_vertices_iterator i = T2[k]->finite_vertices_begin();
	i != T2[k]->finite_vertices_end(); ++i) {
      // old name (up to CGAL 3.4)
      // kept for backwards compatibility but not documented
      T2[k]->finite_incident_vertices(i, std::back_inserter(f_vertices_old));
      // correct name 
      T2[k]->finite_adjacent_vertices(i, std::back_inserter(f_vertices));
      assert(check_all_are_finite(T2[k], f_vertices));
      T2[k]->finite_incident_edges(i, std::back_inserter(f_edges));
      assert(check_all_are_finite(T2[k], f_edges));
      T2[k]->finite_incident_facets(i, std::back_inserter(f_facets));
      assert(check_all_are_finite(T2[k], f_facets));
      T2[k]->finite_incident_cells(i, std::back_inserter(f_cells));
      if(T2[k]->dimension() == 3) { assert(check_all_are_finite(T2[k], f_cells)); }
    }
    unsigned int nb_f_edges = 0;
    Finite_edges_iterator feit = T2[k]->finite_edges_begin();
    while(feit != T2[k]->finite_edges_end()) {
      assert(!T2[k]->is_infinite(*feit));
      ++nb_f_edges;
      ++feit;
    }
    unsigned int nb_f_facets = 0;
    Finite_facets_iterator ffait = T2[k]->finite_facets_begin();
    while(ffait != T2[k]->finite_facets_end()) {
      ++nb_f_facets;
      ++ffait;
    }

    // incidences
    assert(f_edges.size() == f_vertices_old.size());
    assert(f_edges.size() == f_vertices.size());
    assert(2*nb_f_edges == f_edges.size());
    assert(3*nb_f_facets == f_facets.size());
    assert(3*nb_f_facets == f_cells.size());
    
    typename Cls::size_type nb_f_vertices = T2[k]->number_of_vertices();
    
    // Euler relation
    assert(nb_f_vertices - nb_f_edges + nb_f_facets == 1);
    std::cout << "ok\n";
  }
  
       // Finite incident_* to vertex test
  std::cout << "    Testing finite_incident_* in dim 3  "<< std::endl;

  Cls* T3[6];
  T3[0] = &T3_0;
  T3[1] = &T3_1;
  T3[2] = &T3_2_1;
  T3[3] = &T3_2_2;
  T3[4] = &T3_2;
  T3[5] = &T3_3;

  for(int k = 0; k < 6; ++k) {
    std::cout << "      with triangulation " << k + 1 << ": ";

    std::vector<Vertex_handle> f_vertices_old;
    std::vector<Vertex_handle> f_vertices;
    std::vector<Edge> f_edges;
    std::vector<Facet> f_facets;
    std::vector<Cell_handle> f_cells;

    f_vertices.clear();
    f_edges.clear();
    f_facets.clear();
    f_cells.clear();
    
    for(Finite_vertices_iterator i = T3[k]->finite_vertices_begin();
	i != T3[k]->finite_vertices_end(); ++i) {
      // old name (up to CGAL 3.4)
      // kept for backwards compatibility but not documented
      T3[k]->finite_incident_vertices(i, std::back_inserter(f_vertices_old));
      // correct name 
      T3[k]->finite_adjacent_vertices(i, std::back_inserter(f_vertices));
      assert(check_all_are_finite(T3[k], f_vertices));
      T3[k]->finite_incident_edges(i, std::back_inserter(f_edges));
      assert(check_all_are_finite(T3[k], f_edges));
      T3[k]->finite_incident_facets(i, std::back_inserter(f_facets));
      assert(check_all_are_finite(T3[k], f_facets));
      T3[k]->finite_incident_cells(i, std::back_inserter(f_cells));
      if(T3[k]->dimension()==3) { assert(check_all_are_finite(T3[k], f_cells)); }
    }
    unsigned int nb_f_edges = 0;
    Finite_edges_iterator feit = T3[k]->finite_edges_begin();
    while(feit != T3[k]->finite_edges_end()) {
      assert(!T3[k]->is_infinite(*feit));
      ++nb_f_edges;
      ++feit;
    }
    unsigned int nb_f_facets = 0;
    Finite_facets_iterator ffait = T3[k]->finite_facets_begin();
    while(ffait != T3[k]->finite_facets_end()) {
      ++nb_f_facets;
      ++ffait;
    }
    unsigned int nb_f_cells = 0;
    Finite_cells_iterator fcit = T3[k]->finite_cells_begin();
    while(fcit != T3[k]->finite_cells_end()) {
      ++nb_f_cells;
      ++fcit;
    }
    
    // incidences
    assert(f_edges.size() == f_vertices_old.size());
    assert(f_edges.size() == f_vertices.size());
    assert(2*nb_f_edges == f_edges.size());
    assert(3*nb_f_facets == f_facets.size());
    assert(4*nb_f_cells == f_cells.size());
    
    typename Cls::size_type nb_f_vertices = T3[k]->number_of_vertices();
    
    // Euler relation
    assert(nb_f_vertices - nb_f_edges + nb_f_facets - nb_f_cells == 1);
    std::cout << "ok\n";
  }
       // Iterator and circulator test

  Cls T0_1;
  Point p28(1,3,5);
  v0=T0_1.insert(p28);

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
  _test_vertex_iterator(T3_3);
  _test_triangulation_iterator(T3_3); 

  std::cout << "    Testing Circulator  "<< std::endl;
  _test_circulator(T0);
  _test_circulator(T3_1);
  _test_circulator(T3_0);
  _test_circulator(T3_2);
  _test_circulator(T3_3);

  std::cout <<"   Test destructors and quit "<< std::endl;

  T0.clear();
  assert(T0.is_valid());
  T3_0.clear();
  assert(T3_0.is_valid());
  assert(T3_0.dimension()==-1);
  assert(T3_0.number_of_vertices()==0);
}
