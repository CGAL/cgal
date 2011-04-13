// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : include/CGAL/_test_cls_triangulation_3.C
// revision      : 
// revision_date : 

// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <cassert>

#include <iostream>
#include <fstream>

#include <CGAL/triple.h>
#include <list>
#include <vector>
#include "_test_cls_iterator.C"
#include "_test_cls_circulator.C"

#include <CGAL/Random.h>
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

  typedef typename Cls::Vertex_handle        Vertex_handle;
  typedef typename Cls::Cell_handle          Cell_handle; 
  typedef typename Cls::Vertex_iterator      Vertex_iterator;
  typedef typename Cls::Cell_iterator        Cell_iterator;
  //  typedef typename Cls::Point_iterator       Point_iterator;
  typedef typename Cls::Locate_type          Locate_type;
  typedef std::list<Point>                        list_point;


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

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
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

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
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

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
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

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
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

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
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

  // duplication
    
  T1.copy_triangulation(T0);
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


     //setting
     
  T1.set_number_of_vertices(100);
  assert(T1.dimension() == -1);
  assert(T1.number_of_vertices() == 100);
  //     assert(T1.is_valid());

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

  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
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


  //#######################################################################
  std::cout << "  list insertion"<< std::endl;
  Cls T3_2_1;
  T3_2_1.insert(lp2.begin(),lp2.end());
  assert(T3_2_1.is_valid());
  
  assert(T3_2_1.dimension()==3);
  //  assert(T3_2.number_of_vertices()==1000);
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

  Vertex_iterator vit;
  Vertex_handle w;
  for
    (vit=T3_1.finite_vertices_begin();vit!=T3_1.vertices_end();vit++)
    assert(T3_1.is_vertex(vit->point(), w));
         
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
     
  if (! del) { // Delaunay should not be flipped
    // or it will not be Delaunay any longer --> not valid
    std::cout << "  Test flip " << std::endl;
    assert( T3_1.is_valid());
    Cell_iterator cit, cdone = T3_1.cells_end();
    int nbflips=0;
    bool flipped;
    cit = T3_1.finite_cells_begin();
    Cell_iterator next_cell;
    while ( cit != cdone ) {
      // NOTE : cells are deleted during loop
      // the cell_iterator is modified "by hand" (not using ++)
      flipped = false; i=0; j=1;
      next_cell = ++cit; --cit;
      while ( (! flipped) && (i<4) ) {
	if ( (i!=j) ) {
	  flipped = T3_1.flip( &(*cit), i, j ) ;
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
	flipped = T3_1.flip( &(*cit), i );
	if (flipped) {
	  nbflips++;
	  assert(T3_1.is_valid());
	}
      }
    }
    std::cout << nbflips << " flips 2-3" << std::endl;
  }
       // Iterator and circulator test

  Cls T0_1;
  Point p28(1,3,5);
  v0=T0_1.insert(p28);
  if (! del) // to avoid doing the following tests for both Delaunay
    // and non Delaunay triangulations 
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
      _test_vertex_iterator(T3_3);
      _test_triangulation_iterator(T3_3); 
      

      std::cout << "    Testing Circulator  "<< std::endl;
      _test_circulator(T0);
      _test_circulator(T3_1);
      _test_circulator(T3_0);
      _test_circulator(T3_2);
      _test_circulator(T3_3);
    }

  std::cout <<"   Test destructors and quit "<< std::endl;

  T0.clear();
  assert(T0.is_valid());
  T3_0.clear();
  assert(T3_0.is_valid());
  assert(T3_0.dimension()==-1);
  assert(T3_0.number_of_vertices()==0);
        
       

    
}


