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
// file          : test_regular_3.C
// revision      : 
// revision_date : 
// author(s)     : Monique Teillaud (Monique.Teillaud@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>
#include <list>
#include <vector>

#include <CGAL/_test_types.h>
#include <CGAL/triple.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer my_NT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz my_NT;
#else
#include <CGAL/double.h>
typedef double my_NT;
#endif
#endif
// #endif

#include <CGAL/Cartesian.h>

typedef CGAL::Cartesian<my_NT> Test_rep_cartesian;
typedef CGAL::Homogeneous<my_NT> Test_rep_homogeneous;

bool del=true;

int main()
{

  std::cout << " with CGAL::Regular_triangulation_euclidean_traits_3: " << std::endl;
 
  typedef CGAL::Regular_triangulation_euclidean_traits_3<Test_rep_cartesian>  traits;
// works with both geom_traits
 // typedef CGAL::_Triangulation_test_traits_3                         traits;
  typedef CGAL::Triangulation_vertex_base_3<traits>                 Vb;
  typedef CGAL::Triangulation_cell_base_3<traits>                   Fb;
  typedef CGAL::Triangulation_data_structure_3<Vb,Fb>               Tds;
  typedef CGAL::Regular_triangulation_3<traits,Tds>                Cls;

  typedef traits::Bare_point Point;
  typedef traits::Weighted_point Weighted_point;

  typedef list<Weighted_point>                        list_point;

  //  _test_cls_triangulation_3( Cls() );

  // temporary version
  Cls T;

  // 4 first points to ensure that dimension is non degenerate
  T.insert( Weighted_point(Point(0,0,0),1) );
  T.insert( Weighted_point(Point(1,0,0),1) );
  T.insert( Weighted_point(Point(0,1,0),1) );
  T.insert( Weighted_point(Point(0,0,1),1) );


  list_point lp;
  int a, b, d;
  for (a=0;a!=10;a++)
    for (b=0;b!=10;b++)
      for (d=0;d!=10;d++)
	lp.push_back(Weighted_point( Point(a*b-d*a + (a-b)*10 +a ,a-b+d +5*b,
					   a*a-d*d+b),
				     a*b-a*d) );
  list_point::iterator it;
  int count = 0 ;
  std::cout << " number of inserted points : " ;
  for (it=lp.begin(); it!=lp.end();it++){
    count++;
    T.insert(*it);
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

  std::cout << " number of vertices : " << T.number_of_vertices() << std::endl;

  assert(T.is_valid(true));
  assert(T.dimension()==3);

  return 0;
}
