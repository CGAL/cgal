// ============================================================================
//
// Copyright (c) 1998-1999 The CGAL Consortium
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
// file          : demo/Triangulation3/demo.C
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================

// Geomview doesn't work on M$ at the moment, so we don't compile this file.
#if defined(__BORLANDC__) || defined(_MSC_VER)
#include <iostream>
int main()
{
  std::cerr << "Geomview doesn't work on Windows, so this demo doesn't work"
            << std::endl;
  return 0;
}
#else

#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h> 
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Quotient.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#include <iostream>

// an exact number type is needed because we are using constructions
// (circumcenter computations) in this demo, not only predicates 
typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > NT;

typedef CGAL::Simple_cartesian<NT> K;

typedef CGAL::Delaunay_triangulation_3<K> Triangulation;

typedef Triangulation::Point          Point;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 3, 3, 3));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();

  Triangulation T;

  int x,y,z;
  std::cout <<"          Inserting points" << std::endl ;
  for (z=0 ; z<3 ; z++)
    for (y=0 ; y<3 ; y++)
      for (x=0 ; x<3 ; x++) 
	  T.insert(Point(NT(x),NT(y),NT(z)));

  T.is_valid(true);

  std::cout <<"          Visualizing T" << std::endl;
  gv.set_wired(true);
  gv << T;

  std::cout <<"          Visualizing the Voronoi edges" << std::endl;
  gv << CGAL::RED;
  T.draw_dual(gv);

  char ch;
  std::cout << "Enter any character to quit" << std::endl;
  std::cin >> ch;

  return 0;
}

#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
