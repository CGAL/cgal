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
// file          : demo/Triangulation3/demo_color.C
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
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

template < class Traits, class Vb = CGAL::Triangulation_vertex_base_3<Traits> >
class My_vertex_base
  : public Vb
{
public :
  CGAL::Color color;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_base<Traits, Vb2>                    Other;
  };

  My_vertex_base() 
    : color(CGAL::WHITE) {}
};

struct K : public CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > {};

typedef CGAL::Triangulation_data_structure_3<My_vertex_base<K> >  Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                    Delaunay;

typedef Delaunay::Point Point;

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();

  Delaunay T;

  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));  
  T.insert(Point(0,1,0));  
  T.insert(Point(0,0,1));  
  T.insert(Point(2,2,2));  
  T.insert(Point(-1,0,1));  

  // Set the color of finite vertices of degree 6 to red.
  Delaunay::Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    if (T.degree(vit) == 6)
      vit->color = CGAL::RED;

  std::cout << "           Visualization of T" << std::endl;
  gv.set_wired(true);
  gv << T;

  std::cout << "           Vertices of T with their own color" << std::endl
	    << "           red for degree 6 (counting infinite vertex)" 
	    << std::endl 
	    << "           white otherwise" << std::endl;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    gv << vit->color << vit->point();

  char ch;
  std::cout << "Enter any character to quit" << std::endl;
  std::cin >> ch;

  return 0;
}

#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
