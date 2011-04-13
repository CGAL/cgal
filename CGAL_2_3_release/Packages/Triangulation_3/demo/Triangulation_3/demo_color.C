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

template < class Traits >
class My_vertex_base
  : public CGAL::Triangulation_vertex_base_3<Traits>
{
public :
  CGAL::Color color;
  typedef typename Traits::Point_3 Point;

  My_vertex_base() 
    : CGAL::Triangulation_vertex_base_3<Traits>(), color(CGAL::WHITE)
    {}

  My_vertex_base(const Point & p) 
    : CGAL::Triangulation_vertex_base_3<Traits>(p), color(CGAL::WHITE)
    {}

  My_vertex_base(const Point & p, void* c) 
    : CGAL::Triangulation_vertex_base_3<Traits>(p,c), color(CGAL::WHITE)
    {} 

  My_vertex_base(void* c)
    : CGAL::Triangulation_vertex_base_3<Traits>(c), color(CGAL::WHITE)
    {} 
};

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > K;

typedef K::Point_3 Point;

typedef CGAL::Triangulation_cell_base_3<K> Cb;
typedef My_vertex_base<K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;

typedef Delaunay::Vertex_iterator Vertex_iterator;
typedef Delaunay::Vertex_handle Vertex_handle;

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

  Vertex_iterator vit;
  std::set<Vertex_handle> adjacent;
  for (vit = T.finite_vertices_begin(); vit != T.vertices_end(); ++vit) {
    T.incident_vertices( &*vit, adjacent);
    if (adjacent.size() == 6) 
      vit->color = CGAL::RED;
  }

  std::cout << "           Visualization of T" << std::endl;
  gv.set_wired(true);
  gv << T;

  std::cout << "           Vertices of T with their own color" << std::endl
	    << "           red for degree 6 (counting infinite vertex)" 
	    << std::endl 
	    << "           white otherwise" << std::endl;
  for (vit = T.finite_vertices_begin(); vit != T.vertices_end(); ++vit) {
    gv << vit->color;
    gv << vit->point();
  }

  char ch;
  std::cout << "Enter any character to quit" << std::endl;
  std::cin >> ch;

  return 1;
}

#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
