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
// file          : demo/Triangulation3/demo_remove.C
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

#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#include <unistd.h>
#include <vector>

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > K;

typedef CGAL::Triangulation_vertex_base_3<K>            Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
typedef CGAL::Triangulation_cell_base_3<void>            Cb;
typedef CGAL::Triangulation_data_structure_3<Vbh,Cb>     Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>           Dt;
typedef CGAL::Triangulation_hierarchy_3<Dt>              Dh;

typedef Dh::Vertex_iterator Vertex_iterator;
typedef Dh::Vertex_handle Vertex_handle;
typedef Dh::Cell_handle Cell_handle;
typedef K::Point_3 Point;

////////////////////// 
// VISU GEOMVIEW
////////////////////// 
template<class TRIANGULATION>
void visu_cell(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c)
{
  if ( ! T.is_infinite(c) )
    os << T.tetrahedron(c);
  else 
    os << T.triangle(c,c->index(T.infinite_vertex()));
}
template<class TRIANGULATION>
void visu_facet(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c, int i)
{
  if ( ! T.is_infinite(c,i) )
    os << T.triangle(c,i);
}
template<class TRIANGULATION>
void visu_edge(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c, int i, int j)
{
  if ( ! T.is_infinite(c,i,j) )
    os << T.segment(c,i,j);
}
template<class TRIANGULATION>
void visu_vertices(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Vertex_iterator vit = T.finite_vertices_begin();
  Vertex_iterator vdone = T.vertices_end();
  
  if ( vit == vdone ) { std::cout << "no vertex" << std::endl ;}
  else {
    while(vit != vdone) {
      os << vit->point();
      ++vit;
    }
  }
}
template<class TRIANGULATION>
void visu_vertex(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c, int i)
{
  if ( ! T.is_infinite(c->vertex(i)) )
    os << c->vertex(i)->point();
}

////////////////////// 

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 5, 5, 5));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.set_wired(true);
  gv.clear();

  Dh T;

  std::cout <<"          Inserting points" << std::endl ;
  int x,y,z;
  std::vector<Vertex_handle> V(125);
  int i=0;
 
  for (z=0 ; z<5 ; z++)
    for (y=0 ; y<5 ; y++)
      for (x=0 ; x<5 ; x++) 
	  V[i++] = T.insert(Point(x,y,z));

  assert( T.is_valid() );
  assert( T.number_of_vertices() == 125 );
  assert( T.dimension() == 3 );

  std::cout <<"          Visualizing edges" << std::endl;
  gv << T;

  sleep(3);

  std::cout <<"          Removing vertices in random order" << std::endl;
  
  std::random_shuffle(V.begin(), V.end());

  for (i=0; i<125; ++i) {
    assert( T.remove(V[i]) );
    gv.clear();
    gv << T;
  }

  char ch;
  std::cout << "Enter any character to quit" << std::endl;
  std::cin >> ch;

  return 1;
}

#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
