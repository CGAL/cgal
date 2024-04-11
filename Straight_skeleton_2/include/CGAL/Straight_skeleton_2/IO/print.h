// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Fernando Cacciola

#ifndef CGAL_SLS_IO_PRINT_H
#define CGAL_SLS_IO_PRINT_H

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_2.h>

namespace CGAL {
namespace Straight_skeletons_2 {
namespace IO {

template<class K>
void print_point ( CGAL::Point_2<K> const& p )
{
  std::cout << "(" << p.x() << "," << p.y() << ")" ;
}

template<class HDS_V>
void print_vertex ( HDS_V const& v )
{
  if(v->has_infinite_time())
    std::cout << "F" << v->id() << " " ;
  else
    std::cout << "N" << v->id() << " " ;
}

template<class K, class C>
void print_polygon ( CGAL::Polygon_2<K,C> const& poly )
{
  typedef CGAL::Polygon_2<K,C> Polygon ;

  std::cout << "Polygon with " << poly.size() << " vertices" << std::endl ;

  for( typename Polygon::Vertex_const_iterator vi = poly.vertices_begin() ; vi != poly.vertices_end() ; ++ vi )
  {
    print_point(*vi); std::cout << std::endl ;
  }
}

template<class K, class C>
void print_polygons ( std::vector< std::shared_ptr< CGAL::Polygon_2<K,C> > > const& polies )
{
  typedef std::vector< std::shared_ptr< CGAL::Polygon_2<K,C> > > PolygonVector ;

  std::cout << "Polygon list with " << polies.size() << " polygons" << std::endl ;

  for( typename PolygonVector::const_iterator pi = polies.begin() ; pi != polies.end() ; ++ pi )
    print_polygon(**pi);
}

template<class K, class C>
void print_polygon_with_holes ( CGAL::Polygon_with_holes_2<K,C> const& polywh )
{
  typedef CGAL::Polygon_with_holes_2<K,C> PolygonWithHoles ;

  std::cout << "Polygon_with_holes having " << polywh.number_of_holes() << " holes" << std::endl ;

  print_polygon(polywh.outer_boundary());

  for( typename PolygonWithHoles::Hole_const_iterator hi = polywh.holes_begin() ; hi != polywh.holes_end() ; ++ hi )
    print_polygon(*hi);
}

template<class K, class C>
void print_polygons_with_holes ( std::vector< std::shared_ptr< CGAL::Polygon_with_holes_2<K,C> > > const& polies )
{

  typedef std::vector< std::shared_ptr< CGAL::Polygon_with_holes_2<K,C> > > PolygonWithHolesVector ;

  std::cout << "Polygon_with_holes list with " << polies.size() << " element" << std::endl ;

  for( typename PolygonWithHolesVector::const_iterator pi = polies.begin() ; pi != polies.end() ; ++ pi )
    print_polygon_with_holes(**pi);
}


template<class K>
void print_straight_skeleton( CGAL::Straight_skeleton_2<K> const& ss )
{
  typedef CGAL::Straight_skeleton_2<K> Ss ;

  typedef typename Ss::Vertex_const_handle     Vertex_const_handle ;
  typedef typename Ss::Halfedge_const_handle   Halfedge_const_handle ;
  typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

  Halfedge_const_handle null_halfedge ;
  Vertex_const_handle   null_vertex ;

  std::cout << "Straight skeleton with " << ss.size_of_vertices()
            << " vertices, " << ss.size_of_halfedges()
            << " halfedges and " << ss.size_of_faces()
            << " faces" << std::endl ;

  std::cout << "Faces " << std::endl;
  for ( Halfedge_const_iterator h = ss.halfedges_begin(); h != ss.halfedges_end(); ++h )
  {
    if(h->is_inner_bisector() || h->is_bisector())
      continue;

    Halfedge_const_iterator end = h;
    for(;;)
    {
      if(h->is_inner_bisector())
        std::cout << "IBH" << h->id() << " " << std::flush ;
      else if(h->is_bisector())
        std::cout << "BH" << h->id() << " " << std::flush ;
      else
        std::cout << "CH" << h->id() << " " << std::flush ;

      print_vertex(h->vertex());

      h = h->next();
      if(h == end)
        break;

      std::cout << " ==> " << std::flush ;
    }

    std::cout << std::endl;
  }

  std::cout << "All halfedges: " << std::endl;

  for ( Halfedge_const_iterator h = ss.halfedges_begin(); h != ss.halfedges_end(); ++h )
  {
    if(h->is_inner_bisector())
      std::cout << "IBH" << h->id() << " " << std::flush ;
    else if(h->is_bisector())
      std::cout << "BH" << h->id() << " " << std::flush ;
    else
      std::cout << "CH" << h->id() << " " << std::flush ;

    print_vertex(h->prev()->vertex()) ;
    print_point(h->prev()->vertex()->point()) ;
    std::cout << " ==> " << std::flush ;
    print_vertex(h->vertex());
    print_point(h->vertex()->point()) ;
    std::cout << std::endl;
  }
}

} // namespace IO
} // namespace Straight_skeletons_2
} // namespace CGAL

#endif // CGAL_SLS_IO_PRINT_H
