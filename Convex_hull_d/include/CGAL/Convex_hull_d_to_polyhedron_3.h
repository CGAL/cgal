// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel

#ifndef CGAL_CONVEX_HULL_D_TO_POLYHEDRON_3_H
#define CGAL_CONVEX_HULL_D_TO_POLYHEDRON_3_H

#include <CGAL/license/Convex_hull_d.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Convex_hull_d_to_polyhedron_3.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS \
  "The Triangulation package (see http://doc.cgal.org/latest/Triangulation) should be used instead."
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Convex_hull_d.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

namespace CGAL {

template <class _HDS, class _ChullType>
class Build_polyhedron_from_chull : public Modifier_base< _HDS>
{
  public:
    typedef _HDS                           HDS;
    typedef _ChullType                     ChullType;
    typedef typename HDS::Vertex           HdsVertex;
    typedef typename HdsVertex::Point      Point;
    typedef typename ChullType::Facet_const_handle  Facet_handle;
    typedef typename ChullType::Vertex_const_handle Vertex_handle;
    typedef typename ChullType::R R;
    typedef typename ChullType::Point_d Point_d;

Build_polyhedron_from_chull(const ChullType& CH) : ch(CH) {}

Point convert(const Point_d& p) const
{ return Point(p.hx(),p.hy(),p.hz(),p.hw()); }

void operator()(HDS& hds )
{
  CGAL_assertion( ch.current_dimension() == 3); 
  Polyhedron_incremental_builder_3<HDS> B(hds, true);
  B.begin_surface( ch.number_of_vertices(), 
                   ch.number_of_facets() ); 
  // would be nice to have statistical data on
  // Chull available other than print_statistics()
  Unique_hash_map<Vertex_handle, int> index( -1);
  std::list<Facet_handle> L = ch.all_facets();
  typename std::list<Facet_handle>::iterator fit;
  Facet_handle f;
  Vertex_handle v;
  int i = 0;
  for ( fit = L.begin(); fit != L.end(); ++fit) { f = *fit;
    for (int k=0; k < 3; ++k) {
       v = ch.vertex_of_facet(f,k);
       if ( index[v] == -1 ) {
         B.add_vertex( convert(ch.associated_point(v)) );
         index[v] = i++;
       }
    }
  }
  Point_d center = ch.center();
  for ( fit = L.begin(); fit != L.end(); ++fit) { f = *fit;
     B.begin_facet();
     Vertex_handle v0 = ch.vertex_of_facet(f,0);
     Vertex_handle v1 = ch.vertex_of_facet(f,1);
     Vertex_handle v2 = ch.vertex_of_facet(f,2);
     typename R::Orientation_d orientation =
       ch.kernel().orientation_d_object();
     std::vector<Point_d> V(4);
     V[0]= ch.associated_point(v0);
     V[1]= ch.associated_point(v1);
     V[2]= ch.associated_point(v2);
     V[3]= center;
     if ( orientation( V.begin(), V.end() ) == POSITIVE ) {
       B.add_vertex_to_facet( index[v0] );
       B.add_vertex_to_facet( index[v1] );
       B.add_vertex_to_facet( index[v2] );
     } else {
       B.add_vertex_to_facet( index[v0] );
       B.add_vertex_to_facet( index[v2] );
       B.add_vertex_to_facet( index[v1] );
     }
     B.end_facet();
  }
  B.end_surface();
}

private:
  const ChullType& ch;
};

/*{\Mtext \headerline{Low Dimensional Conversion Routine}
include |<CGAL/Convex_hull_d_to_polyhedron_3.h>|
\setopdims{2cm}{3cm}}*/

template <class R, class Polyhedron_3>
void convex_hull_d_to_polyhedron_3(
  const Convex_hull_d<R>& C, Polyhedron_3& P)
/*{\Mfunc converts the convex hull |C| to polyedral surface stored in 
   |P|.\\ \precond |dim == 3| and |dcur == 3|. }*/
{ typedef Convex_hull_d<R> ChullType;

  typedef typename Polyhedron_3::HalfedgeDS  HDS;
  CGAL_assertion_msg(C.dimension()==3&&C.current_dimension()==3,
  "convex_hull_d_to_polyhedron_3: only full manifold can be transformed.");
  Build_polyhedron_from_chull<HDS,ChullType> get_surface(C);
  P.delegate( get_surface );
}

} //namespace CGAL
#endif //CGAL_CONVEX_HULL_D_TO_POLYHEDRON_3_H
