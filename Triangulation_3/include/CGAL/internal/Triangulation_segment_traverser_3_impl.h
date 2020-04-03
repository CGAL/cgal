//A class that follows a straight line through a Delaunay triangulation structure.
//Copyright (C) 2012  Utrecht University
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s): Thijs van Lankveld

#ifndef CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H
#define CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H

#include <boost/array.hpp>

namespace CGAL {

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, Vertex_handle s, Vertex_handle t )
: _tr(tr) {
    CGAL_triangulation_precondition( !_tr->is_infinite(s) );
    CGAL_triangulation_precondition( !_tr->is_infinite(t) );
    CGAL_triangulation_precondition( s->point() != t->point() );
    CGAL_triangulation_precondition( _tr->dimension() >= 2 );

    _source = s->point();
    _target = t->point();
    _s_vertex = s;
    _t_vertex = t;
    
    Cell_handle c = s->cell();
    // If a vertex of an infinite cell, we start inside the convex hull.
    int inf;
    if( c->has_vertex( _tr->infinite_vertex(), inf ) )
        c = c->neighbor(inf);

    _cur = Simplex( c, Tr::VERTEX, c->index(s), -1 );

    jump_to_intersecting_cell();
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, Vertex_handle s, const Point& t )
: _tr(tr) {
    CGAL_triangulation_precondition( !_tr->is_infinite(s) );
    CGAL_triangulation_precondition( s->point() != t );
    CGAL_triangulation_precondition( _tr->dimension() >= 2 );
    CGAL_triangulation_precondition( _tr->dimension() == 3 ||
                                     orientation( *_tr->finite_facets_begin(), t ) == COPLANAR );

    _source = s->point();
    _target = t;
    _s_vertex = s;
    _t_vertex = Vertex_handle();

    Cell_handle c = s->cell();
    // If a vertex of an infinite cell, we start inside the convex hull.
    int inf;
    if( c->has_vertex( _tr->infinite_vertex(), inf ) )
        c = c->neighbor(inf);

    _cur = Simplex( c, Tr::VERTEX, c->index(s), -1 );

    jump_to_intersecting_cell();
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, const Point& s, Vertex_handle t, Cell_handle hint )
: _tr(tr) {
    CGAL_triangulation_precondition( !_tr->is_infinite(t) );
    CGAL_triangulation_precondition( s != t->point() );
    CGAL_triangulation_precondition( _tr->dimension() >= 2 );
    CGAL_triangulation_precondition( _tr->dimension() == 3 ||
                                     orientation( *_tr->finite_facets_begin(), s ) == COPLANAR );

    _source = s;
    _target = t->point();
    _s_vertex = Vertex_handle();
    _t_vertex = t;

    using CGAL::cpp11::get;
    get<0>(_cur) = _tr->locate( s, get<1>(_cur), get<2>(_cur), get<3>(_cur), hint );

    CGAL_triangulation_postcondition( get<0>(_cur) != Cell_handle() );

    jump_to_intersecting_cell();
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, const Point& s, const Point& t, Cell_handle hint )
: _tr(tr) {
    CGAL_triangulation_precondition( s != t );
    CGAL_triangulation_precondition( _tr->dimension() >= 2 );
    CGAL_triangulation_precondition( _tr->dimension() == 3 ||
                                     coplanar( *_tr->finite_facets_begin(), _target ) );

    _source = s;
    _target = t;
    _s_vertex = Vertex_handle();
    _t_vertex = Vertex_handle();

    using CGAL::cpp11::get;
    get<0>(_cur) = _tr->locate( s, get<1>(_cur), get<2>(_cur), get<3>(_cur), hint );

    CGAL_triangulation_postcondition( get<0>(_cur) != Cell_handle() );

    jump_to_intersecting_cell();
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, const Segment& s, Cell_handle hint )
: Triangulation_segment_cell_iterator_3<Tr,Inc>( tr, s.source(), s.target(), hint ) {}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr )
: _tr(tr) {}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>
Triangulation_segment_cell_iterator_3<Tr,Inc>::end() const {
    SCI sci(_tr);
    sci._source = _source;
    sci._target = _target;
    sci._s_vertex = _s_vertex;
    sci._t_vertex = _t_vertex;
    using CGAL::cpp11::get;
    get<0>(sci._cur) = Cell_handle();
    return sci;
}

template < class Tr, class Inc >
inline Triangulation_segment_cell_iterator_3<Tr,Inc>&
Triangulation_segment_cell_iterator_3<Tr,Inc>::operator++() {
    using CGAL::cpp11::get;
    CGAL_triangulation_precondition( get<0>(_cur) != Cell_handle() );
    increment();
    return *this;
}

template < class Tr, class Inc >
inline Triangulation_segment_cell_iterator_3<Tr,Inc>
Triangulation_segment_cell_iterator_3<Tr,Inc>::operator++( int ) {
    SCI tmp( *this );
    ++( *this );
    return tmp;
}

template < class Tr, class Inc >
inline typename Triangulation_segment_cell_iterator_3<Tr,Inc>::Cell_handle
Triangulation_segment_cell_iterator_3<Tr,Inc>::complete() {
    while( has_next() )
        increment();
    using CGAL::cpp11::get;
    return get<0>(_prev);
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator==( const SCI& sci ) const {
    // To be equal, the iterators must traverse the same triangulations
    // along the same line segment and they must have the same current cell.
    // Note that to limit cost, we just compare the triangulation pointers.
    using CGAL::cpp11::get;
    return ( &_tr == &sci._tr &&
             ( _s_vertex == Vertex_handle() ? _source == sci._source : _s_vertex == sci._s_vertex ) &&
             ( _t_vertex == Vertex_handle() ? _target == sci._target : _t_vertex == sci._t_vertex ) &&
             get<0>(_cur) == get<0>(sci._cur) );
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator!=( const SCI& sci ) const {
    return !( *this == sci );
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator==( Nullptr_t CGAL_triangulation_assertion_code(n) ) const {
    CGAL_triangulation_assertion( n == NULL );
    using CGAL::cpp11::get;
    return get<0>(_cur) == Cell_handle();
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator!=( Nullptr_t n ) const {
    return !( *this == n );
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr, Inc>::
jump_to_intersecting_cell()
{
  //copy current simplex
  Cell_handle ch = get<0>(_cur);
  Locate_type lt;
  int li, lj;
  entry(lt, li, lj);

  if (lt == Tr::FACET || lt == Tr::EDGE || lt == Tr::VERTEX)
  {
    //go forward once in the iterator
    Triangulation_segment_cell_iterator_3 cit(*this);
    ++cit;

    // collect the properties of "previous" to be in the right cell
    Cell_handle new_ch = cit.previous();
    if (new_ch == ch) //initialization already went to the good cell
      return;

    if (lt == Tr::VERTEX)
    {
      get<0>(_cur) = new_ch;
      //get<1>(_cur) is Locate_type and unchanged
      get<2>(_cur) = new_ch->index(ch->vertex(li));
      //get<3>(_cur) is lj and unchanged
      CGAL_assertion(get<0>(_cur)->vertex(get<2>(_cur)) == ch->vertex(li));
    }
    else if (lt == Tr::EDGE)
    {
      get<0>(_cur) = new_ch;
      //get<1>(_cur) is Locate_type and unchanged
      get<2>(_cur) = new_ch->index(ch->vertex(li));
      get<3>(_cur) = new_ch->index(ch->vertex(lj));
    }
    else
    {
      get<0>(_cur) = new_ch;
      //get<1>(_cur) is Locate_type and unchanged
      get<2>(_cur) = new_ch->index(ch);
      //get<3>(_cur) is lj and unchanged
    }
  }

  //for other cases (CELL, OUTSIDE_CONVEX_HULL, OUTSIDE_AFFINE_HULL)
  //the entry cell is already the good one, given by `locate()`
  //in the traverser constructor
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next() {
    CGAL_triangulation_precondition( has_next() );
    
    // Check if the target is in the current cell.
    int ti;
    using CGAL::cpp11::get;
    if( get<0>(_cur)->has_vertex( _t_vertex, ti ) ) {
        // The target is inside the cell.
        _prev = Simplex( get<0>(_cur), Tr::VERTEX, ti, -1 );
        get<0>(_cur) = Cell_handle();
        return;
    }

    // Walks to the next cell over a facet intersected by the line from source to target.
    // This method is based on Triangulation_3::locate().
    int inf;
    switch( _tr->dimension() ) {
        case 3: {
            // Infinite cells should be handled differently.
            if( get<0>(_cur)->has_vertex( _tr->infinite_vertex(), inf ) )
                walk_to_next_3_inf( inf );
            else
            {
              const Simplex backup = _cur;
              do {
                std::pair<Simplex, Simplex> p = walk_to_next_3(_prev, _cur);
                _prev = p.first;
                _cur = p.second;

              } while (get<0>(_cur) != Cell_handle()//end
                    && !get<0>(_cur)->has_vertex(_tr->infinite_vertex(), inf)
                    && have_same_entry(backup, _cur));
            }
            break;
        }
        case 2: {
            if( get<0>(_cur)->has_vertex( _tr->infinite_vertex(), inf ) )
                walk_to_next_2_inf( inf );
            else
                walk_to_next_2();
            break;
        }
    }
#ifdef CGAL_TRIANGULATION_3_TRAVERSER_CHECK_INTERSECTION
    if(_tr->dimension() == 3)
    {
      Cell_handle c = get<0>(_cur);
      if (c != Cell_handle() && !_tr->is_infinite(c)) //hard to say anything in this case
      {
        typename Tr::Segment seg(_source, _target);
        bool intersects = false;
        for (int i = 0; i < 4; ++i)
        {
          if (!_tr->is_infinite(c, i)
            && CGAL::do_intersect(_tr->triangle(c, i), seg))
          {
            intersects = true;
            break;
          }
        }
        CGAL_assertion(intersects);
      }
    }
#endif
}

template<class Tr, class Inc>
bool Triangulation_segment_cell_iterator_3<Tr, Inc>::
have_same_entry(const Simplex& s1, const Simplex& s2) const
{
  //type
  using CGAL::cpp11::get;
  if (get<1>(s1) != get<1>(s2))
    return false;
  switch (get<1>(s1))
  {
  case Locate_type::VERTEX:
    return get<0>(s1)->vertex(get<2>(s1)) == get<0>(s2)->vertex(get<2>(s2));
  case Locate_type::EDGE:
  {
    Vertex_handle v1a = get<0>(s1)->vertex(get<2>(s1));
    Vertex_handle v1b = get<0>(s1)->vertex(get<3>(s1));
    Vertex_handle v2a = get<0>(s2)->vertex(get<2>(s2));
    Vertex_handle v2b = get<0>(s2)->vertex(get<3>(s2));
    return (v1a == v2a && v1b == v2b)
        || (v1a == v2b && v1b == v2a);
  }
  case Locate_type::FACET:
    return triangulation()->are_equal(Facet(get<0>(s1), get<2>(s1)),
                                      Facet(get<0>(s2), get<2>(s2)));
  default:
    CGAL_assertion(false);
  };
  return false;
}

template < class Tr, class Inc >
std::pair<typename Triangulation_segment_cell_iterator_3<Tr, Inc>::Simplex,
          typename Triangulation_segment_cell_iterator_3<Tr, Inc>::Simplex >
Triangulation_segment_cell_iterator_3<Tr,Inc>::walk_to_next_3(const Simplex& prev,
                                                              const Simplex& cur) const
{
    using CGAL::cpp11::get;
    boost::array<const Point*, 4> vert
      = {&(get<0>(cur)->vertex(0)->point()),
         &(get<0>(cur)->vertex(1)->point()),
         &(get<0>(cur)->vertex(2)->point()),
         &(get<0>(cur)->vertex(3)->point()) };

    int inside=0,outside=0,regular_case=0,degenerate=0;
    Cell_handle nnext;

    if (get<1>(cur) == Tr::FACET) {
      regular_case = 1;
      int i = get<2>(cur);
      int j0 = Tr::vertex_triple_index(i, 0);
      int j1 = Tr::vertex_triple_index(i, 1);
      int j2 = Tr::vertex_triple_index(i, 2);
      Orientation o0 = _tr->orientation(_source, *vert[i], *vert[j0], _target);
      if (o0 == POSITIVE) {
        Orientation o1 = _tr->orientation(_source, *vert[i], *vert[j1], _target);
        if (o1 != POSITIVE) {
          if (_tr->orientation(*vert[i], *vert[j0], *vert[j1], _target) == POSITIVE) {
            nnext = get<0>(cur)->neighbor(j2);
            outside = j2;
            if (o1 == ZERO) degenerate = 1; //EDGE i j1
          }
          else
            inside = 1;
        }
        else {
          if (_tr->orientation(*vert[i], *vert[j1], *vert[j2], _target) == POSITIVE) {
            nnext = get<0>(cur)->neighbor(j0);
            outside = j0;
          }
          else
            inside = 2;
        }
      }
      else if (o0 == ZERO) {
        Orientation o1 = _tr->orientation(_source, *vert[i], *vert[j1], _target);
        if (o1 == NEGATIVE) {
          if (_tr->orientation(*vert[i], *vert[j0], *vert[j1], _target) == POSITIVE) {
            nnext = get<0>(cur)->neighbor(j2); //EDGE i j0
            degenerate = 2;
            outside = 44;
          }
          else
            inside = 3;
        }
        else if (o1 == ZERO) {
          if (_tr->orientation(*vert[i], *vert[j0], *vert[j2], _target) == POSITIVE)
            inside = 55;
          else
          {
            nnext = get<0>(cur)->neighbor(j2);  //VERTEX i
            degenerate = 3;
            outside = 5;
          }
        }
        else {
          if (_tr->orientation(*vert[i], *vert[j1], *vert[j2], _target) == POSITIVE) {
            nnext = get<0>(cur)->neighbor(j0);
            outside = j0;
          }
          else
            inside = 4;
        }
      }
      else {
        Orientation o2 = _tr->orientation(_source, *vert[i], *vert[j2], _target);
        if (o2 != NEGATIVE) {
          if (_tr->orientation(*vert[i], *vert[j2], *vert[j0], _target) == POSITIVE) {
            nnext = get<0>(cur)->neighbor(j1);
            outside = j1;
            if (o2 == ZERO) degenerate = 4; // EDGE i j2
          }
          else
            inside = 5;
        }
        else {
          if (_tr->orientation(*vert[i], *vert[j1], *vert[j2], _target) == POSITIVE) {
            nnext = get<0>(cur)->neighbor(j0);
            outside = j0;
          }
          else
            inside = 6;
        }
      }

      if ((!degenerate) && (!inside))
      {
        Simplex prev_after_walk(get<0>(cur), Tr::FACET, outside, -1);
        Simplex cur_after_walk( nnext,       Tr::FACET, nnext->index(get<0>(cur)), -1);
        return std::make_pair(prev_after_walk, cur_after_walk);
      }

      if ((!degenerate) && inside)
      {
        Simplex prev_after_walk(get<0>(cur),  Tr::CELL, -1, -1);
        Simplex cur_after_walk(Cell_handle(), Tr::OUTSIDE_AFFINE_HULL, -1, -1);
        return std::make_pair(prev_after_walk, cur_after_walk);
      }
    }


    // We check in which direction the target lies
    // by comparing its position relative to the planes through the
    // source and the edges of the cell.
    Orientation o[6];
    Orientation op[4];
    int pos = 0;
    // We keep track of which orientations are calculated.
    bool calc[6] = { false, false, false, false, false, false };

    if( get<1>(cur) == Tr::VERTEX ) {
        // The three planes through the vertex are set to coplanar.
        for( int j = 0; j < 4; ++j ) {
            if( get<2>(cur) != j ) {
                int ij = edgeIndex( get<2>(cur), j );
                o[ij] = COPLANAR;
                calc[ij] = true;
            }
        }
    }
    else if( get<1>(cur) == Tr::EDGE ) {
        // The plane through the edge is set to coplanar.
        int ij = edgeIndex( get<2>(cur), get<3>(cur) );
        o[ij] = COPLANAR;
        calc[ij] = true;
    }

    // For the remembering stochastic walk, we start trying with a random facet.
    int li = 0;
    CGAL_triangulation_assertion_code( bool incell = true; )
    for( int k = 0; k < 4; ++k, ++li )
    {
        // Skip the previous cell.
        Cell_handle next = get<0>(cur)->neighbor(li);
        if( next == get<0>(prev) )
        {
          op[li] = POSITIVE;
          pos += li;
          continue;
        }
        const Point* backup = vert[li];
        vert[li] = &_target;

        // Check if the target is on the opposite side of the supporting plane.
        op[li] = _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] );
        if( op[li] == POSITIVE )
            pos += li;
        if( op[li] != NEGATIVE ) {
            vert[li] = backup;
            continue;
        }
        CGAL_triangulation_assertion_code( incell = false; )

        // Check if the target is inside the 3-wedge with
        // the source as apex and the facet as an intersection.
          int lj = 0;
        int Or = 0;
        for( int l = 0; l < 4; ++l, ++lj ) {
            if( li == lj )
                continue;

            // We check the orientation of the target compared to the plane
            // Through the source and the edge opposite of ij.
            int oij = 5 - edgeIndex( li, lj );
            if( !calc[oij] ) {
                const Point* backup2 = vert[lj];
                vert[lj] = &_source;
                o[oij] = _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] );
                vert[lj] = backup2;
                calc[oij] = true;
            }

            if( o[oij] == POSITIVE ) {
                // The target is not inside the pyramid.
                // Invert the planes.
                // This can be safely done because either
                // they were not calculated yet,
                // or they will no longer be used.
                for( int j = 0; j < 4; ++j ) {
                    if( li == j ) continue;
                    int oij = 5 - edgeIndex( li, j );
                    o[oij] = -o[oij];
                }
                Or = 0;
                break;
            }
            else
                Or -= o[oij];
        }

        if( Or == 0 ) {
            // Either the target is not inside the pyramid,
            // or the pyramid is degenerate.
            vert[li] = backup;
            continue;
        }

        // The target is inside the pyramid.

        Simplex prev_after_walk;
        Simplex cur_after_walk;

        get<0>(prev_after_walk) = get<0>(cur);
        get<0>(cur_after_walk)  = next;
        switch( Or ) {
            case 3:
                get<1>(prev_after_walk) = Tr::FACET;
                get<2>(prev_after_walk) = li;
                get<1>(cur_after_walk) = Tr::FACET;
                get<2>(cur_after_walk) = get<0>(cur_after_walk)->index(get<0>(prev_after_walk));

                if(regular_case)
                {
                  CGAL_triangulation_assertion( get<0>(cur_after_walk)==nnext );
                  CGAL_triangulation_assertion( li==outside );
                  CGAL_triangulation_assertion( ! inside );
                }
                return std::make_pair(prev_after_walk, cur_after_walk);

            case 2:
                if(regular_case)
                  CGAL_triangulation_assertion(degenerate );

                get<1>(prev_after_walk) = Tr::EDGE;
                get<1>(cur_after_walk)  = Tr::EDGE;
                for( int j = 0; j < 4; ++j ) {
                    if( li != j && o[ 5 - edgeIndex(li, j) ] == COPLANAR) {
                        Edge opp = opposite_edge( get<0>(prev), li, j );
                        get<2>(prev_after_walk) = opp.second;
                        get<3>(prev_after_walk) = opp.third;
                        get<2>(cur_after_walk)
                          = get<0>(cur_after_walk)->index(
                              get<0>(prev_after_walk)->vertex( get<2>(prev_after_walk) ) );
                        get<3>(cur_after_walk)
                          = get<0>(cur_after_walk)->index(
                              get<0>(prev_after_walk)->vertex( get<3>(prev_after_walk) ) );

                        return std::make_pair(prev_after_walk, cur_after_walk);
                    }
                }
                CGAL_triangulation_assertion( false );
                return std::make_pair(prev, cur);
            case 1:
                if(regular_case)
                  CGAL_triangulation_assertion(degenerate );

                get<1>(prev_after_walk) = Tr::VERTEX;
                get<1>(cur_after_walk) = Tr::VERTEX;
                for( int j = 0; j < 4; ++j ) {
                    if( li != j && o[ 5 - edgeIndex(li, j) ] == NEGATIVE ) {
                        get<2>(prev_after_walk) = j;
                        get<2>(cur_after_walk)
                          = get<0>(cur_after_walk)->index(
                              get<0>(prev_after_walk)->vertex(j) );

                        return std::make_pair(prev_after_walk, cur_after_walk);
                    }
                }
                CGAL_triangulation_assertion( false );
                return std::make_pair(prev, cur);
            default:
                CGAL_triangulation_assertion( false );
                return std::make_pair(prev, cur);
        }
    }

    // The target lies inside this cell.
    Simplex prev_after_walk;
    CGAL_triangulation_assertion( incell );
    switch( op[0] + op[1] + op[2] + op[3] ) {
    case 4:
      CGAL_triangulation_assertion( pos == 6 );
      prev_after_walk = Simplex( get<0>(cur), Tr::CELL, -1, -1 );
      CGAL_triangulation_assertion( (! regular_case) || inside );
      break; 

    case 3:
      prev_after_walk = Simplex( get<0>(cur), Tr::FACET, 6-pos, -1 );
      break;
    case 2:
      if( pos < 3 )
        prev_after_walk = Simplex( get<0>(cur), Tr::EDGE, 0, pos+1 );
      else if( pos < 5 )
        prev_after_walk = Simplex( get<0>(cur), Tr::EDGE, 1, pos-1 );
      else
        prev_after_walk = Simplex( get<0>(cur), Tr::EDGE, 2, 3 );
      break;
    case 1:
      prev_after_walk = Simplex( get<0>(cur), Tr::VERTEX, pos, -1 );
      break;
    default:
      prev_after_walk = Simplex( get<0>(cur), Tr::OUTSIDE_AFFINE_HULL, -1, -1 );
      CGAL_triangulation_assertion( false );
    }

    Simplex cur_after_walk(Cell_handle(), Tr::OUTSIDE_AFFINE_HULL, -1, -1);
    return std::make_pair(prev_after_walk, cur_after_walk);
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_3_inf( int inf ) {
    using CGAL::cpp11::get;
    CGAL_triangulation_precondition( _tr->is_infinite( get<0>(_cur)->vertex(inf) ) );

    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = get<0>(_cur)->neighbor(inf);
    if( fin == get<0>(_prev) ) {
        _prev = _cur;
        get<0>(_cur) = Cell_handle();
        return;
    }

    boost::array < Point*, 4> vert;
    for( int i = 0; i != 4; ++i )
        if( i != inf )
            vert[i] = &(get<0>(_cur)->vertex(i)->point());
    vert[inf] = &_target;
    Orientation o[4];

    // Check if the target lies outside the convex hull.
    if( _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] ) == POSITIVE ) {
        // The target lies in an infinite cell.
        // Note that we do not traverse to other infinite cells.
        _prev = Simplex( get<0>(_cur), Tr::OUTSIDE_CONVEX_HULL, -1, -1 );
        get<0>(_cur) = Cell_handle();
        return;
    }

    vert[inf] = &(_source);
    CGAL_triangulation_assertion( _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] ) == POSITIVE );

    int li = 0;
    // Check if the line enters an adjacent infinite cell.
    // This occurs if the target lies on the other side of
    // a plane through one of the finite edges and the source point.
    for( int j = 0; j != 4; ++j, ++li ) {
        if( li == inf ) {
            o[li] = COPLANAR;
            continue;
        }
        
        // Skip the previous cell.
        Cell_handle next = get<0>(_cur)->neighbor(li);
        if( next == get<0>(_prev) ) {
            o[li] = POSITIVE;
            continue;
        }

        Point* backup = vert[li];
        vert[li] = &(_target);
        o[li] = _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] );

        if( o[li] != NEGATIVE ) {
            vert[li] = backup;
            continue;
        }

        // The target lies behind the plane through the source and two finite vertices.
        // Traverse to the incident infinite cell.
        CGAL_triangulation_assertion( _tr->is_infinite( next ) );
        _prev = Simplex( get<0>(_cur), Tr::FACET, li, -1 );
        _cur = Simplex( next, Tr::FACET, next->index( get<0>(_prev) ), -1 );
        return;
    }

    // The line enters the convex hull here (or lies on the finite facet).
    get<0>(_prev) = get<0>(_cur);
    get<0>(_cur) = fin;

    // Check through which simplex the line traverses.
    switch( o[0]+o[1]+o[2]+o[3] ) {
        case 3:
            get<1>(_prev) = Tr::FACET;
            get<2>(_prev) = inf;
            get<1>(_cur) = Tr::FACET;
            get<2>(_cur) = get<0>(_cur)->index(get<0>(_prev));
            return;
        case 2:
            get<1>(_prev) = Tr::EDGE;
            get<1>(_cur) = Tr::EDGE;
            for( int i = 0; i < 4; ++i ) {
                if( o[i] == COPLANAR && i != inf ) {
                    Edge opp = opposite_edge( get<0>(_prev), inf, i );
                    get<2>(_prev) = opp.second;
                    get<3>(_prev) = opp.third;
                    get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                    get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<3>(_prev) ) );
                    return;
                }
            }
            CGAL_triangulation_assertion( false );
            return;
        case 1:
            get<1>(_prev) = Tr::VERTEX;
            get<1>(_cur) = Tr::VERTEX;
            for( int i = 0; i < 4; ++i ) {
                if( o[i] == POSITIVE ) {
                    get<2>(_prev) = i;
                    get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(i) );
                    return;
                }
            }
            CGAL_triangulation_assertion( false );
            return;
        default:
            CGAL_triangulation_assertion( false );
            return;
    }
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_2()
{
    using CGAL::cpp11::get;
    boost::array<Point*, 3> vert
              = { &(get<0>(_cur)->vertex(0)->point()),
                  &(get<0>(_cur)->vertex(1)->point()),
                  &(get<0>(_cur)->vertex(2)->point()) };

    switch( get<1>(_cur) ) {
        case Tr::VERTEX: {
            // First we try the incident edges.
            Orientation ocw = CGAL::coplanar_orientation( *vert[get<2>(_cur)], *vert[_tr->cw(get<2>(_cur))], *vert[_tr->ccw(get<2>(_cur))], _target );
            if( get<0>(_cur)->neighbor( _tr->ccw(get<2>(_cur)) ) != get<0>(_prev) && ocw == NEGATIVE) {
                Cell_handle tmp = get<0>(_cur)->neighbor( _tr->ccw(get<2>(_cur)) );
                _prev = _cur;
                get<0>(_cur) = tmp;
                get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(get<2>(_cur)) );
                return;
            }
            Orientation occw = CGAL::coplanar_orientation( *vert[get<2>(_cur)], *vert[_tr->ccw(get<2>(_cur))], *vert[_tr->cw(get<2>(_cur))], _target );
            if( get<0>(_cur)->neighbor( _tr->cw(get<2>(_cur)) ) != get<0>(_prev) && occw == NEGATIVE) {
                Cell_handle tmp = get<0>(_cur)->neighbor( _tr->cw(get<2>(_cur)) );
                _prev = _cur;
                get<0>(_cur) = tmp;
                get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(get<2>(_cur)) );
                return;
            }

            // Then we try the opposite edge.
            Orientation op = CGAL::coplanar_orientation( *vert[_tr->ccw(get<2>(_cur))], *vert[_tr->cw(get<2>(_cur))], *vert[get<2>(_cur)], _target );
            if( op == NEGATIVE) {
                Cell_handle tmp = get<0>(_cur)->neighbor(get<2>(_cur));
                get<0>(_prev) = get<0>(_cur);
                get<0>(_cur) = tmp;

                switch( ocw+occw ) {
                    case 2:
                        get<1>(_prev) = Tr::EDGE;
                        get<2>(_prev) = _tr->ccw( get<2>(_cur) );
                        get<3>(_prev) = _tr->cw( get<2>(_cur) );
                        get<1>(_cur) = Tr::EDGE;
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<3>(_prev) ) );
                        return;
                    case 1:
                        get<1>(_prev) = Tr::VERTEX;
                        get<1>(_cur) = Tr::VERTEX;
                        if( ocw == COLLINEAR ) get<2>(_prev) = _tr->cw( get<2>(_cur) );
                        else get<2>(_cur) = _tr->ccw( get<2>(_cur) );
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        return;
                    default:
                        // The current vertex is the target.
                        CGAL_triangulation_assertion(false);
                        return;
                }
            }

            // The target lies in this cell.
            switch( ocw+occw+op ) {
            case 3:
                _prev = Simplex( get<0>(_cur), Tr::FACET, 3, -1 );
                break;
            case 2:
                if( ocw == 0 )
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr->ccw(get<2>(_cur)), -1 );
                else if( occw == 0 )
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr->cw(get<2>(_cur)), -1 );
                else
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, get<2>(_cur), -1 );
                break;
            case 1:
                if( ocw == 1 )
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->ccw(get<2>(_cur)), -1 );
                else if( occw == 1 )
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->cw(get<2>(_cur)), -1 );
                else
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, get<2>(_cur), -1 );
                break;
            case 0:
                CGAL_triangulation_assertion(false);
                _prev = Simplex( get<0>(_cur), Tr::OUTSIDE_AFFINE_HULL, -1, -1 );
                break;
            }
            get<0>(_cur) = Cell_handle();
            return;
        }
        case Tr::EDGE: {
            int lk = 3 - get<2>(_cur) - get<3>(_cur);

            if( get<0>(_cur)->neighbor(lk) != get<0>(_prev) ) {
                // Check the edge itself
                switch( CGAL::coplanar_orientation( *vert[get<2>(_cur)], *vert[get<3>(_cur)], *vert[lk], _target ) ) {
                    //_prev = _cur; //code not reached
                    case COLLINEAR:
                        // The target lies in this cell.
                        get<0>(_cur) = Cell_handle();
                        return;
                    case NEGATIVE: {
                        // The target lies opposite of the edge.
                        Cell_handle tmp = get<0>(_cur)->neighbor(lk);
                        get<0>(_cur) = tmp;
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(get<2>(_cur)) );
                        get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(get<3>(_cur)) );
                        return;
                    }
                    default:
                        break;
                }
            }

            Orientation o = CGAL::coplanar_orientation( _source, *vert[lk], *vert[get<2>(_cur)], _target );
            Orientation op;
            switch( o ) {
                case POSITIVE: {
                    // The ray passes through the edge ik.
                    op = CGAL::coplanar_orientation( *vert[lk], *vert[get<2>(_cur)], _source, _target );
                    if( op == NEGATIVE ) {
                        Cell_handle tmp = get<0>(_cur)->neighbor(get<3>(_cur));
                        get<0>(_prev) = get<0>(_cur);
                        get<0>(_cur) = tmp;

                        if( CGAL::collinear( _source, *vert[get<2>(_cur)], _target ) ) {
                            get<1>(_prev) = Tr::VERTEX;
                            get<2>(_prev) = get<2>(_cur);
                            get<1>(_cur) = Tr::VERTEX;
                            get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        }
                        else {
                            get<1>(_prev) = Tr::EDGE;
                            get<2>(_prev) = get<2>(_cur);
                            get<3>(_prev) = lk;
                            get<1>(_cur) = Tr::EDGE;
                            get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                            get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<3>(_prev) ) );
                        }
                        return;
                    }
                    break;
                }
                default: {
                    // The ray passes through the edge jk.
                    op = CGAL::coplanar_orientation( *vert[lk], *vert[get<3>(_cur)], _source, _target );
                    if( op == NEGATIVE ) {
                        Cell_handle tmp = get<0>(_cur)->neighbor(get<2>(_cur));
                        get<0>(_prev) = get<0>(_cur);
                        get<0>(_cur) = tmp;

                        if( CGAL::collinear( _source, *vert[get<3>(_cur)], _target ) ) {
                            get<1>(_prev) = Tr::VERTEX;
                            get<2>(_prev) = get<3>(_cur);
                            get<1>(_cur) = Tr::VERTEX;
                            get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        }
                        else if( o == COLLINEAR ) {
                            get<1>(_prev) = Tr::VERTEX;
                            get<2>(_prev) = lk;
                            get<1>(_cur) = Tr::VERTEX;
                            get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        }
                        else {
                            get<1>(_prev) = Tr::EDGE;
                            get<2>(_prev) = lk;
                            get<3>(_prev) = get<3>(_cur);
                            get<1>(_cur) = Tr::EDGE;
                            get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                            get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        }
                        return;
                    }
                    break;
                }
            }

            // The target lies in this cell.
            if( op == POSITIVE )
                _prev = Simplex( get<0>(_cur), Tr::FACET, 3, -1 );
            else {
                CGAL_triangulation_assertion( op == ZERO );
                switch( o ) {
                case POSITIVE:
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, get<2>(_cur), lk );
                    break;
                case NEGATIVE:
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, get<3>(_cur), lk );
                    break;
                case ZERO:
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, lk, -1 );
                    break;
                }
            }
            get<0>(_cur) = Cell_handle();
            return;
        }
        case Tr::FACET: {
          int li = 0;
          
            Orientation o[3];
            bool calc[3] = { false, false, false };

            for( int j = 0; j != 3; ++j, li = _tr->ccw(li) ) {
                Cell_handle next = get<0>(_cur)->neighbor(li);
                if( next == get<0>(_prev) )
                    continue;

                // The target should lie on the other side of the edge.
                Orientation op = CGAL::coplanar_orientation( *vert[_tr->ccw(li)], *vert[_tr->cw(li)], *vert[li], _target );
                if( op == POSITIVE )
                    continue;

                // The target should lie inside the wedge.
                if( !calc[_tr->ccw(li)] ) {
                    o[_tr->ccw(li)] = CGAL::coplanar_orientation( _source, *vert[_tr->ccw(li)], *vert[_tr->cw(li)], _target );
                    calc[_tr->ccw(li)] = true;
                }
                if( o[_tr->ccw(li)] == NEGATIVE )
                    continue;
                else if( op == COLLINEAR && o[_tr->ccw(li)] == COLLINEAR ) {
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->ccw(li), -1 );
                    get<0>(_cur) = Cell_handle();
                    return;
                }

                if( !calc[_tr->cw(li)] ) {
                    o[_tr->cw(li)] = CGAL::coplanar_orientation( _source, *vert[_tr->cw(li)], *vert[li], _target );
                    calc[_tr->cw(li)] = true;
                }
                if( o[_tr->cw(li)] == POSITIVE )
                    continue;
                else if( op == COLLINEAR && o[_tr->cw(li)] == COLLINEAR ) {
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->cw(li), -1 );
                    get<0>(_cur) = Cell_handle();
                    return;
                }

                get<0>(_prev) = get<0>(_cur);
                get<0>(_cur) = next;

                switch( o[_tr->ccw(li)] + o[_tr->cw(li)] ) {
                    case 2:
                        get<1>(_prev) = Tr::EDGE;
                        get<2>(_prev) = _tr->ccw(li);
                        get<3>(_prev) = _tr->cw(li);
                        get<1>(_cur) = Tr::EDGE;
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( _tr->ccw(li) ) );
                        get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( _tr->cw(li) ) );
                        return;
                    case 1:
                        get<1>(_prev) = Tr::VERTEX;
                        get<1>(_cur) = Tr::VERTEX;
                        if( o[_tr->ccw(li)] == COLLINEAR ) get<2>(_prev) = _tr->ccw(li);
                        else get<2>(_prev) = _tr->cw(li);
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        return;
                    default:
                        CGAL_triangulation_assertion( false );
                        return;
                }
            }

            // The target lies in this cell.
            _prev = Simplex( get<0>(_cur), Tr::FACET, 3, -1 );
            get<0>(_cur) = Cell_handle();
            return;
        }
        default:
        CGAL_triangulation_assertion( false );
    }
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_2_inf( int inf ) {
    using CGAL::cpp11::get;
    CGAL_triangulation_precondition( _tr->is_infinite( get<0>(_cur)->vertex(3) ) );
    CGAL_triangulation_precondition( _tr->is_infinite( get<0>(_cur)->vertex(inf) ) );
	
    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = get<0>(_cur)->neighbor(inf);
    if (fin == get<0>(_prev)) {
        _prev = _cur;
        get<0>(_cur) = Cell_handle();
        return;
    }

    // Check the neighboring cells.
    Orientation occw = CGAL::coplanar_orientation( _source,
      get<0>(_cur)->vertex( _tr->ccw(inf))->point(),
      get<0>(_cur)->vertex(_tr->cw(inf))->point(),
      _target );
    if( occw == NEGATIVE ) {
        Cell_handle tmp = get<0>(_cur)->neighbor(_tr->cw(inf));
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr->ccw(inf), inf );
        _cur = Simplex( tmp, Tr::EDGE, tmp->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), tmp->index( get<0>(_prev)->vertex( get<3>(_prev) ) ) );
        return;
    }
    Orientation ocw = CGAL::coplanar_orientation( _source,
      get<0>(_cur)->vertex( _tr->cw(inf))->point(),
      get<0>(_cur)->vertex(_tr->ccw(inf))->point(),
      _target );
    if( ocw == NEGATIVE ) {
        Cell_handle tmp = get<0>(_cur)->neighbor(_tr->ccw(inf));
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr->cw(inf), inf );
        _cur = Simplex( tmp, Tr::EDGE, tmp->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), tmp->index( get<0>(_prev)->vertex( get<3>(_prev) ) ) );
        return;
    }
    Orientation op = CGAL::coplanar_orientation(
      get<0>(_cur)->vertex( _tr->ccw(inf) )->point(),
      get<0>(_cur)->vertex( _tr->cw(inf) )->point(),
      _source, _target );
    switch( op ) {
    case NEGATIVE:
        if( occw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->ccw(inf), -1 );
            _cur = Simplex( fin, Tr::VERTEX, fin->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), -1 );
            return;
        }
        if( ocw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->cw(inf), -1 );
            _cur = Simplex( fin, Tr::VERTEX, fin->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), -1 );
            return;
        }
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr->ccw(inf), _tr->cw(inf) );
        _cur = Simplex( fin, Tr::EDGE, fin->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), fin->index( get<0>(_prev)->vertex( get<3>(_prev) ) ) );
        return;
    case COLLINEAR:
        if( occw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->ccw(inf), -1 );
            get<0>(_cur) = Cell_handle();
            return;
        }
        if( ocw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr->cw(inf), -1 );
            get<0>(_cur) = Cell_handle();
            return;
        }
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr->ccw(inf), _tr->cw(inf) );
        get<0>(_cur) = Cell_handle();
        return;
    case POSITIVE:
        // The target lies in this infinite cell.
        _prev = Simplex( get<0>(_cur), Tr::OUTSIDE_CONVEX_HULL, -1, -1 );
        get<0>(_cur) = Cell_handle();
        return;
    }
}

template < class Tr, class Inc >
CGAL::Orientation
Triangulation_segment_cell_iterator_3<Tr, Inc>::orientation(
  const Facet& f, const Point& p) const
{
  return _tr->orientation(
    f.first->vertex(Tr::vertex_triple_index(f.second, 0))->point(),
    f.first->vertex(Tr::vertex_triple_index(f.second, 1))->point(),
    f.first->vertex(Tr::vertex_triple_index(f.second, 2))->point(),
    p);
}

template < class Tr, class Inc >
bool
Triangulation_segment_cell_iterator_3<Tr, Inc>::coplanar(
  const Facet& f, const Point& p) const
{
  return orientation(f, p) == COPLANAR;
}

template < class Tr, class Inc >
typename Triangulation_segment_cell_iterator_3<Tr, Inc>::Edge
Triangulation_segment_cell_iterator_3<Tr, Inc>::opposite_edge(
  Cell_handle c, int li, int lj) const
{
  CGAL_triangulation_precondition(li >= 0 && li < 4);
  CGAL_triangulation_precondition(lj >= 0 && lj < 4);
  CGAL_triangulation_precondition(li != lj);

  switch (6 - li - lj) { // i + j + missing indices = 6.
    case 1: return Edge(c, 0, 1);
    case 2: return Edge(c, 0, 2);
    case 3: return (li == 0 || lj == 0) ? Edge(c, 1, 2) : Edge(c, 0, 3);
    case 4: return Edge(c, 1, 3);
    case 5: return Edge(c, 2, 3);
  }

  CGAL_triangulation_assertion(false);
  return Edge();
}

template < class Tr, class Inc >
typename Triangulation_segment_cell_iterator_3<Tr, Inc>::Edge
Triangulation_segment_cell_iterator_3<Tr, Inc>::opposite_edge(
  const Edge& e) const
{
  return opposite_edge(e.first, e.second, e.third);
}


} //end of CGAL namespace

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H
