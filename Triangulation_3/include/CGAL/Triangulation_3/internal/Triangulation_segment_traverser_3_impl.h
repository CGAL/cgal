// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Thijs van Lankveld, Jane Tournois

#ifndef CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H
#define CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/assertions.h>

namespace CGAL {

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, Vertex_handle s, Vertex_handle t )
: _tr(tr) {
    CGAL_precondition( !_tr->is_infinite(s) );
    CGAL_precondition( !_tr->is_infinite(t) );
    CGAL_precondition( s->point() != t->point() );
    CGAL_precondition( _tr->dimension() >= 2 );

    _source = s->point();
    _target = t->point();
    _s_vertex = s;
    _t_vertex = t;

    Cell_handle c = s->cell();
    // If a vertex of an infinite cell, we start inside the convex hull.
    int inf;
    if( c->has_vertex( _tr->infinite_vertex(), inf ) )
        c = c->neighbor(inf);

    _cur = Simplex{ c, Tr::VERTEX, c->index(s), -1 };

    jump_to_intersecting_cell();
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, Vertex_handle s, const Point& t )
: _tr(tr) {
    CGAL_precondition( !_tr->is_infinite(s) );
    CGAL_precondition( s->point() != t );
    CGAL_precondition( _tr->dimension() >= 2 );
    CGAL_precondition( _tr->dimension() == 3 ||
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

    _cur = Simplex{ c, Tr::VERTEX, c->index(s), -1 };

    jump_to_intersecting_cell();
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, const Point& s, Vertex_handle t, Cell_handle hint )
: _tr(tr) {
    CGAL_precondition( !_tr->is_infinite(t) );
    CGAL_precondition( s != t->point() );
    CGAL_precondition( _tr->dimension() >= 2 );
    CGAL_precondition( _tr->dimension() == 3 ||
                                     orientation( *_tr->finite_facets_begin(), s ) == COPLANAR );

    _source = s;
    _target = t->point();
    _s_vertex = Vertex_handle();
    _t_vertex = t;

    cell() = _tr->locate( s, lt(), li(), lj(), hint );

    CGAL_postcondition( cell() != Cell_handle() );

    jump_to_intersecting_cell();
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr* tr, const Point& s, const Point& t, Cell_handle hint )
: _tr(tr) {
    CGAL_precondition( s != t );
    CGAL_precondition( _tr->dimension() >= 2 );
    CGAL_precondition( _tr->dimension() == 3 ||
                                     coplanar( *_tr->finite_facets_begin(), _target ) );

    _source = s;
    _target = t;
    _s_vertex = Vertex_handle();
    _t_vertex = Vertex_handle();

    cell() = _tr->locate( s, lt(), li(), lj(), hint );

    CGAL_postcondition( cell() != Cell_handle() );

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
    sci._cur.cell = Cell_handle();
    return sci;
}

template < class Tr, class Inc >
inline Triangulation_segment_cell_iterator_3<Tr,Inc>&
Triangulation_segment_cell_iterator_3<Tr,Inc>::operator++() {
    CGAL_precondition( cell() != Cell_handle() );
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
    return prev_cell();
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator==( const SCI& sci ) const {
    // To be equal, the iterators must traverse the same triangulations
    // and they must have the same current cell.
    return ( _tr == sci._tr &&
             cell() == sci.cell() );
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator!=( const SCI& sci ) const {
    return !( *this == sci );
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator==( Nullptr_t CGAL_assertion_code(n) ) const {
    CGAL_assertion( n == NULL );
    return cell() == Cell_handle();
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
  Cell_handle ch = cell();
  Locate_type lt;
  int li, clj;
  entry(lt, li, clj);

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
      cell() = new_ch;
      //lt() is Locate_type and unchanged
      this->li() = new_ch->index(ch->vertex(li));
      //lj() is lj and unchanged
      CGAL_assertion(cell()->vertex(this->li()) == ch->vertex(li));
    }
    else if (lt == Tr::EDGE)
    {
      cell() = new_ch;
      //lt() is Locate_type and unchanged
      this->li() = new_ch->index(ch->vertex(li));
      this->lj() = new_ch->index(ch->vertex(clj));
    }
    else
    {
      cell() = new_ch;
      //lt() is Locate_type and unchanged
      this->li() = new_ch->index(ch);
      //lj() is lj and unchanged
    }
  }

  //for other cases (CELL, OUTSIDE_CONVEX_HULL, OUTSIDE_AFFINE_HULL)
  //the entry cell is already the good one, given by `locate()`
  //in the traverser constructor
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next() {
    CGAL_precondition( has_next() );

    // Check if the target is in the current cell.
    int ti;
    if( cell()->has_vertex( _t_vertex, ti ) ) {
        // The target is inside the cell.
        _prev = Simplex{ cell(), Tr::VERTEX, ti, -1 };
        cell() = Cell_handle();
        lt() = Locate_type::VERTEX;
        return;
    }

    // Walks to the next cell over a facet intersected by the line from source to target.
    // This method is based on Triangulation_3::locate().
    int inf;
    switch( _tr->dimension() ) {
        case 3: {
            // Infinite cells should be handled differently.
            if( cell()->has_vertex( _tr->infinite_vertex(), inf ) )
                walk_to_next_3_inf( inf );
            else
            {
              const Simplex backup = _cur;
              do {
                std::pair<Simplex, Simplex> p = walk_to_next_3(_prev, _cur);
                _prev = p.first;
                _cur = p.second;

              } while (cell() != Cell_handle()//end
                    && !cell()->has_vertex(_tr->infinite_vertex(), inf)
                    && have_same_entry(backup, _cur));
            }
            break;
        }
        case 2: {
            if( cell()->has_vertex( _tr->infinite_vertex(), inf ) )
                walk_to_next_2_inf( inf );
            else
                walk_to_next_2();
            break;
        }
    }
#ifdef CGAL_TRIANGULATION_3_TRAVERSER_CHECK_INTERSECTION
    if(_tr->dimension() == 3)
    {
      Cell_handle c = cell();
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
  if (s1.lt != s2.lt)
    return false;
  switch (s1.lt)
  {
  case Locate_type::VERTEX:
    return s1.cell->vertex(s1.li) == s2.cell->vertex(s2.li);
  case Locate_type::EDGE:
  {
    Vertex_handle v1a = s1.cell->vertex(s1.li);
    Vertex_handle v1b = s1.cell->vertex(s1.lj);
    Vertex_handle v2a = s2.cell->vertex(s2.li);
    Vertex_handle v2b = s2.cell->vertex(s2.lj);
    return (v1a == v2a && v1b == v2b)
        || (v1a == v2b && v1b == v2a);
  }
  case Locate_type::FACET:
    return triangulation()->are_equal(Facet(s1.cell, s1.li),
                                      Facet(s2.cell, s2.li));
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
  const auto cur_cell = cur.cell;
  std::array<const Point*, 4> vert = {&(cur_cell->vertex(0)->point()), &(cur_cell->vertex(1)->point()),
                                      &(cur_cell->vertex(2)->point()), &(cur_cell->vertex(3)->point())};

  CGAL_assertion_code(int inside = 0;)
  CGAL_assertion_code(int outside = 0;)
  int regular_case = 0, degenerate = 0;

  if(cur.lt == Tr::FACET && prev.cell != Cell_handle()) {
    // [source, target] entered the cell `cur` via a facet.
    // Note that, if prev.cell == Cell_handle(), that means `source` is *on*
    // the facet, and the block of this `if` cannot be applied.
    Simplex prev_after_walk;
    Simplex cur_after_walk;

    auto case_target_is_inside_cur_cell = [&](int case_nb) {
      CGAL_USE(case_nb);
      CGAL_assertion_code(inside = case_nb;)
      prev_after_walk = {cur_cell, Tr::CELL, -1, -1};
      cur_after_walk = {{}, Tr::CELL, -1, -1};
    };
    auto case_segment_exits_cur_cell_by = [&](int facet_nb, Cell_handle nnext = {}) {
      if(nnext == Cell_handle{}) {
        nnext = cur_cell->neighbor(facet_nb);
      }
      CGAL_assertion_code(outside = facet_nb;)
      prev_after_walk = {cur_cell, Tr::FACET, facet_nb, -1};
      cur_after_walk = {nnext, Tr::FACET, nnext->index(cur_cell), -1};
    };
    regular_case = 1;
    const int i = cur.li;
    const int j0 = Tr::vertex_triple_index(i, 0);
    const int j1 = Tr::vertex_triple_index(i, 1);
    const int j2 = Tr::vertex_triple_index(i, 2);
    Orientation o0 = _tr->orientation(_source, *vert[i], *vert[j0], _target);
    if(o0 == POSITIVE) { // o0 > 0
      Orientation o1 = _tr->orientation(_source, *vert[i], *vert[j1], _target);
      if(o1 != POSITIVE) { // o1 <= 0
        Orientation oi01 = _tr->orientation(*vert[i], *vert[j0], *vert[j1], _target);
        if(oi01 == POSITIVE) {
          case_segment_exits_cur_cell_by(j2);
          if(o1 == ZERO)
            degenerate = 1; // EDGE i j1
        } else {            // o0 > 0, o1 <= 0, oi01 <= 0
          case_target_is_inside_cur_cell(1);
          if(oi01 == ZERO) { // on FACET j2 (i, j0, j1)
            degenerate = 1;
          } // end oi01 == ZERO
        }
      } // end  o1 <= 0
      else
      { // o1 > 0
        Orientation oi12 = _tr->orientation(*vert[i], *vert[j1], *vert[j2], _target);
        if(oi12 == POSITIVE) {
          case_segment_exits_cur_cell_by(j0);
        } else { // o0 > 0, o1 > 0, oi12 <= 0
          case_target_is_inside_cur_cell(2);
          if(oi12 == ZERO) { // on FACET j0 (i, j1, j2)
            degenerate = 1;
          } // end oi12 == ZERO
        }
      }
    } // end o0 > 0
    else if(o0 == ZERO)
    {
      // target is on plane (source, vert[i], vert[j0])
      Orientation o1 = _tr->orientation(_source, *vert[i], *vert[j1], _target);
      if(o1 == NEGATIVE) {
        Orientation oi12 = _tr->orientation(*vert[i], *vert[j0], *vert[j1], _target);
        if(oi12 == POSITIVE) {
          degenerate = 2;
          case_segment_exits_cur_cell_by(44, cur_cell->neighbor(j2)); // EDGE i j0
        } else {
          case_target_is_inside_cur_cell(3);
          if(oi12 == ZERO) { // target is *on* EDGE i j0
            degenerate = 1;
          }
        }
      } else if(o1 == ZERO) {
        // o0 == o1 == 0 -> target is on line source-vert[i]
        if(_tr->orientation(*vert[i], *vert[j0], *vert[j2], _target) == POSITIVE)
          case_target_is_inside_cur_cell(55);
        else {
          degenerate = 3;
          case_segment_exits_cur_cell_by(5, cur_cell->neighbor(j2)); // VERTEX i
        }
      } else { // o0 == 0, o1 > 0
        Orientation oi12 = _tr->orientation(*vert[i], *vert[j1], *vert[j2], _target);
        if(oi12 == POSITIVE) {
          case_segment_exits_cur_cell_by(j0);
        } else {
          case_target_is_inside_cur_cell(4);
          if(oi12 == ZERO) { // on FACET j0 (i, j1, j2)
            degenerate = 1;
          } // end oi12 == ZERO
        }
      }
    } // end o0 == 0
    else
    { // o0 < 0
      Orientation o2 = _tr->orientation(_source, *vert[i], *vert[j2], _target);
      if(o2 != NEGATIVE) {
        // o2 >= 0
        Orientation oi20 = _tr->orientation(*vert[i], *vert[j2], *vert[j0], _target);
        if(oi20 == POSITIVE) {
          case_segment_exits_cur_cell_by(j1);
          if(o2 == ZERO)
            degenerate = 4; // EDGE i j2
        } else {
          case_target_is_inside_cur_cell(5);
          if(oi20 == ZERO) { // on FACET j1 (i, j2, j0)
            degenerate = 1;
          }
        }
      } else {
        Orientation oi12 = _tr->orientation(*vert[i], *vert[j1], *vert[j2], _target);
        if(oi12 == POSITIVE) {
          case_segment_exits_cur_cell_by(j0);
        } else {
          case_target_is_inside_cur_cell(6);
          if(oi12 == ZERO) { // on FACET j0 (i, j1, j2)
            degenerate = 1;
          }
        }
      }
    }

    if(!degenerate) {
        return {prev_after_walk, cur_after_walk};
    }
  }

  // We check in which direction the target lies
  // by comparing its position relative to the planes through the
  // source and the edges of the cell.
  std::array<Orientation, 6> o;
  std::array<Orientation, 4> op;
  int pos = 0;
  // We keep track of which orientations are calculated.
  bool calc[6] = {false, false, false, false, false, false};

  if(cur.lt == Tr::VERTEX) {
    // The three planes through the vertex are set to coplanar.
    for(int j = 0; j < 4; ++j) {
        if(cur.li != j) {
          int ij = edgeIndex(cur.li, j);
          o[ij] = COPLANAR;
          calc[ij] = true;
        }
    }
  } else if(cur.lt == Tr::EDGE) {
    // The plane through the edge is set to coplanar.
    int ij = edgeIndex(cur.li, cur.lj);
    o[ij] = COPLANAR;
    calc[ij] = true;
  }

  // For the remembering stochastic walk, we start trying with a random facet.
  CGAL_assertion_code(bool incell = true;)

  for(int li = 0; li < 4; ++li)
  {
    // Skip the previous cell.
    Cell_handle next = cur_cell->neighbor(li);
    if(next == prev.cell) {
        op[li] = POSITIVE;
        pos += li;
        continue;
    }
    const Point* const backup_vert_li = std::exchange(vert[li], &_target);
    bool op_li_is_null = false;
    if(_t_vertex != Vertex_handle()) {
      for(int i = 0; i < 4; ++i) {
        if(li == i) continue;
        if(cur_cell->vertex(i) == _t_vertex) op_li_is_null = true;
      }
    }
    // Check if the target is on the opposite side of the supporting plane.
    op[li] = op_li_is_null ? ZERO : _tr->orientation(*vert[0], *vert[1], *vert[2], *vert[3]);
    if(op[li] == POSITIVE)
        pos += li;
    if(op[li] != NEGATIVE) {
        vert[li] = backup_vert_li;
        continue;
    }
    CGAL_assertion_code(incell = false;)

    // Check if the target is inside the 3-wedge with
    // the source as apex and the facet as an intersection.
    int Or = 0;
    for(int lj = 0; lj < 4; ++lj) {
      if(li == lj)
        continue;
      // We check the orientation of the target compared to the plane
      // Through the source and the edge opposite of ij.
      const int oij = 5 - edgeIndex(li, lj);
      if(!calc[oij]) {
        const Point* const backup_vert_lj = std::exchange(vert[lj], &_source);
        o[oij] = _tr->orientation(*vert[0], *vert[1], *vert[2], *vert[3]);
        vert[lj] = backup_vert_lj;
        calc[oij] = true;
      }
      if(o[oij] == POSITIVE) {
        // The target is not inside the pyramid.
        // Invert the planes.
        for(int j = 0; j < 4; ++j) {
          if(li == j)
            continue;
          int oij = 5 - edgeIndex(li, j);
          if(calc[oij])
            o[oij] = -o[oij];
        }
        Or = 0;
        break;
      } else
        Or -= o[oij];
    }

    if(Or == 0) {
      // Either the target is not inside the pyramid,
      // or the pyramid is degenerate.
      vert[li] = backup_vert_li;
      continue;
    }

    // The target is inside the pyramid.
    switch(Or) {
    case 3: {
      if(regular_case) {
        CGAL_assertion(li == outside);
        CGAL_assertion(!inside);
      }
      return {{cur_cell, Tr::FACET, li}, {next, Tr::FACET, next->index(cur_cell)}};
    }
    case 2: {
      if(regular_case)
        CGAL_assertion(degenerate);
      for(int j = 0; j < 4; ++j) {
        if(li != j && o[5 - edgeIndex(li, j)] == COPLANAR) {
          Edge opp = opposite_edge(prev.cell, li, j);
          return {
              {cur_cell, Tr::EDGE, opp.second, opp.third},
              {next, Tr::EDGE, next->index(cur_cell->vertex(opp.second)), next->index(cur_cell->vertex(opp.third))}};
        }
      }
      CGAL_unreachable();
      return std::make_pair(prev, cur);
    }
    case 1:
      if(regular_case)
        CGAL_assertion(degenerate);
      for(int j = 0; j < 4; ++j) {
        if(li != j && o[5 - edgeIndex(li, j)] == NEGATIVE) {
          return {{cur_cell, Tr::VERTEX, j}, {next, Tr::VERTEX, next->index(cur_cell->vertex(j))}};
        }
      }
      CGAL_unreachable();
      return std::make_pair(prev, cur);
    default:
        CGAL_unreachable();
        return std::make_pair(prev, cur);
    }
    CGAL_unreachable();
  }

  // The target lies inside this cell.
  CGAL_assertion( incell );
  return {
    std::invoke([&]() -> Simplex {
      switch( op[0] + op[1] + op[2] + op[3] ) {
      case 4:
        CGAL_assertion( pos == 6 );
        CGAL_assertion( (! regular_case) || inside );
        return { cur_cell, Tr::CELL };
        break;
      case 3:
        return { cur_cell, Tr::FACET, 6 - pos };
        break;
      case 2:
        if( pos < 3 ) // first is 0
          return { cur_cell, Tr::EDGE, 0, pos };
        else if( pos < 5 ) { // could be (0, pos), or (1, pos-1)
          if(op[0] == POSITIVE)
            return { cur_cell, Tr::EDGE, 0, pos };
          else
            return { cur_cell, Tr::EDGE, 1, pos-1 };
        }
        else
          return { cur_cell, Tr::EDGE, 2, 3 };
        break;
      case 1:
        return { cur_cell, Tr::VERTEX, pos };
        break;
      default:
        CGAL_unreachable();
      }
    }),
    { Cell_handle() }
  };
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_3_inf( int inf )
{
    CGAL_precondition( _tr->is_infinite( cell()->vertex(inf) ) );

    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = cell()->neighbor(inf);
    if( fin == prev_cell() ) {
        _prev = _cur;
        prev_lt() = Tr::CELL;
        cell() = Cell_handle();
        lt() = Tr::CELL;
        return;
    }

    std::array < Point*, 4> vert;
    for( int i = 0; i != 4; ++i )
        if( i != inf )
            vert[i] = &(cell()->vertex(i)->point());
    vert[inf] = &_target;
    Orientation o[4];

    // Check if the target lies outside the convex hull.
    if( _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] ) == POSITIVE ) {
        // The target lies in an infinite cell.
        // Note that we do not traverse to other infinite cells.
        _prev = Simplex{ cell(), Tr::OUTSIDE_CONVEX_HULL, -1, -1 };
        cell() = Cell_handle();
        return;
    }

    vert[inf] = &(_source);
    CGAL_assertion( _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] ) == POSITIVE );

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
        Cell_handle next = cell()->neighbor(li);
        if( next == prev_cell() ) {
            o[li] = POSITIVE;
            continue;
        }

        Point* backup_vert_li = vert[li];
        vert[li] = &(_target);
        o[li] = _tr->orientation( *vert[0], *vert[1], *vert[2], *vert[3] );

        if( o[li] != NEGATIVE ) {
            vert[li] = backup_vert_li;
            continue;
        }

        // The target lies behind the plane through the source and two finite vertices.
        // Traverse to the incident infinite cell.
        CGAL_assertion( _tr->is_infinite( next ) );
        _prev = Simplex{ cell(), Tr::FACET, li, -1 };
        _cur = Simplex{ next, Tr::FACET, next->index(prev_cell()), -1 };
        return;
    }

    // The line enters the convex hull here (or lies on the finite facet).
    prev_cell() = cell();
    cell() = fin;

    // Check through which simplex the line traverses.
    switch( o[0]+o[1]+o[2]+o[3] ) {
        case 3:
            prev_lt() = Tr::FACET;
            prev_li() = inf;
            lt() = Tr::FACET;
            this->li() = cell()->index(prev_cell());
            return;
        case 2:
            prev_lt() = Tr::EDGE;
            lt() = Tr::EDGE;
            for( int i = 0; i < 4; ++i ) {
                if( o[i] == COPLANAR && i != inf ) {
                    Edge opp = opposite_edge( prev_cell(), inf, i );
                    prev_li() = opp.second;
                    prev_lj() = opp.third;
                    this->li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                    this->lj() = cell()->index( prev_cell()->vertex( prev_lj() ) );
                    return;
                }
            }
            CGAL_unreachable();
            return;
        case 1:
            prev_lt() = Tr::VERTEX;
            lt() = Tr::VERTEX;
            for( int i = 0; i < 4; ++i ) {
                if( o[i] == POSITIVE ) {
                    prev_li() = i;
                    this->li() = cell()->index( prev_cell()->vertex(i) );
                    return;
                }
            }
            CGAL_unreachable();
            return;
        default:
            CGAL_unreachable();
            return;
    }
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_2()
{
    const std::array<const Point*, 3> vert
              = { &(cell()->vertex(0)->point()),
                  &(cell()->vertex(1)->point()),
                  &(cell()->vertex(2)->point()) };

    typename Gt::Coplanar_orientation_3 coplanar_orientation
      = _tr->geom_traits().coplanar_orientation_3_object();

    switch( lt() ) {
        case Tr::VERTEX: {
            // First we try the incident edges.
            Orientation ocw = coplanar_orientation( *vert[li()], *vert[_tr->cw(li())], *vert[_tr->ccw(li())], _target );
            if( cell()->neighbor( _tr->ccw(li()) ) != prev_cell() && ocw == NEGATIVE) {
                Cell_handle tmp = cell()->neighbor( _tr->ccw(li()) );
                _prev = _cur;
                cell() = tmp;
                li() = cell()->index( prev_cell()->vertex(li()) );
                return;
            }
            Orientation occw = coplanar_orientation( *vert[li()], *vert[_tr->ccw(li())], *vert[_tr->cw(li())], _target );
            if( cell()->neighbor( _tr->cw(li()) ) != prev_cell() && occw == NEGATIVE) {
                Cell_handle tmp = cell()->neighbor( _tr->cw(li()) );
                _prev = _cur;
                cell() = tmp;
                li() = cell()->index( prev_cell()->vertex(li()) );
                return;
            }

            // Then we try the opposite edge.
            Orientation op = coplanar_orientation( *vert[_tr->ccw(li())], *vert[_tr->cw(li())], *vert[li()], _target );
            if( op == NEGATIVE) {
                Cell_handle tmp = cell()->neighbor(li());
                prev_cell() = cell();
                cell() = tmp;

                switch( ocw+occw ) {
                    case 2:
                        prev_lt() = Tr::EDGE;
                        prev_li() = _tr->ccw( li() );
                        prev_lj() = _tr->cw( li() );
                        lt() = Tr::EDGE;
                        li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                        lj() = cell()->index( prev_cell()->vertex( prev_lj() ) );
                        return;
                    case 1:
                        prev_lt() = Tr::VERTEX;
                        lt() = Tr::VERTEX;
                        if( ocw == COLLINEAR ) prev_li() = _tr->cw( li() );
                        else li() = _tr->ccw( li() );
                        li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                        return;
                    default:
                        // The current vertex is the target.
                        CGAL_unreachable();
                        return;
                }
            }

            // The target lies in this cell.
            switch( ocw+occw+op ) {
            case 3:
              _prev = Simplex{ cell(), Tr::FACET, 3, -1 };
              break;
            case 2:
                if( ocw == 0 )
                  _prev = Simplex{ cell(), Tr::EDGE, _tr->ccw(li()), -1 };
                else if( occw == 0 )
                  _prev = Simplex{ cell(), Tr::EDGE, _tr->cw(li()), -1 };
                else
                  _prev = Simplex{ cell(), Tr::EDGE, li(), -1 };
                break;
            case 1:
                if( ocw == 1 )
                  _prev = Simplex{ cell(), Tr::VERTEX, _tr->ccw(li()), -1 };
                else if( occw == 1 )
                  _prev = Simplex{ cell(), Tr::VERTEX, _tr->cw(li()), -1 };
                else
                  _prev = Simplex{ cell(), Tr::VERTEX, li(), -1 };
                break;
            case 0:
                CGAL_unreachable();
                break;
            }
            cell() = Cell_handle();
            return;
        }
        case Tr::EDGE: {
            int lk = 3 - li() - lj();

            if( cell()->neighbor(lk) != prev_cell() ) {
                // Check the edge itself
                switch(coplanar_orientation( *vert[li()], *vert[lj()], *vert[lk], _target ) ) {
                    //_prev = _cur; //code not reached
                    case COLLINEAR:
                        // The target lies in this cell.
                        cell() = Cell_handle();
                        return;
                    case NEGATIVE: {
                        // The target lies opposite of the edge.
                        Cell_handle tmp = cell()->neighbor(lk);
                        cell() = tmp;
                        li() = cell()->index( prev_cell()->vertex(li()) );
                        lj() = cell()->index( prev_cell()->vertex(lj()) );
                        return;
                    }
                    default:
                        break;
                }
            }

            typename Gt::Collinear_3 collinear
              = _tr->geom_traits().collinear_3_object();
            Orientation o = coplanar_orientation( _source, *vert[lk], *vert[li()], _target );
            Orientation op;
            switch( o ) {
                case POSITIVE: {
                    // The ray passes through the edge ik.
                    op = coplanar_orientation( *vert[lk], *vert[li()], _source, _target );
                    if( op == NEGATIVE ) {
                        Cell_handle tmp = cell()->neighbor(lj());
                        prev_cell() = cell();
                        cell() = tmp;

                        if( collinear( _source, *vert[li()], _target ) ) {
                            prev_lt() = Tr::VERTEX;
                            prev_li() = li();
                            lt() = Tr::VERTEX;
                            li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                        }
                        else {
                            prev_lt() = Tr::EDGE;
                            prev_li() = li();
                            prev_lj() = lk;
                            lt() = Tr::EDGE;
                            li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                            lj() = cell()->index( prev_cell()->vertex( prev_lj() ) );
                        }
                        return;
                    }
                    break;
                }
                default: {
                    // The ray passes through the edge jk.
                    op = coplanar_orientation( *vert[lk], *vert[lj()], _source, _target );
                    if( op == NEGATIVE ) {
                        Cell_handle tmp = cell()->neighbor(li());
                        prev_cell() = cell();
                        cell() = tmp;

                        if( collinear( _source, *vert[lj()], _target ) ) {
                            prev_lt() = Tr::VERTEX;
                            prev_li() = lj();
                            lt() = Tr::VERTEX;
                            li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                        }
                        else if( o == COLLINEAR ) {
                            prev_lt() = Tr::VERTEX;
                            prev_li() = lk;
                            lt() = Tr::VERTEX;
                            li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                        }
                        else {
                            prev_lt() = Tr::EDGE;
                            prev_li() = lk;
                            prev_lj() = lj();
                            lt() = Tr::EDGE;
                            li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                            lj() = cell()->index( prev_cell()->vertex( prev_li() ) );
                        }
                        return;
                    }
                    break;
                }
            }

            // The target lies in this cell.
            if( op == POSITIVE )
              _prev = Simplex{ cell(), Tr::FACET, 3, -1 };
            else {
                CGAL_assertion( op == ZERO );
                switch( o ) {
                case POSITIVE:
                  _prev = Simplex{ cell(), Tr::EDGE, li(), lk };
                  break;
                case NEGATIVE:
                  _prev = Simplex{ cell(), Tr::EDGE, lj(), lk };
                  break;
                case ZERO:
                  _prev = Simplex{ cell(), Tr::VERTEX, lk, -1 };
                  break;
                }
            }
            cell() = Cell_handle();
            return;
        }
        case Tr::FACET: {
          int li = 0;

            Orientation o[3];
            bool calc[3] = { false, false, false };

            for( int j = 0; j != 3; ++j, li = _tr->ccw(li) ) {
                Cell_handle next = cell()->neighbor(li);
                if( next == prev_cell() )
                    continue;

                // The target should lie on the other side of the edge.
                Orientation op = coplanar_orientation( *vert[_tr->ccw(li)], *vert[_tr->cw(li)], *vert[li], _target );
                if( op == POSITIVE )
                    continue;

                // The target should lie inside the wedge.
                if( !calc[_tr->ccw(li)] ) {
                    o[_tr->ccw(li)] = coplanar_orientation( _source, *vert[_tr->ccw(li)], *vert[_tr->cw(li)], _target );
                    calc[_tr->ccw(li)] = true;
                }
                if( o[_tr->ccw(li)] == NEGATIVE )
                    continue;
                else if( op == COLLINEAR && o[_tr->ccw(li)] == COLLINEAR ) {
                    _prev = Simplex{ cell(), Tr::VERTEX, _tr->ccw(li), -1 };
                    cell() = Cell_handle();
                    return;
                }

                if( !calc[_tr->cw(li)] ) {
                    o[_tr->cw(li)] = coplanar_orientation( _source, *vert[_tr->cw(li)], *vert[li], _target );
                    calc[_tr->cw(li)] = true;
                }
                if( o[_tr->cw(li)] == POSITIVE )
                    continue;
                else if( op == COLLINEAR && o[_tr->cw(li)] == COLLINEAR ) {
                    _prev = Simplex{ cell(), Tr::VERTEX, _tr->cw(li), -1 };
                    cell() = Cell_handle();
                    return;
                }

                prev_cell() = cell();
                cell() = next;

                switch( o[_tr->ccw(li)] + o[_tr->cw(li)] ) {
                    case 2:
                        prev_lt() = Tr::EDGE;
                        prev_li() = _tr->ccw(li);
                        prev_lj() = _tr->cw(li);
                        lt() = Tr::EDGE;
                        this->li() = cell()->index( prev_cell()->vertex( _tr->ccw(li) ) );
                        this->lj() = cell()->index( prev_cell()->vertex( _tr->cw(li) ) );
                        return;
                    case 1:
                        prev_lt() = Tr::VERTEX;
                        this->lt() = Tr::VERTEX;
                        if( o[_tr->ccw(li)] == COLLINEAR ) prev_li() = _tr->ccw(li);
                        else prev_li() = _tr->cw(li);
                        this->li() = cell()->index( prev_cell()->vertex( prev_li() ) );
                        return;
                    default:
                        CGAL_unreachable();
                        return;
                }
            }

            // The target lies in this cell.
            _prev = Simplex{ cell(), Tr::FACET, 3, -1 };
            cell() = Cell_handle();
            return;
        }
        default:
        CGAL_unreachable();
    }
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_2_inf( int inf )
{
    CGAL_precondition( _tr->is_infinite( cell()->vertex(3) ) );
    CGAL_precondition( _tr->is_infinite( cell()->vertex(inf) ) );

    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = cell()->neighbor(inf);
    if (fin == prev_cell()) {
        _prev = _cur;
        cell() = Cell_handle();
        return;
    }

    typename Gt::Coplanar_orientation_3 coplanar_orientation
      = _tr->geom_traits().coplanar_orientation_3_object();

    // Check the neighboring cells.
    Orientation occw = coplanar_orientation( _source,
      cell()->vertex( _tr->ccw(inf))->point(),
      cell()->vertex(_tr->cw(inf))->point(),
      _target );
    if( occw == NEGATIVE ) {
        Cell_handle tmp = cell()->neighbor(_tr->cw(inf));
        _prev = Simplex{ cell(), Tr::EDGE, _tr->ccw(inf), inf };
        _cur = Simplex{ tmp, Tr::EDGE, tmp->index(prev_cell()->vertex(prev_li())), tmp->index(prev_cell()->vertex(prev_lj())) };
        return;
    }
    Orientation ocw = coplanar_orientation( _source,
      cell()->vertex( _tr->cw(inf))->point(),
      cell()->vertex(_tr->ccw(inf))->point(),
      _target );
    if( ocw == NEGATIVE ) {
        Cell_handle tmp = cell()->neighbor(_tr->ccw(inf));
        _prev = Simplex{ cell(), Tr::EDGE, _tr->cw(inf), inf };
        _cur = Simplex{ tmp, Tr::EDGE, tmp->index(prev_cell()->vertex(prev_li())), tmp->index(prev_cell()->vertex(prev_lj())) };
        return;
    }
    Orientation op = coplanar_orientation(
      cell()->vertex( _tr->ccw(inf) )->point(),
      cell()->vertex( _tr->cw(inf) )->point(),
      _source, _target );
    switch( op ) {
    case NEGATIVE:
        if( occw == COLLINEAR ) {
          _prev = Simplex{ cell(), Tr::VERTEX, _tr->ccw(inf), -1 };
          _cur = Simplex{ fin, Tr::VERTEX, fin->index(prev_cell()->vertex(prev_li())), -1 };
          return;
        }
        if( ocw == COLLINEAR ) {
          _prev = Simplex{ cell(), Tr::VERTEX, _tr->cw(inf), -1 };
          _cur = Simplex{ fin, Tr::VERTEX, fin->index(prev_cell()->vertex(prev_li())), -1 };
          return;
        }
        _prev = Simplex{ cell(), Tr::EDGE, _tr->ccw(inf), _tr->cw(inf) };
        _cur = Simplex{ fin, Tr::EDGE, fin->index(prev_cell()->vertex(prev_li())), fin->index(prev_cell()->vertex(prev_lj())) };
        return;
    case COLLINEAR:
        if( occw == COLLINEAR ) {
          _prev = Simplex{ cell(), Tr::VERTEX, _tr->ccw(inf), -1 };
          cell() = Cell_handle();
          return;
        }
        if( ocw == COLLINEAR ) {
          _prev = Simplex{ cell(), Tr::VERTEX, _tr->cw(inf), -1 };
          cell() = Cell_handle();
          return;
        }
        _prev = Simplex{ cell(), Tr::EDGE, _tr->ccw(inf), _tr->cw(inf) };
        cell() = Cell_handle();
        return;
    case POSITIVE:
        // The tarstd::std::get lies in this infinite cell.
      _prev = Simplex{ cell(), Tr::OUTSIDE_CONVEX_HULL, -1, -1 };
      cell() = Cell_handle();
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
  CGAL_precondition(li >= 0 && li < 4);
  CGAL_precondition(lj >= 0 && lj < 4);
  CGAL_precondition(li != lj);

  switch (6 - li - lj) { // i + j + missing indices = 6.
    case 1: return Edge(c, 0, 1);
    case 2: return Edge(c, 0, 2);
    case 3: return (li == 0 || lj == 0) ? Edge(c, 1, 2) : Edge(c, 0, 3);
    case 4: return Edge(c, 1, 3);
    case 5: return Edge(c, 2, 3);
  }

  CGAL_unreachable();
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
