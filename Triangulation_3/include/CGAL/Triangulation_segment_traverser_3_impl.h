//A class that follows a straight line through a Delaunay triangulation structure.
//Copyright (C) 2012  Utrecht University
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s): Thijs van Lankveld

#ifndef CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H
#define CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H

namespace CGAL {

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr& tr, Vertex_handle s, Vertex_handle t )
: _tr(tr) {
    CGAL_triangulation_precondition( !_tr.is_infinite(s) );
    CGAL_triangulation_precondition( !_tr.is_infinite(t) );
    CGAL_triangulation_precondition( s->point() != t->point() );
    CGAL_triangulation_precondition( _tr.dimension() >= 2 );

    _source = s->point();
    _target = t->point();
    _s_vertex = s;
    _t_vertex = t;
    _s_vert = s;
    _t_vert = t;
    
    Cell_handle c = s->cell();
    // If a vertex of an infinite cell, we start inside the convex hull.
    int inf;
    if( c->has_vertex( _tr.infinite_vertex(), inf ) )
        c = c->neighbor(inf);

    _cur = Simplex( c, Tr::VERTEX, c->index(s), -1 );
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr& tr, Vertex_handle s, const Point& t )
: _tr(tr) {
    CGAL_triangulation_precondition( !_tr.is_infinite(s) );
    CGAL_triangulation_precondition( s->point() != t );
    CGAL_triangulation_precondition( _tr.dimension() >= 2 );
    CGAL_triangulation_precondition( _tr.dimension() == 3 ||
                                     _tr.orientation( *_tr.finite_facets_begin(), t ) == COPLANAR );

    _source = s->point();
    _target = t;
    _s_vertex = s;
    _t_vertex = Vertex_handle();
    _s_vert = s;
    _t_vert = _tds2.create_vertex( Vertex(_target) );

    Cell_handle c = s->cell();
    // If a vertex of an infinite cell, we start inside the convex hull.
    int inf;
    if( c->has_vertex( _tr.infinite_vertex(), inf ) )
        c = c->neighbor(inf);

    _cur = Simplex( c, Tr::VERTEX, c->index(s), -1 );
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr& tr, const Point& s, Vertex_handle t, Cell_handle hint )
: _tr(tr) {
    CGAL_triangulation_precondition( !_tr.is_infinite(t) );
    CGAL_triangulation_precondition( s != t->point() );
    CGAL_triangulation_precondition( _tr.dimension() >= 2 );
    CGAL_triangulation_precondition( _tr.dimension() == 3 ||
                                     _tr.orientation( *_tr.finite_facets_begin(), s ) == COPLANAR );

    _source = s;
    _target = t->point();
    _s_vertex = Vertex_handle();
    _t_vertex = t;
    _s_vert = _tds2.create_vertex( Vertex(_source) );
    _t_vert = t;
    
    get<0>(_cur) = _tr.locate( s, get<1>(_cur), get<2>(_cur), get<3>(_cur), hint );

    CGAL_triangulation_postcondition( get<0>(_cur) != Cell_handle() );
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr& tr, const Point& s, const Point& t, Cell_handle hint )
: _tr(tr) {
    CGAL_triangulation_precondition( s != t );
    CGAL_triangulation_precondition( _tr.dimension() >= 2 );
    CGAL_triangulation_precondition( _tr.dimension() == 3 ||
                                     _tr.coplanar( *_tr.finite_facets_begin(), _target ) );

    _source = s;
    _target = t;
    _s_vertex = Vertex_handle();
    _t_vertex = Vertex_handle();
    _s_vert = _tds2.create_vertex( Tr::Vertex(_source) );
    _t_vert = _tds2.create_vertex( Tr::Vertex(_target) );

    get<0>(_cur) = _tr.locate( s, get<1>(_cur), get<2>(_cur), get<3>(_cur), hint );

    CGAL_triangulation_postcondition( get<0>(_cur) != Cell_handle() );
}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr& tr, const Segment& s, Cell_handle hint )
: Triangulation_segment_cell_iterator_3<Tr,Inc>( tr, s.source(), s.target(), hint ) {}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>::
Triangulation_segment_cell_iterator_3( const Tr& tr )
: _tr(tr) {}

template < class Tr, class Inc >
Triangulation_segment_cell_iterator_3<Tr,Inc>
Triangulation_segment_cell_iterator_3<Tr,Inc>::end() const {
    SCI sci(_tr);
    sci._source = _source;
    sci._target = _target;
    sci._s_vertex = _s_vertex;
    sci._t_vertex = _t_vertex;
    get<0>(sci._cur) = Cell_handle();
    return sci;
}

template < class Tr, class Inc >
inline Triangulation_segment_cell_iterator_3<Tr,Inc>&
Triangulation_segment_cell_iterator_3<Tr,Inc>::operator++() {
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
    return get<0>(_prev);
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator==( const SCI& sci ) const {
    // To be equal, the iterators must traverse the same triangulations
    // along the same line segment and they must have the same current cell.
    // Note that to limit cost, we just compare the triangulation pointers.
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
    return get<0>(_cur) == Cell_handle();
}

template < class Tr, class Inc >
inline bool Triangulation_segment_cell_iterator_3<Tr,Inc>::
operator!=( Nullptr_t n ) const {
    return !( *this == n );
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next() {
    CGAL_triangulation_precondition( has_next() );
    
    // Check if the target is in the current cell.
    int ti;
    if( get<0>(_cur)->has_vertex( _t_vertex, ti ) ) {
        // The target is inside the cell.
        _prev = Simplex( get<0>(_cur), Tr::VERTEX, ti, -1 );
        get<0>(_cur) = Cell_handle();
        return;
    }

    // Walks to the next cell over a facet intersected by the line from source to target.
    // This method is based on Triangulation_3::locate().
    int inf;
    switch( _tr.dimension() ) {
        case 3: {
            // Infinite cells should be handled differently.
            if( get<0>(_cur)->has_vertex( _tr.infinite_vertex(), inf ) )
                walk_to_next_3_inf( inf );
            else
                walk_to_next_3();
            break;
        }
        case 2: {
            if( get<0>(_cur)->has_vertex( _tr.infinite_vertex(), inf ) )
                walk_to_next_2_inf( inf );
            else
                walk_to_next_2();
            break;
        }
    }
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_3() {
    Vertex_handle vert[4] = { get<0>(_cur)->vertex(0),
                              get<0>(_cur)->vertex(1),
                              get<0>(_cur)->vertex(2),
                              get<0>(_cur)->vertex(3) };

    // We check in which direction the target lies
    // by comparing its position relative to the planes through the
    // source and the edges of the cell.
    Orientation o[6];
    Orientation op[4]; int pos = 0;
    // We keep track of which orientations are calculated.
    bool calc[6] = { false, false, false, false, false, false };

    if( get<1>(_cur) == Tr::VERTEX ) {
        // The three planes through the vertex are set to coplanar.
        for( int j = 0; j < 4; ++j ) {
            if( get<2>(_cur) != j ) {
                int ij = edgeIndex( get<2>(_cur), j );
                o[ij] = COPLANAR;
                calc[ij] = true;
            }
        }
    }
    else if( get<1>(_cur) == Tr::EDGE ) {
        // The plane through the edge is set to coplanar.
        int ij = edgeIndex( get<2>(_cur), get<3>(_cur) );
        o[ij] = COPLANAR;
        calc[ij] = true;
    }

    // For the remembering stochastic walk, we start trying with a random facet.
    int li = rng.template get_bits<2>();
    CGAL_triangulation_assertion_code( bool incell = true; )
    for( int k = 0; k < 4; ++k, li = _tr.increment_index(li) ) {
        // Skip the previous cell.
        Cell_handle next = get<0>(_cur)->neighbor(li);
        if( next == get<0>(_prev) )
            continue;
  
        Vertex_handle backup = vert[li];
        vert[li] = _t_vert;
       
        // Check if the target is on the opposite side of the supporting plane.
        op[li] = _tr.orientation( vert[0], vert[1], vert[2], vert[3] );
        if( op[li] == POSITIVE )
            pos += li;
        if( op[li] != NEGATIVE ) {
            vert[li] = backup;
            continue;
        }
        CGAL_triangulation_assertion_code( incell = false; )

        // Check if the target is inside the 3-wedge with
        // the source as apex and the facet as an intersection.
        int lj = rng.template get_bits<2>();
        int Or = 0;
        for( int l = 0; l < 4; ++l, lj = _tr.increment_index(lj) ) {
            if( li == lj )
                continue;

            // We check the orientation of the target compared to the plane
            // Through the source and the edge opposite of ij.
            int oij = 5 - edgeIndex( li, lj );
            if( !calc[oij] ) {
                Vertex_handle backup2 = vert[lj];
                vert[lj] = _s_vert;
                o[oij] = _tr.orientation( vert[0], vert[1], vert[2], vert[3] );
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
        get<0>(_prev) = get<0>(_cur);
        get<0>(_cur) = next;
        switch( Or ) {
            case 3:
                get<1>(_prev) = Tr::FACET;
                get<2>(_prev) = li;
                get<1>(_cur) = Tr::FACET;
                get<2>(_cur) = get<0>(_cur)->index(get<0>(_prev));
                return;
            case 2:
                get<1>(_prev) = Tr::EDGE;
                get<1>(_cur) = Tr::EDGE;
                for( int j = 0; j < 4; ++j ) {
                    if( li != j && o[ 5 - edgeIndex(li, j) ] == COPLANAR) {
                        Edge opp = _tr.opposite_edge( get<0>(_prev), li, j );
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
                for( int j = 0; j < 4; ++j ) {
                    if( li != j && o[ 5 - edgeIndex(li, j) ] == NEGATIVE ) {
                        get<2>(_prev) = j;
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(j) );
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

    // The target lies inside this cell.
    CGAL_triangulation_assertion( incell );
    switch( op[0] + op[1] + op[2] + op[3] ) {
    case 4:
        CGAL_triangulation_assertion( pos == 6 );
        _prev = Simplex( get<0>(_cur), Tr::CELL, -1, -1 );
        break;
    case 3:
        _prev = Simplex( get<0>(_cur), Tr::FACET, 6-pos, -1 );
        break;
    case 2:
        if( pos < 3 )
            _prev = Simplex( get<0>(_cur), Tr::EDGE, 0, pos+1 );
        else if( pos < 5 )
            _prev = Simplex( get<0>(_cur), Tr::EDGE, 1, pos-1 );
        else
            _prev = Simplex( get<0>(_cur), Tr::EDGE, 2, 3 );
        break;
    case 1:
        _prev = Simplex( get<0>(_cur), Tr::VERTEX, pos, -1 );
        break;
    default:
        _prev = Simplex( get<0>(_cur), Tr::OUTSIDE_AFFINE_HULL, -1, -1 );
        CGAL_triangulation_assertion( false );
    }
    get<0>(_cur) = Cell_handle();
    return;
}

template < class Tr, class Inc >
void Triangulation_segment_cell_iterator_3<Tr,Inc>::
walk_to_next_3_inf( int inf ) {
    CGAL_triangulation_precondition( _tr.is_infinite( get<0>(_cur)->vertex(inf) ) );

    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = get<0>(_cur)->neighbor(inf);
    if( fin == get<0>(_prev) ) {
        _prev = _cur;
        get<0>(_cur) = Cell_handle();
        return;
    }

    Vertex_handle vert[4];
    for( int i = 0; i != 4; ++i )
        if( i != inf )
            vert[i] = get<0>(_cur)->vertex(i);
    vert[inf] = _t_vert;
    Orientation o[4];

    // Check if the target lies outside the convex hull.
    if( _tr.orientation( vert[0], vert[1], vert[2], vert[3] ) == POSITIVE ) {
        // The target lies in an infinite cell.
        // Note that we do not traverse to other infinite cells.
        _prev = Simplex( get<0>(_cur), Tr::OUTSIDE_CONVEX_HULL, -1, -1 );
        get<0>(_cur) = Cell_handle();
        return;
    }

    vert[inf] = _s_vert;
    CGAL_triangulation_assertion( _tr.orientation( vert[0], vert[1], vert[2], vert[3] ) == POSITIVE );

    // For the remembering stochastic walk, we start trying with a random index:
    int li = rng.template get_bits<2>();

    // Check if the line enters an adjacent infinite cell.
    // This occurs if the target lies on the other side of
    // a plane through one of the finite edges and the source point.
    for( int j = 0; j != 4; ++j, li = _tr.increment_index(li) ) {
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

        Vertex_handle backup = vert[li];
        vert[li] = _t_vert;
        o[li] = _tr.orientation( vert[0], vert[1], vert[2], vert[3] );

        if( o[li] != NEGATIVE ) {
            vert[li] = backup;
            continue;
        }

        // The target lies behind the plane through the source and two finite vertices.
        // Traverse to the incident infinite cell.
        CGAL_triangulation_assertion( _tr.is_infinite( next ) );
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
                    Edge opp = _tr.opposite_edge( get<0>(_prev), inf, i );
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
walk_to_next_2() {
    Vertex_handle vert[3] = { get<0>(_cur)->vertex(0),
                              get<0>(_cur)->vertex(1),
                              get<0>(_cur)->vertex(2) };

    switch( get<1>(_cur) ) {
        case Tr::VERTEX: {
            // First we try the incident edges.
            Orientation ocw = _tr.coplanar_orientation( vert[get<2>(_cur)], vert[_tr.cw(get<2>(_cur))], vert[_tr.ccw(get<2>(_cur))], _t_vert );
            if( get<0>(_cur)->neighbor( _tr.ccw(get<2>(_cur)) ) != get<0>(_prev) && ocw == NEGATIVE) {
                Cell_handle tmp = get<0>(_cur)->neighbor( _tr.ccw(get<2>(_cur)) );
                _prev = _cur;
                get<0>(_cur) = tmp;
                get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(get<2>(_cur)) );
                return;
            }
            Orientation occw = _tr.coplanar_orientation( vert[get<2>(_cur)], vert[_tr.ccw(get<2>(_cur))], vert[_tr.cw(get<2>(_cur))], _t_vert );
            if( get<0>(_cur)->neighbor( _tr.cw(get<2>(_cur)) ) != get<0>(_prev) && occw == NEGATIVE) {
                Cell_handle tmp = get<0>(_cur)->neighbor( _tr.cw(get<2>(_cur)) );
                _prev = _cur;
                get<0>(_cur) = tmp;
                get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex(get<2>(_cur)) );
                return;
            }

            // Then we try the opposite edge.
            Orientation op = _tr.coplanar_orientation( vert[_tr.ccw(get<2>(_cur))], vert[_tr.cw(get<2>(_cur))], vert[get<2>(_cur)], _t_vert );
            if( op == NEGATIVE) {
                Cell_handle tmp = get<0>(_cur)->neighbor(get<2>(_cur));
                get<0>(_prev) = get<0>(_cur);
                get<0>(_cur) = tmp;

                switch( ocw+occw ) {
                    case 2:
                        get<1>(_prev) = Tr::EDGE;
                        get<2>(_prev) = _tr.ccw( get<2>(_cur) );
                        get<3>(_prev) = _tr.cw( get<2>(_cur) );
                        get<1>(_cur) = Tr::EDGE;
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<2>(_prev) ) );
                        get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( get<3>(_prev) ) );
                        return;
                    case 1:
                        get<1>(_prev) = Tr::VERTEX;
                        get<1>(_cur) = Tr::VERTEX;
                        if( ocw == COLLINEAR ) get<2>(_prev) = _tr.cw( get<2>(_cur) );
                        else get<2>(_cur) = _tr.ccw( get<2>(_cur) );
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
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr.ccw(get<2>(_cur)), -1 );
                else if( occw == 0 )
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr.cw(get<2>(_cur)), -1 );
                else
                    _prev = Simplex( get<0>(_cur), Tr::EDGE, get<2>(_cur), -1 );
                break;
            case 1:
                if( ocw == 1 )
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.ccw(get<2>(_cur)), -1 );
                else if( occw == 1 )
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.cw(get<2>(_cur)), -1 );
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
                switch( _tr.coplanar_orientation( vert[get<2>(_cur)], vert[get<3>(_cur)], vert[lk], _t_vert ) ) {
                    _prev = _cur;
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

            Orientation o = _tr.coplanar_orientation( _s_vert, vert[lk], vert[get<2>(_cur)], _t_vert );
            Orientation op;
            switch( o ) {
                case POSITIVE: {
                    // The ray passes through the edge ik.
                    op = _tr.coplanar_orientation( vert[lk], vert[get<2>(_cur)], _s_vert, _t_vert );
                    if( op == NEGATIVE ) {
                        Cell_handle tmp = get<0>(_cur)->neighbor(get<3>(_cur));
                        get<0>(_prev) = get<0>(_cur);
                        get<0>(_cur) = tmp;

                        if( _tr.collinear( _s_vert, vert[get<2>(_cur)], _t_vert ) ) {
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
                    op = _tr.coplanar_orientation( vert[lk], vert[get<3>(_cur)], _s_vert, _t_vert );
                    if( op == NEGATIVE ) {
                        Cell_handle tmp = get<0>(_cur)->neighbor(get<2>(_cur));
                        get<0>(_prev) = get<0>(_cur);
                        get<0>(_cur) = tmp;

                        if( _tr.collinear( _s_vert, vert[get<3>(_cur)], _t_vert ) ) {
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
            // We test its edges in a random order until we find a neighbor to go further
            int li = rng.get_int(0, 3);

            Orientation o[3];
            bool calc[3] = { false, false, false };

            for( int j = 0; j != 3; ++j, li = _tr.ccw(li) ) {
                Cell_handle next = get<0>(_cur)->neighbor(li);
                if( next == get<0>(_prev) )
                    continue;

                // The target should lie on the other side of the edge.
                Orientation op = _tr.coplanar_orientation( vert[_tr.ccw(li)], vert[_tr.cw(li)], vert[li], _t_vert );
                if( op == POSITIVE )
                    continue;

                // The target should lie inside the wedge.
                if( !calc[_tr.ccw(li)] ) {
                    o[_tr.ccw(li)] = _tr.coplanar_orientation( _s_vert, vert[_tr.ccw(li)], vert[_tr.cw(li)], _t_vert );
                    calc[_tr.ccw(li)] = true;
                }
                if( o[_tr.ccw(li)] == NEGATIVE )
                    continue;
                else if( op == COLLINEAR && o[_tr.ccw(li)] == COLLINEAR ) {
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.ccw(li), -1 );
                    get<0>(_cur) = Cell_handle();
                    return;
                }

                if( !calc[_tr.cw(li)] ) {
                    o[_tr.cw(li)] = _tr.coplanar_orientation( _s_vert, vert[_tr.cw(li)], vert[li], _t_vert );
                    calc[_tr.cw(li)] = true;
                }
                if( o[_tr.cw(li)] == POSITIVE )
                    continue;
                else if( op == COLLINEAR && o[_tr.cw(li)] == COLLINEAR ) {
                    _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.cw(li), -1 );
                    get<0>(_cur) = Cell_handle();
                    return;
                }

                get<0>(_prev) = get<0>(_cur);
                get<0>(_cur) = next;

                switch( o[_tr.ccw(li)] + o[_tr.cw(li)] ) {
                    case 2:
                        get<1>(_prev) = Tr::EDGE;
                        get<2>(_prev) = _tr.ccw(li);
                        get<3>(_prev) = _tr.cw(li);
                        get<1>(_cur) = Tr::EDGE;
                        get<2>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( _tr.ccw(li) ) );
                        get<3>(_cur) = get<0>(_cur)->index( get<0>(_prev)->vertex( _tr.cw(li) ) );
                        return;
                    case 1:
                        get<1>(_prev) = Tr::VERTEX;
                        get<1>(_cur) = Tr::VERTEX;
                        if( o[_tr.ccw(li)] == COLLINEAR ) get<2>(_prev) = _tr.ccw(li);
                        else get<2>(_prev) = _tr.cw(li);
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
    CGAL_triangulation_precondition( _tr.is_infinite( get<0>(_cur)->vertex(3) ) );
    CGAL_triangulation_precondition( _tr.is_infinite( get<0>(_cur)->vertex(inf) ) );
	
    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = get<0>(_cur)->neighbor(inf);
    if (fin == get<0>(_prev)) {
        _prev = _cur;
        get<0>(_cur) = Cell_handle();
        return;
    }

    // Check the neighboring cells.
    Orientation occw = _tr.coplanar_orientation( _s_vert, get<0>(_cur)->vertex( _tr.ccw(inf)), get<0>(_cur)->vertex(_tr.cw(inf)), _t_vert );
    if( occw == NEGATIVE ) {
        Cell_handle tmp = get<0>(_cur)->neighbor(_tr.cw(inf));
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr.ccw(inf), inf );
        _cur = Simplex( tmp, Tr::EDGE, tmp->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), tmp->index( get<0>(_prev)->vertex( get<3>(_prev) ) ) );
        return;
    }
    Orientation ocw = _tr.coplanar_orientation( _s_vert, get<0>(_cur)->vertex( _tr.cw(inf)), get<0>(_cur)->vertex(_tr.ccw(inf)), _t_vert );
    if( ocw == NEGATIVE ) {
        Cell_handle tmp = get<0>(_cur)->neighbor(_tr.ccw(inf));
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr.cw(inf), inf );
        _cur = Simplex( tmp, Tr::EDGE, tmp->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), tmp->index( get<0>(_prev)->vertex( get<3>(_prev) ) ) );
        return;
    }
    Orientation op = _tr.coplanar_orientation( get<0>(_cur)->vertex( _tr.ccw(inf) ), get<0>(_cur)->vertex( _tr.cw(inf) ), _s_vert, _t_vert );
    switch( op ) {
    case NEGATIVE:
        if( occw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.ccw(inf), -1 );
            _cur = Simplex( fin, Tr::VERTEX, fin->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), -1 );
            return;
        }
        if( ocw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.cw(inf), -1 );
            _cur = Simplex( fin, Tr::VERTEX, fin->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), -1 );
            return;
        }
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr.ccw(inf), _tr.cw(inf) );
        _cur = Simplex( fin, Tr::EDGE, fin->index( get<0>(_prev)->vertex( get<2>(_prev) ) ), fin->index( get<0>(_prev)->vertex( get<3>(_prev) ) ) );
        return;
    case COLLINEAR:
        if( occw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.ccw(inf), -1 );
            get<0>(_cur) = Cell_handle();
            return;
        }
        if( ocw == COLLINEAR ) {
            _prev = Simplex( get<0>(_cur), Tr::VERTEX, _tr.cw(inf), -1 );
            get<0>(_cur) = Cell_handle();
            return;
        }
        _prev = Simplex( get<0>(_cur), Tr::EDGE, _tr.ccw(inf), _tr.cw(inf) );
        get<0>(_cur) = Cell_handle();
        return;
    case POSITIVE:
        // The target lies in this infinite cell.
        _prev = Simplex( get<0>(_cur), Tr::OUTSIDE_CONVEX_HULL, -1, -1 );
        get<0>(_cur) = Cell_handle();
        return;
    }
}
	
} //end of CGAL namespace

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H
