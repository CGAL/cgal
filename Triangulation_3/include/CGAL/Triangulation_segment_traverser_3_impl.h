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

template <class Tr>
Triangulation_segment_traverser_3<Tr>::
Triangulation_segment_traverser_3( const Tr& tr, Vertex_handle s, const Point& t )
: _tr(tr), _pos(), _prev(), _lj(-1), _lt(Tr::VERTEX), _done(false) {
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

    _pos = s->cell();
    int inf;
    if( _pos->has_vertex( _tr.infinite_vertex(), inf ) )
        _pos = _pos->neighbor(inf);
    _li = _pos->index(s);
	
    CGAL_triangulation_postcondition( _pos != Cell_handle() );
}

template < class Tr>
Triangulation_segment_traverser_3<Tr>::
Triangulation_segment_traverser_3( const Tr& tr, Vertex_handle s, Vertex_handle t )
: _tr(tr), _pos(), _prev(), _lj(-1), _lt(Tr::VERTEX), _done(false) {
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

    _pos = s->cell();
    int inf;
    if( _pos->has_vertex( _tr.infinite_vertex(), inf ) )
        _pos = _pos->neighbor(inf);
    _li = _pos->index(s);	
}


template <class Tr>
Triangulation_segment_traverser_3<Tr>::
Triangulation_segment_traverser_3( const Tr& tr, const Point& s, const Point& t, Cell_handle hint )
: _tr(tr), _pos(), _prev(), _li(-1), _lj(-1), _done(false) {
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

    _pos = _tr.locate( s, _lt, _li, _lj, hint );

    CGAL_triangulation_postcondition( _pos != Cell_handle() );
}

template <class Tr>
Triangulation_segment_traverser_3<Tr>::
Triangulation_segment_traverser_3( const Tr& tr, const Segment& s, Cell_handle hint )
: Triangulation_segment_traverser_3<Tr>( tr, s.source(), s.target(), hint ) {}

template <class Tr>
inline Triangulation_segment_traverser_3<Tr>&
Triangulation_segment_traverser_3<Tr>::operator++() {
    CGAL_triangulation_precondition( _pos != Cell_handle() );
    increment();
    return *this;
}

template <class Tr>
inline Triangulation_segment_traverser_3<Tr>
Triangulation_segment_traverser_3<Tr>::operator++( int ) {
    TST tmp( *this );
    ++( *this );
    return tmp;
}

template <class Tr>
inline typename Triangulation_segment_traverser_3<Tr>::Cell_handle
Triangulation_segment_traverser_3<Tr>::traverse() {
    while( !_done )
        increment();
    return _pos;
}

template <class Tr>
inline bool Triangulation_segment_traverser_3<Tr>::
operator==( const TST& ct ) const {
    CGAL_triangulation_precondition( _pos != Cell_handle() );
    CGAL_triangulation_precondition( ct._pos != Cell_handle() );
    return ( _pos == ct._pos &&
             _prev == ct._prev &&
             _tr == ct._tr &&
             _li == ct._li &&
             _lj == ct._lj &&
             _lt == ct._lt &&
             _source == ct._source &&
             _target == ct._target &&
             _s_vertex == ct._s_vertex &&
             _t_vertex == ct._t_vertex &&
             _done == ct._done );
}

template <class Tr>
inline bool Triangulation_segment_traverser_3<Tr>::
operator!=( const TST& ct ) const {
    return !( *this == ct );
}

template <class Tr>
inline bool Triangulation_segment_traverser_3<Tr>::
operator==( Nullptr_t CGAL_triangulation_assertion_code(n) ) const {
    CGAL_triangulation_assertion( n == NULL );
    return _pos == Cell_handle();
}

template <class Tr>
inline bool Triangulation_segment_traverser_3<Tr>::
operator!=( Nullptr_t n ) const {
    return !( *this == n );
}

template <class Tr>
void Triangulation_segment_traverser_3<Tr>::
increment() {
    CGAL_triangulation_precondition( !_done );
    // Walks to the next cell over a facet intersected by the line from source to target.
    // This method is based on Triangulation_3::locate().
    int inf;
    switch( _tr.dimension() ) {
        case 3: {
            // Infinite cells should be handled differently.
            if( _pos->has_vertex( _tr.infinite_vertex(), inf ) )
                increment_3_inf( inf );
            else
                increment_3();
            break;
        }
        case 2: {
            if( _pos->has_vertex( _tr.infinite_vertex(), inf ) )
                increment_2_inf( inf );
            else
                increment_2();
            break;
        }
    }
}

template <class Tr>
void Triangulation_segment_traverser_3<Tr>::
increment_3() {
    Vertex_handle vert[4] = { _pos->vertex(0),
                              _pos->vertex(1),
                              _pos->vertex(2),
                              _pos->vertex(3) };

    // We check in which direction the target lies
    // by comparing its position relative to the planes through the
    // source and the edges of the cell.
    Orientation o[6];
    // We keep track of which orientations are calculated.
    bool calc[6] = { false, false, false, false, false, false };

    if( _lt == Tr::VERTEX ) {
        // The three planes through the vertex are set to coplanar.
        for( int j = 0; j < 4; ++j ) {
            if( _li != j ) {
                int ij = edgeIndex( _li, j );
                o[ij] = COPLANAR;
                calc[ij] = true;
            }
        }
    }
    else if( _lt == Tr::EDGE ) {
        // The plane through the edge is set to coplanar.
        int ij = edgeIndex( _li, _lj );
        o[ij] = COPLANAR;
        calc[ij] = true;
    }

    // For the remembering stochastic walk, we start trying with a random facet.
    int li = rng.template get_bits<2>();

    CGAL_triangulation_assertion_code( bool incell = true; )
    for( int k = 0; k < 4; ++k, li = _tr.increment_index(li) ) {
        Cell_handle next = _pos->neighbor(li);
        if( next == _prev )
            continue;
  
        if( _t_vertex == _pos->vertex(0) || _t_vertex == _pos->vertex(1) ||
            _t_vertex == _pos->vertex(2) || _t_vertex == _pos->vertex(3) ) {
            // The target is inside the cell.
            _done = true;
            return;
        }
  
        // Check if the target is outside the cell.
        Vertex_handle backup = vert[li];
        vert[li] = _t_vert;
       
        if( _tr.orientation( vert[0], vert[1], vert[2], vert[3] ) != NEGATIVE ) {
            vert[li] = backup;
            continue;
        }
        CGAL_triangulation_assertion_code( incell = false; )

        // Check if the target is inside the pyramid.
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

            Or -= o[oij];
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
        }

        if( Or == 0 ) {
            // Either the target is not inside the pyramid,
            // or the pyramid is degenerate.
            vert[li] = backup;
            continue;
        }

        // The target is inside the pyramid.
        _prev = _pos;
        _pos = next;

        switch( Or ) {
            case 3:
                _lt = Tr::FACET;
                _li = _pos->index(_prev);
                return;
            case 2:
                _lt = Tr::EDGE;
                for( int j = 0; j < 4; ++j ) {
                    if( li != j && o[ 5 - edgeIndex(li, j) ] == COPLANAR) {
                        Edge opp = _tr.opposite_edge( _prev, li, j );
                        _li = _pos->index( _prev->vertex(opp.second) );
                        _lj = _pos->index( _prev->vertex(opp.third) );
                        return;
                    }
                }
                CGAL_triangulation_assertion( false );
                return;
            case 1:
                _lt = Tr::VERTEX;
                for( int j = 0; j < 4; ++j ) {
                    if( li != j && o[ 5 - edgeIndex(li, j) ] == NEGATIVE ) {
                        _li = _pos->index( _prev->vertex(j) );
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
    _done = true;
    return;
}

template <class Tr>
void Triangulation_segment_traverser_3<Tr>::
increment_3_inf( int inf ) {
    CGAL_triangulation_precondition( _tr.is_infinite( _pos->vertex(inf) ) );

    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = _pos->neighbor(inf);
    if( fin == _prev ) {
        _done = true;
        // Here, I do not change the _lt etc. because they were already set.
        return;
    }

    Vertex_handle vert[4];
    for( int i = 0; i != 4; ++i )
        if( i != inf )
            vert[i] = _pos->vertex(i);
    vert[inf] = _t_vert;
    Orientation o[4];

    // For the remembering stochastic walk, we start trying with a random index:
    int li = rng.template get_bits<2>();

    // Check if the target lies outside the convex hull.
    if( _tr.orientation( vert[0], vert[1], vert[2], vert[3] ) == POSITIVE ) {
        // The target lies in an infinite cell.
        // Note that we do not traverse to other infinite cells.
        _done = true;
        return;
    }

    vert[inf] = _s_vert;
    CGAL_triangulation_assertion( _tr.orientation( vert[0], vert[1], vert[2], vert[3] ) == POSITIVE );

    // Check if the line enters an adjacent infinite cell.
    // This occurs if the target lies on the other side of
    // a plane through one of the finite edges and the source point.
    for( int j = 0; j != 4; ++j, li = _tr.increment_index(li) ) {
        if( li == inf ) {
            o[li] = COPLANAR;
            continue;
        }

        Cell_handle next = _pos->neighbor(li);
        if( next == _prev ) {
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
        _prev = _pos;
        _pos = next;
        CGAL_triangulation_assertion( _tr.is_infinite(_pos) );

        _lt = Tr::FACET;
        _li = _pos->index(_prev);
        return;
    }

    // The line enters the convex hull here (or lies on the finite facet).
    _prev = _pos;
    _pos = fin;

    // Check through which simplex the line traverses.
    switch( o[0]+o[1]+o[2]+o[3] ) {
        case 3:
            _lt = Tr::FACET;
            _li = _pos->index(_prev);
            return;
        case 2:
            _lt = Tr::EDGE;
            for( int i = 0; i < 4; ++i ) {
                if( o[i] == COPLANAR && i != inf ) {
                    Edge opp = _tr.opposite_edge( _prev, inf, i );
                    _li = _pos->index( _prev->vertex(opp.second) );
                    _lj = _pos->index( _prev->vertex(opp.third) );
                    return;
                }
            }
            CGAL_triangulation_assertion( false );
            return;
        case 1:
            _lt = Tr::VERTEX;
            for( int i = 0; i < 4; ++i ) {
                if( o[i] == POSITIVE ) {
                    _li = _pos->index( _prev->vertex(i) );
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

template <class Tr>
void Triangulation_segment_traverser_3<Tr>::
increment_2() {
    Vertex_handle vert[3] = { _pos->vertex(0),
                              _pos->vertex(1),
                              _pos->vertex(2) };
			
    switch( _lt ) {
        case Tr::VERTEX: {
            // First we try the incident edges.
            Orientation ocw = _tr.coplanar_orientation( vert[_li], vert[_tr.cw(_li)], vert[_tr.ccw(_li)], _t_vert );
            if( _pos->neighbor( _tr.ccw(_li) ) != _prev && ocw == NEGATIVE) {
                Cell_handle tmp = _pos->neighbor( _tr.ccw(_li) );
                _prev = _pos;
                _pos = tmp;
                _li = _pos->index( _prev->vertex(_li) );
                return;
            }
            Orientation occw = _tr.coplanar_orientation( vert[_li], vert[_tr.ccw(_li)], vert[_tr.cw(_li)], _t_vert );
            if( _pos->neighbor( _tr.cw(_li) ) != _prev && occw == NEGATIVE) {
                Cell_handle tmp = _pos->neighbor( _tr.cw(_li) );
                _prev = _pos;
                _pos = tmp;
                _li = _pos->index( _prev->vertex(_li) );
                return;
            }

            // Then we try the opposite edge.
            if( _tr.coplanar_orientation( vert[_tr.ccw(_li)], vert[_tr.cw(_li)], vert[_li], _t_vert ) == NEGATIVE) {
                Cell_handle tmp = _pos->neighbor(_li);
                _prev = _pos;
                _pos = tmp;

                switch( ocw+occw ) {
                    case 2:
                        _lt = Tr::EDGE;
                        _li = _pos->index( _prev->vertex( _tr.ccw(_li) ) );
                        _lj = _pos->index( _prev->vertex( _tr.cw(_li) ) );
                        return;
                    case 1:
                        _lt = Tr::VERTEX;
                        if( ocw == COLLINEAR ) _li = _pos->index( _prev->vertex( _tr.cw(_li) ) );
                        else _li = _pos->index( _prev->vertex( _tr.ccw(_li) ) );
                        return;
                    default:
                        CGAL_triangulation_assertion(false);
                        return;
                }
            }

            // The target lies in this cell.
            _done = true;
            return;
        }
        case Tr::EDGE: {
            int lk = 3-_li-_lj;

            if( _pos->neighbor(lk) != _prev ) {
                // Check the edge itself
                switch( _tr.coplanar_orientation( vert[_li], vert[_lj], vert[lk], _t_vert ) ) {
                    case COLLINEAR:
                        // The target lies in this cell.
                        _done = true;
                        return;
                    case NEGATIVE: {
                        // The target lies opposite of the edge.
                        Cell_handle tmp = _pos->neighbor(lk);
                        _prev = _pos;
                        _pos = tmp;
                        _li = _pos->index( _prev->vertex(_li) );
                        _lj = _pos->index( _prev->vertex(_lj) );
                        return;
                    }
                    default:
                        break;
                }
            }

            Orientation o = _tr.coplanar_orientation( _s_vert, vert[lk], vert[_li], _t_vert );
            switch( o ) {
                case POSITIVE: {
                    // The ray passes through the edge ik.
                    if( _tr.coplanar_orientation( vert[lk], vert[_li], _s_vert, _t_vert ) == NEGATIVE ) {
                        Cell_handle tmp = _pos->neighbor(_lj);
                        _prev = _pos;
                        _pos = tmp;

                        if( _tr.collinear( _s_vert, vert[_li], _t_vert ) ) {
                            _lt = Tr::VERTEX;
                            _li = _pos->index( _prev->vertex(_li) );
                        }
                        else {
                            _lt = Tr::EDGE;
                            _li = _pos->index( _prev->vertex(_li) );
                            _lj = _pos->index( _prev->vertex(lk) );
                        }
                        return;
                    }
                    break;
                }
                default: {
                    // The ray passes through the edge jk.
                    if( _tr.coplanar_orientation( vert[lk], vert[_lj], _s_vert, _t_vert ) == NEGATIVE ) {
                        Cell_handle tmp = _pos->neighbor(_li);
                        _prev = _pos;
                        _pos = tmp;

                        if( _tr.collinear( _s_vert, vert[_lj], _t_vert ) ) {
                            _lt = Tr::VERTEX;
                            _li = _pos->index( _prev->vertex(_lj) );
                        }
                        else if( o == COLLINEAR ) {
                            _lt = Tr::VERTEX;
                            _li = _pos->index( _prev->vertex(lk) );
                        }
                        else {
                            _lt = Tr::EDGE;
                            _li = _pos->index( _prev->vertex(lk) );
                            _lj = _pos->index( _prev->vertex(_lj) );
                        }
                        return;
                    }
                    break;
                }
            }

            // The target lies in this cell.
            _done = true;
            return;
        }
        case Tr::FACET: {
            // We test its edges in a random order until we find a neighbor to go further
            int li = rng.get_int(0, 3);

            Orientation o[3];
            bool calc[3] = { false, false, false };

            for( int j = 0; j != 3; ++j, li = _tr.ccw(li) ) {
                Cell_handle next = _pos->neighbor(li);
                if( next == _prev )
                    continue;

                // The target should lie on the other side of the edge.
                if( _tr.coplanar_orientation( vert[_tr.ccw(li)], vert[_tr.cw(li)], vert[li], _t_vert ) != NEGATIVE )
                    continue;

                // The target should lie inside the wedge.
                if( !calc[_tr.ccw(li)] ) {
                    o[_tr.ccw(li)] = _tr.coplanar_orientation( _s_vert, vert[_tr.ccw(li)], vert[_tr.cw(li)], _t_vert );
                    calc[_tr.ccw(li)] = true;
                }
                if( o[_tr.ccw(li)] == NEGATIVE ) continue;

                if( !calc[_tr.cw(li)] ) {
                    o[_tr.cw(li)] = _tr.coplanar_orientation( _s_vert, vert[_tr.cw(li)], vert[_tr.ccw(li)], _t_vert );
                    calc[_tr.cw(li)] = true;
                }
                if( o[_tr.cw(li)] == POSITIVE ) continue;

                _prev = _pos;
                _pos = next;

                switch( o[_tr.ccw(li)] + o[_tr.cw(li)] ) {
                    case 2:
                        _lt = Tr::EDGE;
                        _li = _pos->index( _prev->vertex( _tr.ccw(li) ) );
                        _lj = _pos->index( _prev->vertex( _tr.cw(li) ) );
                        return;
                    case 1:
                        _lt = Tr::VERTEX;
                        if( o[_tr.ccw(li)] == COLLINEAR ) _li = _pos->index( _prev->vertex( _tr.ccw(li) ) );
                        else _li = _pos->index( _prev->vertex( _tr.cw(li) ) );
                        return;
                    default:
                        CGAL_triangulation_assertion( false );
                        return;
                }
            }

            // The target lies in this cell.
            _done = true;
            return;
        }
        default:
        CGAL_triangulation_assertion( false );
    }
}

template <class Tr>
void Triangulation_segment_traverser_3<Tr>::
increment_2_inf( int inf ) {
    CGAL_triangulation_precondition( _tr.is_infinite( _pos->vertex(3) ) );
    CGAL_triangulation_precondition( _tr.is_infinite( _pos->vertex(inf) ) );
	
    // If this cell was reached by traversal from a finite one, it must be the final cell.
    Cell_handle fin = _pos->neighbor(inf);
    if (fin == _prev) {
        _done = true;
        return;
    }

    // Check the neighboring cells.
    Gt::Coplanar_orientation_3 coplanar_orientation = Gt().coplanar_orientation_3_object();
    if( _tr.coplanar_orientation( _s_vert, _pos->vertex( _tr.ccw(inf)), _pos->vertex(_tr.cw(inf)), _t_vert ) == NEGATIVE ) {
        Cell_handle tmp = _pos->neighbor(_tr.cw(inf));
        _prev = _pos;
        _pos = tmp;
        return;
    }
    if( _tr.coplanar_orientation( _s_vert, _pos->vertex( _tr.cw(inf)), _pos->vertex(_tr.ccw(inf)), _t_vert ) == NEGATIVE ) {
        Cell_handle tmp = _pos->neighbor(_tr.ccw(inf));
        _prev = _pos;
        _pos = tmp;
        return;
    }
    if( _tr.coplanar_orientation( _pos->vertex( _tr.ccw(inf)), _pos->vertex(_tr.cw(inf)), _s_vert, _t_vert ) != POSITIVE ) {
        _prev = _pos;
        _pos = fin;
        return;
    }

    // The target lies in this infinite cell.
    _done = true;
    return;
}
	
} //end of CGAL namespace

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_IMPL_H
