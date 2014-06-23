//A method to construct a surface.
//Copyright (C) 2014  INRIA - Sophia Antipolis
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
// Author(s):      Thijs van Lankveld

#ifndef CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_IMPL_H
#define CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_IMPL_H

#include <omp.h>

#include <CGAL/Scale_space_reconstruction_3/Weighted_PCA_projection_3.h>

namespace CGAL {
    
template < class GeomTraits, class FixedScale, class Shells >
class
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
Finite_point_iterator: public Finite_vertices_iterator {
    typedef Finite_vertices_iterator    Base;
    typedef Finite_point_iterator       Self;

public:
    Finite_point_iterator(): Base() {}
    Finite_point_iterator( const Base& b ): Base(b) {}

    operator Point&() const { return this->point(); }
}; // class Finite_point_iterator


template < class GeomTraits, class FixedScale, class Shells >
class
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
In_surface_tester {
    const Shape* _s;

public:
    In_surface_tester() {}
    In_surface_tester( const Shape* s ): _s(s) {}

    inline bool operator()( const Cell_iterator& c ) const {
        return false;
    }

    inline bool operator()( const Edge_iterator& e ) const {
        typename Shape::Classification_type cl = _s->classify( e );
        return cl == Shape::REGULAR || cl == Shape::SINGULAR;
    }

    inline bool operator()( const Facet_iterator& f ) const {
        typename Shape::Classification_type cl = _s->classify( f );
        return cl == Shape::REGULAR || cl == Shape::SINGULAR;
    }

    inline bool operator()( const Vertex_iterator& v ) const {
        return !_s->is_infinite( v );
    }
}; // In_surface_tester


template < class GeomTraits, class FixedScale, class Shells >
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
Scale_space_surface_reconstruction_3( unsigned int neighbors, unsigned int samples )
: _mean_neighbors(neighbors), _samples(samples), _squared_radius(-1), _shape(0) {}

template < class GeomTraits, class FixedScale, class Shells >
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
Scale_space_surface_reconstruction_3( FT sq_radius )
: _mean_neighbors(30), _samples(200), _squared_radius(sq_radius), _shape(0) {}

template < class GeomTraits, class FixedScale, class Shells >
inline bool
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
is_handled( Cell_handle c, unsigned int li ) const {
    switch( li ) {
    case 0: return ( c->info()&1 ) != 0;
    case 1: return ( c->info()&2 ) != 0;
    case 2: return ( c->info()&4 ) != 0;
    case 3: return ( c->info()&8 ) != 0;
    }
    return false;
}

template < class GeomTraits, class FixedScale, class Shells >
inline void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
mark_handled( Cell_handle c, unsigned int li ) {
    switch( li ) {
    case 0: c->info() |= 1; return;
    case 1: c->info() |= 2; return;
    case 2: c->info() |= 4; return;
    case 3: c->info() |= 8; return;
    }
}

template < class GeomTraits, class FixedScale, class Shells >
inline typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::Triple
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
ordered_facet_indices( const Facet& f ) const {
    if( (f.second&1) == 0 )
        return make_array<unsigned int>( f.first->vertex( (f.second+2)&3 )->info(),
                                         f.first->vertex( (f.second+1)&3 )->info(),
                                         f.first->vertex( (f.second+3)&3 )->info() );
    else
        return make_array<unsigned int>( f.first->vertex( (f.second+1)&3 )->info(),
                                         f.first->vertex( (f.second+2)&3 )->info(),
                                         f.first->vertex( (f.second+3)&3 )->info() );
}

template < class GeomTraits, class FixedScale, class Shells >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
collect_shell( Cell_handle c, unsigned int li ) {
    // Collect one surface mesh from the alpha-shape in a fashion similar to ball-pivoting.
    // Invariant: the facet is regular or singular.

    // To stop stack overflows: use own stack.
    std::stack<Facet> stack;
    stack.push( Facet(c, li) );

    Facet f;
    Cell_handle n, p;
    int ni, pi;
    Vertex_handle a;
    Classification_type cl;
    bool processed = false;
    while( !stack.empty() ) {
        f = stack.top();
        stack.pop();

        // Check if the cell was already handled.
        // Note that this is an extra check that in many cases is not necessary.
        if( is_handled(f) )
            continue;

        // The facet is part of the surface.
        CGAL_assertion( !_shape->is_infinite(f) );
        mark_handled(f);

        // Output the facet as a triple of indices.
        _surface.push_back( ordered_facet_indices(f) );
        if( !processed ) {
            _shells.push_back( --_surface.end() );
            processed = true;
        }
		
        // Pivot over each of the facet's edges and continue the surface at the next regular or singular facet.
        for( unsigned int i = 0; i < 4; ++i ) {
            // Skip the current facet.
            if( i == f.second || is_handled(f.first, i) )
                continue;

            // Rotate around the edge (starting from the shared facet in the current cell) until a regular or singular facet is found.
            n = f.first;
            ni = i;
            a = f.first->vertex( f.second );
            cl = _shape->classify( Facet(n, ni) );
            while( cl != Shape::REGULAR && cl != Shape::SINGULAR ) {
                p = n;
                n = n->neighbor(ni);
                ni = n->index(a);
                pi = n->index(p);
                a = n->vertex(pi);
                cl = _shape->classify( Facet(n, ni) );
            }

            // Continue the surface at the next regular or singular facet.
            stack.push( Facet(n, ni) );
        }
    }
}

template < class GeomTraits, class FixedScale, class Shells >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
collect_facets( Tag_true ) {
    // Collect all surface meshes from the alpha-shape in a fashion similar to ball-pivoting.
    // Reset the facet handled markers.
    for( All_cells_iterator cit = _shape->all_cells_begin(); cit != _shape->all_cells_end(); ++cit )
        cit->info() = 0;

    // We check each of the facets: if it is not handled and either regular or singular,
    // we start collecting the next surface from here.
    for( Finite_facets_iterator fit = _shape->finite_facets_begin(); fit != _shape->finite_facets_end(); ++fit ) {
        switch( _shape->classify( *fit ) ) {
            case Shape::REGULAR:
                // Build a surface from the outer cell.
                if( _shape->classify(fit->first) == Shape::EXTERIOR )
                    collect_shell( *fit );
                else
                    collect_shell( _shape->mirror_facet( *fit ) );
                break;
            case Shape::SINGULAR:
                // Build a surface from both incident cells.
                collect_shell( *fit );
                collect_shell( _shape->mirror_facet( *fit ) );
                break;
        }
    }
}

template < class GeomTraits, class FixedScale, class Shells >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
collect_facets( Tag_false ) {
    // Collect all facets from the alpha-shape in an unordered fashion.
    for( Finite_facets_iterator fit = _shape->finite_facets_begin(); fit != _shape->finite_facets_end(); ++fit ) {
        switch( _shape->classify( *fit ) ) {
            case Shape::REGULAR:
                // Collect the outer cell.
                if( _shape->classify(fit->first) == Shape::EXTERIOR )
                    _surface.push_back( ordered_facet_indices( *fit ) );
                else
                    _surface.push_back( ordered_facet_indices( _shape->mirror_facet( *fit ) ) );
                break;
            case Shape::SINGULAR:
                // Collect both incident cells.
                _surface.push_back( ordered_facet_indices( *fit ) );
                _surface.push_back( ordered_facet_indices( _shape->mirror_facet( *fit ) ) );
            break;
        }
    }
}

template < class GeomTraits, class FixedScale, class Shells >
const typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::Shape&
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
shape() const {
	if( !has_shape() )
        _shape = Shape_construction_3().construct( _shape, _squared_radius );
    return *_shape;
}

template < class GeomTraits, class FixedScale, class Shells >
typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::FT
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
estimate_neighborhood_radius( unsigned int neighbors, unsigned int samples ) {
    typename GeomTraits::Compute_squared_distance_3 squared_distance = GeomTraits().compute_squared_distance_3_object();

    unsigned int handled = 0;
    unsigned int checked = 0;
    FT radius = 0;
        
    if( !_tree.is_built() )
        _tree.build();

    for( Const_point_iterator it = _tree.begin(); it != _tree.end(); ++it ) {
        unsigned int left = (unsigned int)( _tree.size() - handled );
        if( samples >= left || _generator.get_double() < double(samples - checked) / left ) {
            // The neighborhood should contain the point itself as well.
            Static_search search( _tree, *it, neighbors+1 );
            radius += sqrt( to_double( squared_distance( *it, ( search.end()-1 )->first ) ) );
            ++checked;
        }
        ++handled;
    }
    radius /= double(checked);

    set_neighborhood_squared_radius( radius * radius );
    return radius;
}

// Doxygen has a bug where it cannot link the declaration and implementation
// of methods with a templated parameter.
#ifndef DOXYGEN_RUNNING
template < class GeomTraits, class FixedScale, class Shells >
template < class InputIterator >
typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::FT
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
estimate_neighborhood_radius( InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples ) {
    clear();
	insert_points( begin, end );
	return estimate_neighborhood_radius( neighbors, samples );
}
#endif // DOXYGEN_RUNNING

template < class GeomTraits, class FixedScale, class Shells >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
advance_scale_space( unsigned int iterations ) {
    typedef std::vector< unsigned int >					CountVec;
    typedef std::map<Point, size_t>			            PIMap;
    typedef std::vector<Point>                          Pointset;

    typedef Weighted_PCA_projection_3< GeomTraits >     WPCAP;

    typedef internal::_ENVIRONMENT::s_ptr_type			p_size_t;
		
    // This method must be called after filling the point collection.
    if( iterations == 0 || _tree.empty() ) return;
        
    if( !has_neighborhood_radius() )
        estimate_neighborhood_radius();

    // In order to reduce memory consumption, we maintain two data structures:
    // a local tree (because openMP cannot work on global members), and
    // a vector for the points after smoothing.
    Pointset points;
    points.assign( _tree.begin(), _tree.end() );
    _tree.clear();
    Search_tree tree;

    // Construct a search tree of the points.
    // Note that the tree has to be local for openMP.
    for( unsigned int iter = 0; iter < iterations; ++iter ) {
        tree.insert( points.begin(), points.end() );
        if( !tree.is_built() )
            tree.build();

        // Collect the number of neighbors of each point.
        // This can be done parallel.
        CountVec neighbors( tree.size(), 0 );
        typename GeomTraits::Compare_squared_distance_3 compare;
        p_size_t count = tree.size(); // openMP can only use signed variables.
        const FT squared_radius = _squared_radius; // openMP can only use local variables.
#ifdef _OPENMP
#pragma omp parallel for shared(count,tree,points,squared_radius,neighbors) firstprivate(compare)
#endif
        for( p_size_t i = 0; i < count; ++i ) {
            // Iterate over the neighbors until the first one is found that is too far.
            Dynamic_search search( tree, points[i] );
            for( typename Dynamic_search::iterator nit = search.begin(); nit != search.end() && compare( points[i], nit->first, squared_radius ) != LARGER; ++nit )
                ++neighbors[i];
        }
            
        // Construct a mapping from each point to its index.
        PIMap indices;
        p_size_t index = 0;
        for( typename Pointset::const_iterator pit = points.begin(); pit != points.end(); ++pit, ++index)
            indices[ *pit ] = index;

        // Compute the tranformed point locations.
        // This can be done parallel.
#ifdef _OPENMP
#pragma omp parallel for shared(count,neighbors,tree,squared_radius) firstprivate(compare)
#endif
        for( p_size_t i = 0; i < count; ++i ) {
            // If the neighborhood is too small, the vertex is not moved.
            if( neighbors[i] < 4 )
                continue;

            // Collect the vertices within the ball and their weights.
            Dynamic_search search( tree, points[i] );
            WPCAP pca( neighbors[i] );
            unsigned int column = 0;
            for( typename Dynamic_search::iterator nit = search.begin(); nit != search.end() && column < neighbors[i]; ++nit, ++column ) {
                pca.set_point( column, nit->first, 1.0 / neighbors[ indices[nit->first] ] );
            }
            CGAL_assertion( column == neighbors[i] );

            // Compute the weighted least-squares planar approximation of the point set.
            if( !pca.approximate() )
                continue;

            // The vertex is moved by projecting it onto the plane
            // through the barycenter and orthogonal to the Eigen vector with smallest Eigen value.
            points[i] = pca.project( points[i] );
        }
            
        tree.clear();
    }

    _tree.insert( points.begin(), points.end() );
}

template < class GeomTraits, class FixedScale, class Shells >
template < class InputIterator >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
construct_shape( InputIterator begin, InputIterator end ) {
    deinit_shape();
    if( !has_neighborhood_radius() )
        estimate_neighborhood_radius();
    _shape = Shape_construction_3().construct( begin, end, _squared_radius );
}

template < class GeomTraits, class FixedScale, class Shells >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
collect_surface() {
    _surface.clear();
    if( !has_shape() )
        construct_shape();
    collect_facets();
}

template < class GeomTraits, class FixedScale, class Shells >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
reconstruct_surface( unsigned int iterations ) {
    // Smooth the scale space.
    advance_scale_space( iterations );

    // Mesh the perturbed points.
    collect_surface();
}

template < class GeomTraits, class FixedScale, class Shells >
template < class InputIterator >
void
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::
reconstruct_surface( InputIterator begin, InputIterator end, unsigned int iterations ) {
    // Compute the radius for which the mean ball would contain the required number of neighbors.
    clear();
    insert_points( begin, end );

    // Smooth the scale space.
    advance_scale_space( iterations );

    // Mesh the perturbed points.
    collect_surface();
}

template < class GeomTraits, class FixedScale, class Shells >
typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::Const_triple_iterator
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::shell_begin( std::size_t shell ) const {
    CGAL_assertion( Shells::value == true );
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    return _shells[ shell ];
}

template < class GeomTraits, class FixedScale, class Shells >
typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::Triple_iterator
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::shell_begin( std::size_t shell ) {
    CGAL_assertion( Shells::value == true );
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    return _shells[ shell ];
}

template < class GeomTraits, class FixedScale, class Shells >
typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::Const_triple_iterator
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::shell_end( std::size_t shell ) const {
    CGAL_assertion( Shells::value == true );
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    if( shell == _shells.size()-1 )
        return _surface.end();
    return _shells[ shell+1 ];
}

template < class GeomTraits, class FixedScale, class Shells >
typename Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::Triple_iterator
Scale_space_surface_reconstruction_3<GeomTraits,FixedScale,Shells>::shell_end( std::size_t shell ) {
    CGAL_assertion( Shells::value == true );
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    if( shell == _shells.size()-1 )
        return _surface.end();
    return _shells[ shell+1 ];
}

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_IMPL_H