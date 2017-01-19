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

#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_IMPL_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_IMPL_H

#include <CGAL/license/Scale_space_reconstruction_3.h>


//#include <omp.h>
#ifdef CGAL_LINKED_WITH_TBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#endif // CGAL_LINKED_WITH_TBB

#include <sstream>
#include <fstream>

#include <boost/function_output_iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <string>



namespace CGAL {

  
template < class Gt, class FS, class wA, class Ct >
class
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
Finite_point_iterator: public Finite_vertices_iterator {
    typedef Finite_vertices_iterator    Base;
    typedef Finite_point_iterator       Self;

public:
    Finite_point_iterator(): Base() {}
    Finite_point_iterator( const Base& b ): Base(b) {}

    operator Point&() const { return this->point(); }
}; // class Finite_point_iterator


template < class Gt, class FS, class wA, class Ct >
class
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
In_surface_tester {
    const Shape* _s;

public:
    In_surface_tester() {}
    In_surface_tester( const Shape* s ): _s(s) {}

    inline bool operator()( const Cell_iterator& /* c */ ) const {
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


namespace internal {
  namespace Scale_space {

    struct Inc {
      unsigned int * i;
      
      Inc(unsigned int& i)
        : i(&i)
      {}
      
      template <typename T>
      void operator()(const T&) const
      {
        ++(*i);
      }
      
    };

    template <typename T>
    struct operator_less
    {
      bool operator() (const T& a, const T& b) const
      {
	return &*a < &*b;
      }
    };

  }
}

// Compute the number of neighbors of a point that lie within a fixed radius.
template < class Gt, class FS, class wA, class Ct >
class
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
ComputeNN {
public:
    typedef std::vector< Point >        Pointset;
    typedef std::vector< unsigned int >	CountVec;

private:
    typename Gt::Compare_squared_distance_3 compare;

    const Pointset&     _pts;
    const Search_tree&  _tree;
    const FT           _sq_rd;
    CountVec&           _nn;

public:
    ComputeNN(const Pointset& points, const Search_tree&  tree,
              const FT& sq_radius, CountVec& nn)
      : _pts(points), _tree(tree), _sq_rd(sqrt(sq_radius)), _nn(nn) {}
    
#ifdef CGAL_LINKED_WITH_TBB
    void operator()( const tbb::blocked_range< std::size_t >& range ) const {
        for( std::size_t i = range.begin(); i != range.end(); ++i )
            (*this)( i );
    }
#endif // CGAL_LINKED_WITH_TBB
    void operator()( const std::size_t& i ) const {

      Sphere sp(_pts[i], _sq_rd);

      internal::Scale_space::Inc inc(_nn[i]);
      _tree.search(boost::make_function_output_iterator(inc),sp);
    }
}; // class ComputeNN


// Advance a point to a coarser scale.
template < class Gt, class FS, class wA, class Ct >
class
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
AdvanceSS {
public:
    typedef std::vector< Point >            Pointset;
    typedef std::vector< unsigned int >	    CountVec;

private:
    const Search_tree&  _tree;
    const CountVec&     _nn;
    Pointset&           _pts;
    
public:
    AdvanceSS(const Search_tree& tree, const CountVec& nn, Pointset& points)
      : _tree(tree), _nn(nn),_pts(points) {}
    
#ifdef CGAL_LINKED_WITH_TBB
    void operator()( const tbb::blocked_range< std::size_t >& range ) const {
        for( std::size_t i = range.begin(); i != range.end(); ++i )
            (*this)( i );
    }
#endif // CGAL_LINKED_WITH_TBB
    void operator()( const std::size_t& i ) const {
        // If the neighborhood is too small, the vertex is not moved.
        if( _nn[i] < 4 )
            return;

        Static_search search(_tree, _pts[i], _nn[i]);

	Point barycenter (0., 0., 0.);
	FT weight_sum = 0.;
        unsigned int column = 0;
	// Compute total weight
        for( typename Static_search::iterator nit = search.begin();
             nit != search.end() && column < _nn[i];
             ++nit, ++column )
	  weight_sum += (1.0 / _nn[boost::get<1>(nit->first)]);

	column = 0;
	// Compute barycenter
        for( typename Static_search::iterator nit = search.begin();
             nit != search.end() && column < _nn[i];
             ++nit, ++column )
	  {
	    Vector v (CGAL::ORIGIN, boost::get<0>(nit->first));
	    barycenter = barycenter + ((1.0 / _nn[boost::get<1>(nit->first)]) / weight_sum) * v;
	  }
	
	CGAL::cpp11::array<FT, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};
	column = 0;
	// Compute covariance matrix of Weighted PCA
        for( typename Static_search::iterator nit = search.begin();
             nit != search.end() && column < _nn[i];
             ++nit, ++column )
	  {
	    Vector v (barycenter, boost::get<0>(nit->first));
	    FT w = (1.0 / _nn[boost::get<1>(nit->first)]);
	    v = w*v;
	    covariance[0] += w * v.x () * v.x ();
	    covariance[1] += w * v.x () * v.y ();
	    covariance[2] += w * v.x () * v.z ();
	    covariance[3] += w * v.y () * v.y ();
	    covariance[4] += w * v.y () * v.z ();
	    covariance[5] += w * v.z () * v.z ();
	  }

        // Compute the weighted least-squares planar approximation of the point set.
	CGAL::cpp11::array<FT, 9> eigenvectors = {{ 0., 0., 0.,
						    0., 0., 0.,
						    0., 0., 0. }};
	CGAL::cpp11::array<FT, 3> eigenvalues = {{ 0., 0., 0. }};
	wA::diagonalize_selfadjoint_covariance_matrix
	  (covariance, eigenvalues, eigenvectors);

        // The vertex is moved by projecting it onto the plane
        // through the barycenter and orthogonal to the Eigen vector with smallest Eigen value.
	Vector norm (eigenvectors[0], eigenvectors[1], eigenvectors[2]);
	Vector b2p (barycenter, _pts[i]);

	_pts[i] = barycenter + b2p - ((norm * b2p) * norm);
    }
}; // class AdvanceSS


template < class Gt, class FS, class wA, class Ct >
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
Scale_space_surface_reconstruction_3( unsigned int neighbors, unsigned int samples)
  : _generator (0),_mean_neighbors(neighbors), _samples(samples), _squared_radius(-1), _shape(0) {}

template < class Gt, class FS, class wA, class Ct >
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
Scale_space_surface_reconstruction_3( FT sq_radius )
: _mean_neighbors(0), _samples(0), _squared_radius(sq_radius),_shape(0) {
    CGAL_precondition( sq_radius >= 0 );
}

template < class Gt, class FS, class wA, class Ct >
inline bool
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
is_handled( Cell_handle c, unsigned int li ) const {
    switch( li ) {
    case 0: return ( c->info()&1 ) != 0;
    case 1: return ( c->info()&2 ) != 0;
    case 2: return ( c->info()&4 ) != 0;
    case 3: return ( c->info()&8 ) != 0;
    }
    return false;
}

template < class Gt, class FS, class wA, class Ct >
inline void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
mark_handled( Cell_handle c, unsigned int li ) {
    switch( li ) {
    case 0: c->info() |= 1; return;
    case 1: c->info() |= 2; return;
    case 2: c->info() |= 4; return;
    case 3: c->info() |= 8; return;
    }
}

template < class Gt, class FS, class wA, class Ct >
inline void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
mark_opposite_handled( Facet f )  {

  Classification_type cl = _shape->classify (f);

  // If cell is singular, simply mark mirror facet as handled
  if (cl == Shape::SINGULAR)
    {
      Facet mirror = _shape->mirror_facet (f);
      mark_handled (mirror);
    }
  // If cell is regular, get corresponding bubble and mark
  // facets of the other layer of the bubble as handled
  else if (cl == Shape::REGULAR)
    {
      Facet fac = (_shape->classify (f.first) == Shape::EXTERIOR)
	? f
	: _shape->mirror_facet (f);

      typename std::map<Facet, std::size_t>::iterator
	search = _map_f2b.find (fac);

      if (search == _map_f2b.end ())
	return;
      
      unsigned int layer = (_bubbles[search->second][0].find (fac) == _bubbles[search->second][0].end ())
	? 0 : 1;

      typename std::set<Facet>::iterator it = _bubbles[search->second][layer].begin ();

      // If bubble has already been handled, no need to do it again
      if (is_handled (*it))
	return;
      
      for (;it != _bubbles[search->second][layer].end (); ++ it)
	{
	  _garbage.push_back (ordered_facet_indices (*it));
	  mark_handled (*it);
	}
      
    }


}


  template < class Gt, class FS, class wA, class Ct >
inline typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::Triple
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
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

template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
collect_shell( Cell_handle c, unsigned int li, bool separate_shells, bool force_manifold ) {
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
	  if (separate_shells || _shells.size () == 0)
	    {
	      _shells.push_back( --_surface.end() );
	      _index ++;
	    }
            processed = true;
        }

	// Save in which shell the facet is stored
	_map_f2s[f] = _index-1;

	// If the surface is required to be manifold,
	// the opposite layer should be ignored
	if (force_manifold)
	  mark_opposite_handled (f);
		
        // Pivot over each of the facet's edges and continue the surface at the next regular or singular facet.
        for( int i = 0; i < 4; ++i ) {
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

template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
collect_facets (bool separate_shells, bool force_manifold) {


  // We check each of the facets: if it is not handled and either regular or singular,
  // we start collecting the next surface from here.
  for( Finite_facets_iterator fit = _shape->finite_facets_begin(); fit != _shape->finite_facets_end(); ++fit ) {
    switch( _shape->classify( *fit ) ) {
    case Shape::REGULAR:
      // Build a surface from the outer cell.
      if( _shape->classify(fit->first) == Shape::EXTERIOR )
	collect_shell( *fit, separate_shells, force_manifold );
      else
	collect_shell( _shape->mirror_facet( *fit ), separate_shells, force_manifold );
      break;
    case Shape::SINGULAR:
      // Build a surface from both incident cells.
      collect_shell( *fit, separate_shells, force_manifold );
      collect_shell( _shape->mirror_facet( *fit ), separate_shells, force_manifold );
      break;
    default:
      break;
    }
  }

}

template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
collect_facets_quick( ) {
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
            default:
                break;
        }
    }
}

template < class Gt, class FS, class wA, class Ct >
const typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::Shape&
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
shape() const {
	if( !has_shape() )
        _shape = Shape_construction_3().construct( _shape, _squared_radius );
    return *_shape;
}

template < class Gt, class FS, class wA, class Ct >
typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::FT
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
estimate_neighborhood_squared_radius( unsigned int neighbors, unsigned int samples ) {
    typename Gt::Compute_squared_distance_3 squared_distance = Gt().compute_squared_distance_3_object();

    unsigned int handled = 0;
    unsigned int checked = 0;
    FT radius = 0;
        
    if( !_tree.is_built() )
        _tree.build();

    for(typename Search_tree::const_iterator it = _tree.begin(); it != _tree.end(); ++it ) {
        unsigned int left = (unsigned int)( _tree.size() - handled );
        if( samples >= left || _generator.get_double() < double(samples - checked) / left ) {
            // The neighborhood should contain the point itself as well.
          Static_search search( _tree, boost::get<0>(*it), neighbors+1 );
          radius += sqrt( to_double( squared_distance( boost::get<0>(*it), boost::get<0>(( search.end()-1 )->first ) ) ) );
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
template < class Gt, class FS, class wA, class Ct >
template < class InputIterator >
typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::FT
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
estimate_neighborhood_squared_radius( InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples,
                              typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* ) {
    clear();
	insert( begin, end );
	return estimate_neighborhood_squared_radius( neighbors, samples );
}
#endif // DOXYGEN_RUNNING

template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
increase_scale( unsigned int iterations ) {
    typedef std::vector< unsigned int >		CountVec;

    // This method must be called after filling the point collection.
    if( iterations == 0 || _tree.empty() ) return;

    if( !has_neighborhood_squared_radius() )
        estimate_neighborhood_squared_radius();

    for( unsigned int iter = 0; iter < iterations; ++iter ) {

        _tree.clear();
	_tree.reserve (_points.size());

	for (std::size_t i = 0; i < _points.size (); ++ i)
	  _tree.insert (boost::make_tuple(_points[i], i));

        if( !_tree.is_built() )
            _tree.build();

        // Collect the number of neighbors of each point.
        // This can be done concurrently.
        CountVec neighbors( _tree.size(), 0 );
        try_parallel( ComputeNN( _points, _tree, _squared_radius, neighbors ), 0, _tree.size() );
       
        // Compute the transformed point locations.
        // This can be done concurrently.
        try_parallel( AdvanceSS( _tree, neighbors, _points ), 0, _tree.size() );

    }

}

template < class Gt, class FS, class wA, class Ct >
template< class F >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
try_parallel( const F& func, std::size_t begin, std::size_t end, Sequential_tag ) const {
    for( std::size_t i = begin; i < end; ++i ) func( i );
}
    
template < class Gt, class FS, class wA, class Ct >
template< class F >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
try_parallel( const F& func, std::size_t begin, std::size_t end, Parallel_tag ) const {
#ifdef CGAL_LINKED_WITH_TBB
    tbb::parallel_for( tbb::blocked_range< std::size_t >( begin, end ), func );
#else // CGAL_LINKED_WITH_TBB
    try_parallel( func, begin, end, Sequential_tag() );
#endif // CGAL_LINKED_WITH_TBB
}

template < class Gt, class FS, class wA, class Ct >
template < class InputIterator >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
construct_shape( InputIterator begin, InputIterator end,
                 typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* ) {
    deinit_shape();
    if( !has_neighborhood_squared_radius() )
        estimate_neighborhood_squared_radius();
    _shape = Shape_construction_3().construct( begin, end, _squared_radius );
}

template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
collect_surface (bool separate_shells, bool force_manifold, FT border_angle) {

    clear_surface();
    if( !has_shape() )
        construct_shape();

    if( !has_neighborhood_squared_radius() )
      estimate_neighborhood_squared_radius();

    // If shells are not separated and no manifold constraint is given,
    // then the quick collection of facets can be applied
    if (!separate_shells && !force_manifold)
      {
	collect_facets_quick ();
	_shells.push_back (_surface.begin ());
      }
    else
      {
	// Init shell index
	_index = 0;

	// Collect all surface meshes from the alpha-shape in a fashion similar to ball-pivoting.
	// Reset the facet handled markers.
	for( All_cells_iterator cit = _shape->all_cells_begin(); cit != _shape->all_cells_end(); ++cit )
	  cit->info() = 0;

	if (force_manifold)
	  {
	    // If manifold surface is wanted, small volumes (bubbles) are first detected in
	    // order to know which layer to ignore when collecting facets
	    detect_bubbles(border_angle);
	  }

	collect_facets (separate_shells, force_manifold);

	if (force_manifold)
	  {
	    // Even when taking into account facets, some nonmanifold features might remain
	    fix_nonmanifold_edges();
	    fix_nonmanifold_vertices();
	  }
      }
}

template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
find_two_other_vertices(const Facet& f, Vertex_handle v,
			Vertex_handle& v1, Vertex_handle& v2) {
  Vertex_handle vother = f.first->vertex (f.second);
  bool v1found = false;
  
  for (unsigned int i = 0; i < 4; ++ i)
    {
      Vertex_handle vi = f.first->vertex (i);
      if (vi != v && vi != vother)
	{
	  if (v1found)
	    {
	      v2 = vi;
	      return;
	    }
	  else
	    {
	      v1 = vi;
	      v1found = true;
	    }
	}
    }
}




template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
detect_bubbles(FT border_angle) {

  std::set<Cell_handle> done;

  unsigned int nb_facets_removed = 0;
  
  unsigned int nb_skipped = 0;
  for (Cell_iterator cit = _shape->cells_begin (); cit != _shape->cells_end (); ++ cit)
    {
      if (_shape->is_infinite (cit))
	continue;
      if (_shape->classify (cit) != Shape::INTERIOR)
	continue;
      if (done.find (cit) != done.end ())
	continue;

      std::set<VEdge> borders;
      std::vector<Cell_handle> cells;
      std::stack<Cell_handle> todo;
      todo.push (cit);

      // Get all cells of volume and all borders
      while (!(todo.empty ()))
	{
	  Cell_handle c = todo.top ();
	  todo.pop ();

	  if (!(done.insert (c).second))
	    continue;

	  cells.push_back (c);

	  for (unsigned int i = 0; i < 4; ++ i)
	    {
	      if (_shape->classify (c->neighbor (i)) == Shape::INTERIOR)
		todo.push (c->neighbor (i));
	      else
		{
		  // Test if edge is proper border
		  for (unsigned int j = 0; j < 3; ++ j)
		    {
		      unsigned int i0 = (i + j + 1)%4;
		      unsigned int i1 = (i + (j+1)%3 + 1)%4;
                      CGAL_assertion (i0 != i && i1 != i);
		      Edge edge (c, i0, i1);

		      if (_shape->classify (edge) != Shape::REGULAR)
			continue;

		      VEdge vedge = (c->vertex (i0) < c->vertex (i1))
			? std::make_pair (c->vertex (i0), c->vertex (i1))
			: std::make_pair (c->vertex (i1), c->vertex (i0));

		      if (borders.find (vedge) != borders.end ())
			continue;

		      Facet_circulator start = _shape->incident_facets (edge);
		      Facet_circulator circ = start;
		      unsigned int cnt = 0;
		      do
			{
			  if (_shape->classify (*circ) == Shape::SINGULAR
			      || _shape->classify (*circ) == Shape::REGULAR)
			    ++ cnt;
			  ++ circ;
			}
		      while (circ != start);

		      // If edge is non-manifold, use as border
		      if (cnt > 2)
			{
			  borders.insert (vedge);
			  continue;
			}

		      // Else, if facets in cell are regular and angle is
		      // under _border_angle limit, use as border
		      Facet f0 (c, i);
		      Facet f1 (c, (i + (j+2)%3 + 1)%4);

		      if (_shape->classify (f0) != Shape::REGULAR
			  || _shape->classify (f1) != Shape::REGULAR)
			continue;
		      
		      double angle = Gt().compute_approximate_dihedral_angle_3_object()(vedge.first->point (),
                                                                       vedge.second->point (),
                                                                       c->vertex (i)->point (),
                                                                       c->vertex ((i + (j+2)%3 + 1)%4)->point ());

		      if (-border_angle < angle && angle < border_angle)
			{
			  borders.insert (vedge);
			}
		    }

		}
	    }
	}

      int layer = -1;

      // Try to generate bubble from the volume found
      _bubbles.push_back (Bubble());
      std::set<Facet> done;
      for (unsigned int c = 0; c < cells.size (); ++ c)
	{
	  for (unsigned int ii = 0; ii < 4; ++ ii)
	    {
	      Facet start = _shape->mirror_facet (Facet (cells[c], ii));

	      if (_shape->classify (start) != Shape::REGULAR)
		continue;

	      if (done.find (start) != done.end ())
		continue;

	      ++ layer;
	      
	      std::stack<Facet> stack;
	      stack.push (start);

	      Facet f;
	      Cell_handle n, p;
	      int ni, pi;
	      Vertex_handle a;
	      Classification_type cl;

	      // A bubble is well formed is the border contains one loop that
	      // separates two layers.
	      // If the number of layers is different than 2, the volume is completely ignored.
	      while( !stack.empty() )
		{
		  f = stack.top();
		  stack.pop();

		  if (!(done.insert (f).second))
		    continue;

		  if (_shape->classify (f.first) == Shape::EXTERIOR)
		    {
		      if (layer < 2)
			{
			  _bubbles.back ()[layer].insert (f);
			  _map_f2b[f] = _bubbles.size () - 1;
			}
		      else
			{
			  nb_facets_removed ++;
			  mark_handled (f);
			  _garbage.push_back (ordered_facet_indices (f));
			}
		    }
		  else
		    {
		      if (layer < 2)
			{
			  _bubbles.back ()[layer].insert (_shape->mirror_facet (f));
			  _map_f2b[_shape->mirror_facet (f)] = _bubbles.size () - 1;
			}
		      else
			{
			  nb_facets_removed ++;
			  mark_handled (_shape->mirror_facet (f));
			  _garbage.push_back (ordered_facet_indices (_shape->mirror_facet (f)));
			}
 		    }


		  for( int i = 0; i < 4; ++i )
		    {
		      // Skip the current facet.
		      if( i == f.second)
			continue;

		      n = f.first;
		      ni = i;
		      a = f.first->vertex( f.second );
		      cl = _shape->classify( Facet(n, ni) );

		      int n0 = -1, n1 = -1;
		      bool n0found = false;
		      for (int j = 0; j < 4; ++ j)
			{
			  if (j != ni && j != f.second)
			    {
			      if (n0found)
				{
				  n1 = j;
				  break;
				}
			      else
				{
				  n0 = j;
				  n0found = true;
				}
			    }
			}

		      VEdge vedge = (n->vertex (n0) < n->vertex (n1))
			? std::make_pair (n->vertex (n0), n->vertex (n1))
			: std::make_pair (n->vertex (n1), n->vertex (n0));

		      // If the edge is a border, propagation stops in this direction.
		      if (borders.find (vedge) != borders.end ())
			continue;
		      
		      while( cl != Shape::REGULAR && cl != Shape::SINGULAR ) {
			p = n;
			n = n->neighbor(ni);
			ni = n->index(a);
			pi = n->index(p);
			a = n->vertex(pi);
			cl = _shape->classify( Facet(n, ni) );
		      }
		      
		      stack.push (Facet (n, ni));

		    }

		}
	    }

	}
            
      // If number of layers is != 2, ignore volume and discard bubble
      if (layer != 1)
	{
	  nb_skipped ++;
	  for (unsigned int i = 0; i < 2; ++ i)
	    for (typename std::set<Facet>::iterator fit = _bubbles.back()[i].begin ();
		   fit != _bubbles.back()[i].end (); ++ fit)
	      {
		mark_handled (*fit);
		_map_f2b.erase (*fit);
		_garbage.push_back (ordered_facet_indices (*fit));
		nb_facets_removed ++;
	      }
	  _bubbles.pop_back ();
	}
      
    }
}
  

template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
fix_nonmanifold_edges() {

  typedef std::map<std::pair<VEdge, unsigned int>, std::set<Triple> > Edge_shell_map_triples;
  typedef typename Edge_shell_map_triples::iterator Edge_shell_map_triples_iterator;

  unsigned int nb_facets_removed = 0;

  unsigned int nb_nm_edges = 0;

  // Store for each pair edge/shell the incident facets
  Edge_shell_map_triples eshell_triples;
  std::map<Triple, Facet> map_t2f;
  
  for (typename Map_facet_to_shell::iterator fit = _map_f2s.begin ();
       fit != _map_f2s.end (); ++ fit)
    {
      Facet f = fit->first;
      Triple t = ordered_facet_indices (f);
      map_t2f[t] = f;
	      
      for (unsigned int k = 0; k < 3; ++ k)
	{
	  Vertex_handle v0 = f.first->vertex ((f.second + k + 1)%4);
	  Vertex_handle v1 = f.first->vertex ((f.second + (k+1)%3 + 1)%4);
	  VEdge vedge = (v0 < v1) ? std::make_pair (v0, v1) : std::make_pair (v1, v0);
		  
	  std::pair<Edge_shell_map_triples_iterator, bool>
	    search = eshell_triples.insert (std::make_pair (std::make_pair (vedge, fit->second),
							    std::set<Triple>()));

	  search.first->second.insert (t);
	}
    }

  for (Edge_shell_map_triples_iterator eit = eshell_triples.begin ();
       eit != eshell_triples.end (); ++ eit)
    {
      // If an edge has more than 2 incident facets for one shell, it is non-manifold
      if (eit->second.size () < 3)
	continue;

      ++ nb_nm_edges;

      Triple_iterator tit = _shells[eit->first.second];
      Triple_iterator end = (eit->first.second == _shells.size () - 1)
	? _surface.end () : _shells[eit->first.second + 1];

      // Remove facets until the edge is manifold in this shell
      while (tit != end && eit->second.size () > 2)
	{
	  Triple_iterator current = tit ++;

	  typename std::set<Triple>::iterator search = eit->second.find (*current);

	  if (search != eit->second.end ())
	    {
	      if (current == _shells[eit->first.second])
		_shells[eit->first.second] = tit;

	      _garbage.push_back (*current);
	      _map_f2s.erase (map_t2f[*current]);
	      _surface.erase (current);

	      ++ nb_facets_removed;
	      eit->second.erase (search);
	    }

	}
	  
    }
}



  template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
fix_nonmanifold_vertices() {

  typedef ::CGAL::Union_find<Facet> UF;
  typedef typename UF::handle UF_handle;


  typedef std::map<std::pair<Vertex_handle, unsigned int>, std::vector<Facet> > Vertex_shell_map_facets;
  typedef typename Vertex_shell_map_facets::iterator Vertex_shell_map_facet_iterator;

  // For faster facet removal, we sort the triples of each shell as a preprocessing
  for (unsigned int i = 0; i < _shells.size (); ++ i)
    {
      Triple_iterator begin = _shells[i];
      Triple_iterator end = (i+1 == _shells.size ()) ? _surface.end () : _shells[i+1];
      
      Tripleset tmp;
      tmp.splice (tmp.end(), _surface, begin, end);
      
      tmp.sort();
      _shells[i] = tmp.begin ();      
      _surface.splice(end, tmp, tmp.begin(), tmp.end());
    }

  unsigned int nb_facets_removed = 0;
  unsigned int nb_nm_vertices = 0;
  // Removing facets to fix non-manifold vertices might make some other vertices
  // become non-manifold, therefore we iterate until no facet needs to be removed.
  do
    {
      nb_nm_vertices = 0;
      nb_facets_removed = 0;

      // Store for each pair vertex/shell the incident facets
      Vertex_shell_map_facets vshell_facets;

      for (typename Map_facet_to_shell::iterator fit = _map_f2s.begin ();
	   fit != _map_f2s.end (); ++ fit)
	{
	  Facet f = fit->first;
	  
	  for (unsigned int k = 0; k < 3; ++ k)
	    {
	      Vertex_handle v = f.first->vertex ((f.second+k+1)%4);

	      std::pair<Vertex_shell_map_facet_iterator, bool>
		search = vshell_facets.insert (std::make_pair (std::make_pair (v, fit->second),
							       std::vector<Facet>()));
	      search.first->second.push_back (f);

	    }
	  
	}

      for (Vertex_shell_map_facet_iterator fit = vshell_facets.begin ();
	   fit != vshell_facets.end (); ++ fit)
	{
	  if (fit->second.size () < 2)
	    continue;

	  Vertex_handle vit = fit->first.first;
	  unsigned int shell = fit->first.second;

	  UF uf;
	  std::map<Facet, UF_handle> map_f2h;
	  
	  for (unsigned int i = 0; i < fit->second.size (); ++ i)
	    map_f2h.insert (std::make_pair (fit->second[i], uf.make_set (fit->second[i])));

	  std::map<Vertex_handle, Facet> map_v2f;
	    
	  for (unsigned int i = 0; i < fit->second.size (); ++ i)
	    {
	      Vertex_handle v1, v2;
	      find_two_other_vertices (fit->second[i], vit, v1, v2);
	      std::pair<typename std::map<Vertex_handle, Facet>::iterator, bool>
		insertion1 = map_v2f.insert (std::make_pair (v1, fit->second[i]));
	      if (!(insertion1.second))
		uf.unify_sets (map_f2h[fit->second[i]], map_f2h[insertion1.first->second]);
	      std::pair<typename std::map<Vertex_handle, Facet>::iterator, bool>
		insertion2 = map_v2f.insert (std::make_pair (v2, fit->second[i]));
	      if (!(insertion2.second))
		uf.unify_sets (map_f2h[fit->second[i]], map_f2h[insertion2.first->second]);
	    }

	  if (uf.number_of_sets () > 1)
	    {
	      ++ nb_nm_vertices;

	      typedef std::map<UF_handle, std::vector<Facet>, internal::Scale_space::operator_less<UF_handle> > Map_uf_sets;
	      Map_uf_sets map_h2f;
	      for (unsigned int i = 0; i < fit->second.size (); ++ i)
		{
		  UF_handle handle = uf.find (map_f2h[fit->second[i]]);

		  std::pair<typename Map_uf_sets::iterator, bool>
		    insertion = map_h2f.insert (std::make_pair (handle, std::vector<Facet>()));

		  insertion.first->second.push_back (fit->second[i]);
		}

	      typename Map_uf_sets::iterator largest = map_h2f.end ();
	      std::size_t nb_largest = 0;
	      for (typename Map_uf_sets::iterator ufit = map_h2f.begin (); ufit != map_h2f.end (); ++ ufit)
		{
		  std::size_t size = ufit->second.size ();
		  if (size > nb_largest)
		    {
		      nb_largest = size;
		      largest = ufit;
		    }
		}

	      std::vector<Triple> triples;

	      for (typename Map_uf_sets::iterator ufit = map_h2f.begin (); ufit != map_h2f.end (); ++ ufit)
		{
		  if (ufit == largest)
		    continue;
		  for (unsigned int i = 0; i < ufit->second.size (); ++ i)
		    {
		      _map_f2s.erase (ufit->second[i]);
		      triples.push_back (ordered_facet_indices (ufit->second[i]));
		    }
		}
	      std::sort (triples.begin (), triples.end ());

	      Triple_iterator tit = _shells[shell];
	      Triple_iterator end = (shell == _shells.size () - 1)
		? _surface.end () : _shells[shell + 1];

	      unsigned int tindex = 0;
	      
	      while (tit != end && tindex < triples.size ())
		{
		  Triple_iterator current = tit ++;

		  if (*current == triples[tindex])
		    {
		      if (current == _shells[shell])
			_shells[shell] = tit;

		      _garbage.push_back (*current);
		      _surface.erase (current);

		      ++ nb_facets_removed;
		      ++ tindex;
		    }

		}
	    }

	}

    }
  while (nb_nm_vertices != 0);

}



  template < class Gt, class FS, class wA, class Ct >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
  reconstruct_surface( unsigned int iterations, bool separate_shells,
		       bool force_manifold, FT border_angle) {

    // Smooth the scale space.
    increase_scale( iterations );

    // Mesh the perturbed points.
    collect_surface (separate_shells, force_manifold, border_angle);

}

/// \cond internal_doc
template < class Gt, class FS, class wA, class Ct >
template < class InputIterator >
void
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::
reconstruct_surface( InputIterator begin, InputIterator end, unsigned int iterations,
		     bool separate_shells, bool force_manifold, FT border_angle,
                     typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type*) {
  
    // Compute the radius for which the mean ball would contain the required number of neighbors.
      clear();
    insert( begin, end );

    reconstruct_surface (iterations, separate_shells, force_manifold, border_angle);
}
/// \endcond

template < class Gt, class FS, class wA, class Ct >
typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::Triple_const_iterator
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::shell_begin( std::size_t shell ) const {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    return _shells[ shell ];
}

template < class Gt, class FS, class wA, class Ct >
typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::Triple_iterator
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::shell_begin( std::size_t shell ) {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    return _shells[ shell ];
}

template < class Gt, class FS, class wA, class Ct >
typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::Triple_const_iterator
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::shell_end( std::size_t shell ) const {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    if( shell == _shells.size()-1 )
        return _surface.end();
    return _shells[ shell+1 ];
}

template < class Gt, class FS, class wA, class Ct >
typename Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::Triple_iterator
Scale_space_surface_reconstruction_3<Gt,FS,wA,Ct>::shell_end( std::size_t shell ) {
    CGAL_assertion( shell >= 0 && shell < _shells.size() );
    if( shell == _shells.size()-1 )
        return _surface.end();
    return _shells[ shell+1 ];
}

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_IMPL_H
