//A method to construct a surface.
//Copyright (C) 2013  INRIA - Sophia Antipolis
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


#ifndef CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTER_H
#define CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTER_H

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <omp.h>

#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/internal/Shape_type.h>

#include <CGAL/utility.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Random.h>

#include <CGAL/internal/check3264.h>

#include <Eigen/Dense>


namespace CGAL {

/// Compute a smoothed surface mesh from a collection of points.
/** An appropriate neighborhood size is estimated, followed by a
 *  number of smoothing iterations. Finally, the surface of the
 *  smoothed point set is computed.
 *
 *  The order of the point set remains the same, meaning that
 *  the smoothed surface can be transposed back to its unsmoothed
 *  version by overwriting the smoothed point collection with its
 *  unsmoothed version.
 *  \tparam Kernel the geometric traits class. This class
 *  specifies, amongst other things, the number types and
 *  predicates used.
 *  \tparam Fixed_shape indicates whether shape of the object
 *  should be constructed for a fixed neighborhood size.
 *
 *  Generally, constructing for a fixed neighborhood size is more
 *  efficient. This is not the case if the surface should be
 *  constructed for different neighborhood sizes without changing
 *  the point set or recomputing the scale-space.
 *  \tparam Shells indicates whether to collect the surface per shell.
 *  
 *  A shell is a connected component of the surface where connected
 *  facets are locally oriented towards the same side of the surface.
 *  \sa Mean_neighborhood_ball.
 *  \sa Scale_space_transform.
 *  \sa Surface_mesher.
 */
template < class Kernel, class Fixed_shape = CGAL::Tag_true, class Shells = CGAL::Tag_true >
class Scale_space_surface_reconstructer_3 {
    // Searching for neighbors.
	typedef typename CGAL::Search_traits_3< Kernel >    Search_traits;
	typedef typename CGAL::Orthogonal_k_neighbor_search< Search_traits >
                                                        Static_search;
	typedef typename CGAL::Orthogonal_incremental_neighbor_search< Search_traits >
                                                        Dynamic_search;
	typedef typename Dynamic_search::Tree               Search_tree;

	typedef CGAL::Random                                Random;

    // Constructing the surface.
    typedef internal::Shape_type< Kernel, Fixed_shape > Shape_generator;

 	typedef typename Shape_generator::Shape             Shape;
	typedef typename Shape::Vertex_handle				Vertex_handle;
	typedef typename Shape::Cell_handle					Cell_handle;
	typedef typename Shape::Facet						Facet;

	typedef typename Shape::All_cells_iterator          All_cells_iterator;
	typedef typename Shape::Finite_facets_iterator      Finite_facets_iterator;
	typedef typename Shape::Classification_type         Classification_type;

public:
    typedef Shells                                      Collect_per_shell;      ///< Whether to collect the surface per shell.

	typedef typename Kernel::FT                         Scalar;                 ///< The number type.

	typedef typename Kernel::Point_3                    Point;                  ///< The point type.
	typedef typename Kernel::Triangle_3                 Triangle;               ///< The triangle type.

    typedef Triple< unsigned int, unsigned int, unsigned int >
                                                        Triple;                 ///< A triangle of the surface.
    typedef std::list< Triple >                         Tripleset;              ///< A collection of triples.

    typedef typename Search_tree::iterator              Point_iterator;         ///< An iterator over the points.
    typedef typename Search_tree::const_iterator        Point_const_iterator;   ///< A constant iterator over the points.

private:
	typedef typename std::vector<Point>                 Pointset;               ///< A collection of points.

private:
    Search_tree _tree;              // To quickly search for nearest neighbors.

	Random      _generator;         // For sampling random points.

    unsigned int _mean_neighbors;   // The number of nearest neighbors in the mean neighborhood.
    unsigned int _samples;          // The number of sample points for estimating the mean neighborhood.

    Scalar      _squared_radius;    // The squared mean neighborhood radius.

    // The shape must be a pointer, because the alpha of
    // a Fixed_alpha_shape_3 can only be set at
    // construction and its assignment operator is private.
	Shape* _shape;

    // The surface. If the surface is collected per shell,
    // consecutive triples belong to the same shell and
    // different shells are separated by a (0,0,0) triple.
    Tripleset       _surface;

private:
    void clear_tree() { _tree.clear(); }

public:
    /// Default constructor.
    /** \param sq_radius the squared radius of the
     *  neighborhood size. If this value is negative when
     *  the point set is smoothed or when the surface is computed,
     *  the neighborhood size will be computed automatically.
     */
	Scale_space_surface_reconstructer_3(unsigned int neighbors = 30, unsigned int samples = 200, Scalar sq_radius = -1 ): _mean_neighbors(neighbors), _samples(samples), _squared_radius( sq_radius ), _shape(0) {}
	~Scale_space_surface_reconstructer_3() { deinit_shape(); }

private:
    void deinit_shape() { if( _shape != 0 ) { delete _shape; _shape = 0; } }

public:
    Point_iterator scale_space_begin() { return _tree.begin(); }
    Point_iterator scale_space_end() { return _tree.end(); }
    
    Point_const_iterator scale_space_begin() const { return _tree.begin(); }
    Point_const_iterator scale_space_end() const { return _tree.begin(); }

    /// Check whether the spatial structure has been constructed.
    /** \return true if the structure exists and false otherwise.
     *  \sa construct_shape(InputIterator start, InputIterator end).
     */
	bool has_surface() const { return _shape != 0; }

	void clear_surface() {
        if( has_surface() ) {
            _shape->clear();
        }
    }
    
    /// Insert a collection of points.
    /** \tparam InputIterator an iterator over a collection of points.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa compute_surface(InputIterator start, InputIterator end).
     */
	template < class InputIterator >
	void insert_points( InputIterator start, InputIterator end ) {
		_tree.insert( start, end );
	}

    void clear() {
		clear_tree();
        clear_surface();
    }
    
    /// Estimate the mean neighborhood size based on a number of sample points.
    /** A neighborhood size is expressed as the radius of the smallest
     *  ball centered on a point such that the ball contains at least
     *  a specified number of points.
     *
     *  The mean neighborhood size is then the mean of these radii,
     *  taken over a number of point samples.
     * \return the estimated mean neighborhood size.
     *  \sa handlePoint(const Point& p).
     *  \sa operator()(InputIterator start, InputIterator end).
     */
	Scalar estimate_mean_neighborhood( unsigned int neighbors = 30, unsigned int samples = 200 ) {
        Kernel::Compute_squared_distance_3 squared_distance = Kernel().compute_squared_distance_3_object();

        _mean_neighbors = neighbors;
        _samples = samples;

		unsigned int handled = 0;
		unsigned int checked = 0;
		Scalar radius = 0;
        
		if( !_tree.is_built() )
            _tree.build();

		for (typename Search_tree::const_iterator it = _tree.begin(); it != _tree.end(); ++it) {
			if (samples >= (_tree.size() - handled) || _generator.get_double() < double(samples - checked) / (_tree.size() - handled)) {
                // The neighborhood should contain the point itself as well.
		        Static_search search( _tree, *it, _mean_neighbors+1 );
		        radius += CGAL::sqrt( CGAL::to_double( squared_distance( *it, ( search.end()-1 )->first ) ) );
				++checked;
			}
			++handled;
		}
		radius /= double(checked);

        set_mean_neighborhood( radius );
		return radius;
	}

    
	template < class InputIterator >
	Scalar estimate_mean_neighborhood(InputIterator start, InputIterator end, unsigned int neighbors = 30, unsigned int samples = 200) {
		insert_points(start, end);
		return estimate_mean_neighborhood(neighbors, samples);
	}

    // -1 if not yet set.
    Scalar get_squared_mean_neighborhood() const { return _squared_radius; }
    void set_mean_neighborhood( const Scalar& radius ) {
        _squared_radius = radius * radius;
		if( has_surface() )
            _shape = Shape_generator().construct( _shape, _squared_radius );
    }
    bool has_squared_mean_neighborhood() const {
        return sign( _squared_radius ) == POSITIVE;
    }

    unsigned int get_mean_neighbors() const { return _mean_neighbors; }
    unsigned int get_number_samples() const { return _samples; }
    void set_mean_neighbors(unsigned int neighbors) { _mean_neighbors = neighbors;}
    void set_number_samples(unsigned int samples) { _samples = samples; }


    /// Compute a number of iterations of scale-space transforming.
    /** If earlier iterations have been computed, calling smooth_scale_space()
     *  will add more iterations.
     *
     *  If the mean neighborhood is negative, it will be computed first.
     *  \param iterations the number of iterations to perform.
     */
	void smooth_scale_space(unsigned int iterations = 1) {
	
		typedef std::vector<unsigned int>					CountVec;
		typedef typename std::map<Point, size_t>			PIMap;

		typedef Eigen::Matrix<double, 3, Eigen::Dynamic>	Matrix3D;
		typedef Eigen::Array<double, 1, Eigen::Dynamic>		Array1D;
		typedef Eigen::Matrix3d								Matrix3;
		typedef Eigen::Vector3d								Vector3;
		typedef Eigen::SelfAdjointEigenSolver<Matrix3>		EigenSolver;

		typedef _ENV::s_ptr_type							p_size_t;
		
		// This method must be called after filling the point collection.
		CGAL_assertion(!_tree.empty());

        if( !has_squared_mean_neighborhood() )
            estimate_mean_neighborhood( _mean_neighbors, _samples );

        Pointset points;
		points.assign( _tree.begin(), _tree.end() );
        _tree.clear();

		// Construct a search tree of the points.
        // Note that the tree has to be local for openMP.
		for (unsigned int iter = 0; iter < iterations; ++iter) {
            Search_tree tree( points.begin(), points.end() );
		    if(!tree.is_built())
			    tree.build();

		    // Collect the number of neighbors of each point.
		    // This can be done parallel.
		    CountVec neighbors(tree.size(), 0);
		    Kernel::Compare_squared_distance_3 compare;
		    p_size_t count = tree.size(); // openMP can only use signed variables.
            const Scalar squared_radius = _squared_radius; // openMP can only use local variables.
#pragma omp parallel for shared(count,tree,points,squared_radius,neighbors) firstprivate(compare)
		    for (p_size_t i = 0; i < count; ++i) {
			    // Iterate over the neighbors until the first one is found that is too far.
			    Dynamic_search search(tree, points[i]);
			    for (Dynamic_search::iterator nit = search.begin(); nit != search.end() && compare(points[i], nit->first, squared_radius) != CGAL::LARGER; ++nit)
				    ++neighbors[i];
		    }
		
		    // Construct a mapping from each point to its index.
		    PIMap indices;
            p_size_t index = 0;
		    for( Pointset::const_iterator pit = points.begin(); pit != points.end(); ++pit, ++index)
			    indices[ *pit ] = index;

		    // Compute the tranformed point locations.
		    // This can be done parallel.
#pragma omp parallel for shared(count,neighbors,tree,squared_radius) firstprivate(compare)
		    for (p_size_t i = 0; i < count; ++i) {
			    // If the neighborhood is too small, the vertex is not moved.
			    if (neighbors[i] < 4)
				    continue;

			    // Collect the vertices within the ball and their weights.
			    Dynamic_search search(tree, points[i]);
			    Matrix3D pts(3, neighbors[i]);
			    Array1D wts(1, neighbors[i]);
			    int column = 0;
			    for (Dynamic_search::iterator nit = search.begin(); nit != search.end() && compare(points[i], nit->first, squared_radius) != CGAL::LARGER; ++nit, ++column) {
				    pts(0, column) = CGAL::to_double(nit->first[0]);
				    pts(1, column) = CGAL::to_double(nit->first[1]);
				    pts(2, column) = CGAL::to_double(nit->first[2]);
				    wts(column) = 1.0 / neighbors[indices[nit->first]];
			    }

			    // Construct the barycenter.
			    Vector3 bary = (pts.array().rowwise() * wts).rowwise().sum() / wts.sum();
			
			    // Replace the points by their difference with the barycenter.
			    pts = (pts.colwise() - bary).array().rowwise() * wts;

			    // Construct the weighted covariance matrix.
			    Matrix3 covariance = Matrix3::Zero();
			    for (column = 0; column < pts.cols(); ++column)
				    covariance += wts.matrix()(column) * pts.col(column) * pts.col(column).transpose();

			    // Construct the Eigen system.
			    EigenSolver solver(covariance);

			    // If the Eigen system does not converge, the vertex is not moved.
			    if (solver.info() != Eigen::Success)
				    continue;

			    // Find the Eigen vector with the smallest Eigen value.
			    std::ptrdiff_t index;
			    solver.eigenvalues().minCoeff(&index);
			    if (solver.eigenvectors().col(index).imag() != Vector3::Zero()) {
				    // This should actually never happen!
				    CGAL_assertion(false);
				    continue;
			    }
			    Vector3 n = solver.eigenvectors().col(index).real();

			    // The vertex is moved by projecting it onto the plane
			    // through the barycenter and orthogonal to the Eigen vector with smallest Eigen value.
			    Vector3 bv = Vector3(CGAL::to_double(points[i][0]), CGAL::to_double(points[i][1]), CGAL::to_double(points[i][2])) - bary;
			    Vector3 per = bary + bv - (n.dot(bv) * n);

			    points[i] = Point(per(0), per(1), per(2));
		    }
        }

        _tree.insert( points.begin(), points.end() );
	}

    
    /// Compute a number of transform iterations on a collection of points.
    /** This method is equivalent to running [insert_points(start, end)](\ref insert_points)
     *  followed by [smooth_scale_space(iterations)](\ref smooth_scale_space).
     *  \tparam InputIterator an iterator over a collection of points.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \param iterations the number of iterations to perform.
     *  \sa operator()(InputIterator start, InputIterator end).
     *  \sa assign(InputIterator start, InputIterator end).
     *  \sa iterate(unsigned int iterations).
     */
	template < class InputIterator >
	void smooth_scale_space(InputIterator start, InputIterator end, unsigned int iterations = 1) {
		insert_points(start, end);
		smooth_scale_space(iterations);
	}

private:
	// Once a facet is added to the surface, it is marked as handled.
	inline bool is_handled( Cell_handle c, unsigned int li ) const {
		switch( li ) {
		case 0: return ( c->info()&1 ) != 0;
		case 1: return ( c->info()&2 ) != 0;
		case 2: return ( c->info()&4 ) != 0;
		case 3: return ( c->info()&8 ) != 0;
		}
		return false;
	}
	inline bool is_handled( const Facet& f ) const { return is_handled( f.first, f.second ); }
	inline void set_handled( Cell_handle c, unsigned int li ) {
		switch( li ) {
		case 0: c->info() |= 1; return;
		case 1: c->info() |= 2; return;
		case 2: c->info() |= 4; return;
		case 3: c->info() |= 8; return;
		}
	}
	inline void set_handled( Facet f ) { set_handled( f.first, f.second ); }

public:

    // make new construct_shape() method for when pts already known...
    void construct_shape() {
		construct_shape( scale_space_begin(), scale_space_end() );
    }

    /*// For if you already have a Delaunay triangulation of the points.
    void construct_shape(Shape_generator::Structure& tr ) {
        deinit_shape();
        if( !has_squared_mean_neighborhood() )
            estimate_mean_neighborhood( _mean_neighbors, _samples );
        _shape = Shape_generator().construct( *tr, r2 );
	}*/

    // If you don't want to smooth the point set.
    /// Construct a new spatial structure on a set of sample points.
    /** \tparam InputIterator an iterator over the point sample.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa is_constructed().
     */
	template < class InputIterator >
	void construct_shape(InputIterator start, InputIterator end) {
		deinit_shape();
        if( !has_squared_mean_neighborhood() )
            estimate_mean_neighborhood( _mean_neighbors, _samples );
        _shape = Shape_generator().construct( start, end, _squared_radius );
	}

private:
    Triple ordered_facet_indices( const Facet& f ) const {
		if( (f.second&1) == 0 )
            return Triple( f.first->vertex( (f.second+2)&3 )->info(),
                           f.first->vertex( (f.second+1)&3 )->info(),
                           f.first->vertex( (f.second+3)&3 )->info() );
		else
            return Triple( f.first->vertex( (f.second+1)&3 )->info(),
                           f.first->vertex( (f.second+2)&3 )->info(),
                           f.first->vertex( (f.second+3)&3 )->info() );
    }

    /// Collect the triangles of one shell of the surface.
    /** A shell is a maximal collection of triangles that are
     *  connected and locally facing towards the same side of the
     *  surface.
     *  \tparam OutputIterator an output iterator for a collection
     *  of triangles. The iterator must point to a Triangle type.
     *  \param c a cell touching a triangle of the shell from the
     *  outside of the object.
     *  \param li the index of the vertex of the cell opposite to
     *  the triangle touched.
     *  \param out an iterator to place to output the triangles.
     *  \sa collect_shell(const Facet& f, OutputIterator out).
     *  \sa collect_surface(OutputIterator out).
     */
	void collect_shell( Cell_handle c, unsigned int li ) {
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
			CGAL_triangulation_assertion( !_shape->is_infinite(f) );
			set_handled(f);

            // Output the facet as a triple of indices.
			_surface.push_back( ordered_facet_indices(f) );
		
            if( Shells::value ) {
			    // Pivot over each of the facet's edges and continue the surface at the next regular or singular facet.
			    for( unsigned int i = 0; i < 4; ++i ) {
				    // Skip the current facet.
				    if( i == f.second || is_handled(f.first, i) )
					    continue;

				    // Rotate around the edge (starting from the shared facet in the current cell) until a regular or singular facet is found.
				    n = f.first;
				    ni = i;
				    a = f.first->vertex(f.second);
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
                processed = true;
            }
		}
        
        // We indicate the end of a shell by the (0,0,0) triple.
        if( Shells::value && processed )
            _surface.push_back( Triple(0,0,0) );
	}

    /// Collect the triangles of one shell of the surface.
    /** A shell is a maximal collection of triangles that are
     *  connected and locally facing towards the same side of the
     *  surface.
     *  \tparam OutputIterator an output iterator for a collection
     *  of triangles. The iterator must point to a Triangle type.
     *  \param f a facet of the shell as seen from the outside.
     *  \param out an iterator to place to output the triangles.
     *  \sa collect_shell(Cell_handle c, unsigned int li, OutputIterator out).
     *  \sa collect_surface(OutputIterator out).
     */
	void collect_shell( const Facet& f ) {
		collect_shell( f.first, f.second );
	}

public:
    /// Compute a surface mesh from a collection of points.
    /** The surface mesh indicates the boundary of a sampled
     *  object. The point sample must be dense enough, such
     *  that any region that fits an empty ball with a given
     *  radius can be assumed to be ouside the object.
     *
     *  The resulting surface need not be manifold and in
     *  many cases where the points sample the surface of
     *  an object, the surface computed will contain both
     *  an 'outward-facing' and an 'inward-facing' surface.
     *
     *  However, this surface cannot have holes and its
     *  triangles are all oriented towards the same side of
     *  the surface. This orientation is expressed using the
     *  'right-hand rule' on the ordered vertices of each
     *  triangle.
     *
     *  The computed triangles will be returned ordered such
     *  that a connected 'shell' contains a set of consecutive
     *  triangles in the output. For this purpose, a 'shell'
     *  is a maximal collection of triangles that are connected
     *  and locally facing towards the same side of the surface.
     *  \tparam Kernel the geometric traits class. This class
     *  specifies, amongst other things, the number types and
     *  predicates used.
     *  \tparam Vb the vertex base used in the internal data
     *  structure that spatially orders the point set.
     *  \tparam Cb the cell base used in the internal data
     *  structure that spatially orders the point set.
     */
    /// Collect the triangles of the complete surface.
    /** These triangles are given ordered per shell.
     *  A shell is a maximal collection of triangles that are
     *  connected and locally facing towards the same side of the
     *  surface.
     *  \tparam OutputIterator an output iterator for a collection
     *  of triangles. The iterator must point to a Triangle type.
     *  \param out an iterator to place to output the triangles.
     *  \sa collect_shell(const Facet& f, OutputIterator out).
     *  \sa collect_shell(Cell_handle c, unsigned int li, OutputIterator out).
     */
	template < class OutputIterator >
	OutputIterator collect_surface( OutputIterator out ) {
        if( !has_squared_mean_neighborhood() )
            estimate_mean_neighborhood( _mean_neighbors, _samples );

		// Collect all surface meshes from the alpha-shape in a fashion similar to ball-pivoting.
		// Reset the facet handled markers.
		for( All_cells_iterator cit = _shape->all_cells_begin(); cit != _shape->all_cells_end(); ++cit )
			cit->info() = 0;

		// We check each of the facets: if it is not handled and either regular or singular,
		// we start collecting the next surface from here.
		Facet m;
		int ns = 0;
		for( Finite_facets_iterator fit = _shape->finite_facets_begin(); fit != _shape->finite_facets_end(); ++fit ) {
			m = _shape->mirror_facet( *fit );
			switch( _shape->classify( *fit ) ) {
			case Shape::REGULAR:
				if( !is_handled(*fit) && !is_handled(m) )
					++ns;
				// Build a surface from the outer cell.
				if( _shape->classify(fit->first) == Shape::EXTERIOR )
					collect_shell( *fit );
				else
					collect_shell( m );
				break;
			case Shape::SINGULAR:
				if( !is_handled( *fit ) )
					++ns;
				// Build a surface from both incident cells.
				collect_shell( *fit );
				if( !is_handled(m) )
					++ns;
				collect_shell( m );
				break;
			}
		}

		return out;
	}

    /// Construct the surface triangles from a set of sample points.
    /** These triangles are given ordered per shell.
     *  A shell is a maximal collection of triangles that are
     *  connected and locally facing towards the same side of the
     *  surface.
     *
     *  This method is equivalent to running [construct_shape(start, end)](\ref construct_shape)
     *  followed by [collect_surface(out)](\ref collect_surface).
     *  \tparam InputIterator an iterator over the point sample.
     *  The iterator must point to a Point type.
     *  \tparam OutputIterator an output iterator for a collection
     *  of triangles. The iterator must point to a Triangle type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \param out an iterator to place to output the triangles.
     *  \sa construct_shape(InputIterator start, InputIterator end).
     *  \sa collect_surface(OutputIterator out).
     */
	template < class InputIterator, class OutputIterator >
	OutputIterator collect_surface( InputIterator start, InputIterator end, OutputIterator out ) {
		construct_shape( start, end );
		return collect_surface( out );
	}

	template < class InputIterator >
	void collect_surface( InputIterator start, InputIterator end ) {
		construct_shape( start, end );
		collect_surface( std::back_inserter(_surface) );
	}

    void collect_surface() {
		construct_shape();
		collect_surface( std::back_inserter(_surface) );
    }

    const Tripleset& surface() const { return _surface; }
    Tripleset& surface() { return _surface; }


	
	/// Compute a smoothed surface mesh from a collection of points.
    /** This is equivalent to calling insert_points(),
     *  iterate(), and compute_surface().
     *
     *  After this computation, the surface can be accessed using
     *  surface().
     *  \tparam InputIterator an iterator over the point collection.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  
     *  \sa surface().
     *  \sa surface() const.
     */
	template < class InputIterator >
	void compute_surface(InputIterator start, InputIterator end, unsigned int neighbors = 30, unsigned int iterations = 4, unsigned int samples = 200) {
		// Compute the radius for which the mean ball would contain the required number of neighbors.
        clear();
		insert_points( start, end );

        _mean_neighbors = neighbors;
        _samples = samples;

		// Smooth the scale space.
		smooth_scale_space( iterations );

		// Mesh the perturbed points.
        collect_surface();
	}
}; // class Scale_space_surface_reconstructer_3


template< typename T1, typename T2, typename T3 >
std::ostream&
operator<<( std::ostream& os, const Triple< T1, T2, T3 >& t ) {
    return os << t.first << " " << t.second << " " << t.third;
}

template< typename T1, typename T2, typename T3 >
std::istream&
operator>>( std::istream& is, Triple< T1, T2, T3 >& t ) {
    return is >> t.first >> t.second >> t.third;
}

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTER_H