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

#ifndef CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTER_H
#define CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTER_H

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <omp.h>

#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/utility.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Random.h>

#include <CGAL/Scale_space_reconstruction_3/internal/check3264.h>
#include <CGAL/Scale_space_reconstruction_3/Shape_of_points_3.h>

#include <Eigen/Dense>


namespace CGAL {

/// Compute a triangular surface mesh from a collection of points.
/** An appropriate neighborhood size is estimated, followed by a
 *  number of smoothing iterations on the point cloud. Finally, the surface of the
 *  smoothed point set is computed and its connectivity is stored
 *  as triples of indices to the point set. These tasks can either
 *  be performed individually, or as a batch process.
 *
 *  The order of the point set remains the same, meaning that
 *  the surface on the smoothed point set can interpolate the unsmoothed
 *  point set by applying the indices of the surface triangles to the
 *  unsmoothed point set. There are no guarantees on the correctness
 *  of the resulting surface. However, depending on the geometry of the point set,
 *  the number of iterations and the neighborhood size, the surface will generally
 *  not self-intersect.
 *
 *  The shape can be constructed either at a fixed predefined scale,
 *  or at a dynamic scale. The first option is faster when constructing
 *  a single shape. It is undefined whether a shape with fixed scale may
 *  have its scale changed, but if so, this will likely require more time
 *  than changing a dynamic scale. In either case, it is possible to change
 *  the point set while maintaining the same scale.
 *
 *  The surface can be given either as an unordered collection of triangles,
 *  or as a collection of shells. A shell is a connected component of the
 *  surface where connected facets are locally oriented towards the same
 *  side of the surface. Shells are separated by a (0,0,0) triple.
 *  
 *  \tparam Kernel the geometric traits class. It should be a model of
 *  DelaunayTriangulationTraits_3.
 *  \tparam Fixed_scale whether the shape is constructed for a fixed scale.
 *  It should be a model of Boolean_tags. The default value is Tag_true.
 *  \tparam Shells whether to collect the surface per shell. It should be
 *  a model of Boolean_tags. The default value is Tag_true.
 */
#ifdef DOXYGEN_RUNNING
template < class Kernel, class Fixed_scale, class Shells >
#else
template < class Kernel, class Fixed_scale = Tag_true, class Shells = Tag_true >
#endif
class Scale_space_surface_reconstructer_3 {
    // Searching for neighbors.
    typedef typename Search_traits_3< Kernel >      Search_traits;
    typedef typename Orthogonal_k_neighbor_search< Search_traits >
                                                    Static_search;
    typedef typename Orthogonal_incremental_neighbor_search< Search_traits >
                                                    Dynamic_search;
    typedef typename Dynamic_search::Tree           Search_tree;

    typedef Random                                  Random;

    // Constructing the surface.
    typedef Shape_of_points_3< Kernel, Fixed_scale >    Shape_of_points_3;

    typedef typename Shape_of_points_3::Shape       Shape;

    typedef typename Shape::Vertex_handle           Vertex_handle;
    typedef typename Shape::Cell_handle             Cell_handle;
    typedef typename Shape::Facet                   Facet;
    
    typedef typename Shape::Vertex_iterator         Vertex_iterator;
    typedef typename Shape::Cell_iterator           Cell_iterator;
    typedef typename Shape::Facet_iterator          Facet_iterator;
    typedef typename Shape::Edge_iterator           Edge_iterator;

    typedef typename Shape::All_cells_iterator      All_cells_iterator;
    typedef typename Shape::Finite_facets_iterator  Finite_facets_iterator;
    typedef typename Shape::Classification_type     Classification_type;

public:
/// \name Types
/// \{
    typedef Shells                                  Collect_per_shell;      ///< Whether to collect the surface per shell.

	typedef typename Kernel::FT                     FT;                     ///< The number field type.

	typedef typename Kernel::Point_3                Point;                  ///< The point type.
	typedef typename Kernel::Triangle_3             Triangle;               ///< The triangle type.

#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                        Point_iterator;         ///< An iterator over the points.
    typedef const unspecified_type                  Const_point_iterator;   ///< A constant iterator over the points.
#else
    typedef typename Search_tree::iterator          Point_iterator;
    typedef typename Search_tree::const_iterator    Const_point_iterator;
#endif // DOXYGEN_RUNNING

    typedef Triple< unsigned int, unsigned int, unsigned int >
                                                    Triple;                 ///< A triple indicating a triangle of the surface.
private:
    typedef std::list< Triple >                     Tripleset;              ///< A collection of triples.

public:
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                        Triple_iterator;        ///< An iterator over the triples.
    typedef const unspecified_type                  Const_triple_iterator;  ///< A constant iterator over the triples.
#else
    typedef Tripleset::iterator                     Triple_iterator;
    typedef Tripleset::const_iterator               Const_triple_iterator;
#endif // DOXYGEN_RUNNING
/// \}

private:
    class Insurface_tester;
    typedef Filter_iterator< Facet_iterator, Insurface_tester >
                                                    Surface_facets_iterator;

private:
    Search_tree     _tree;              // To quickly search for nearest neighbors.

	Random          _generator;         // For sampling random points.

    unsigned int    _mean_neighbors;    // The number of nearest neighbors in the mean neighborhood.
    unsigned int    _samples;           // The number of sample points for estimating the mean neighborhood.

    FT              _squared_radius;    // The squared mean neighborhood radius.

    // The shape must be a pointer, because the alpha of
    // a Fixed_alpha_shape_3 can only be set at
    // construction and its assignment operator is private.
    // We want to be able to set the alpha after constructing
    // the scale-space reconstructer object.
	Shape*          _shape;

    // The surface. If the surface is collected per shell,
    // consecutive triples belong to the same shell and
    // different shells are separated by a (0,0,0) triple.
    Tripleset       _surface;

public:
/// \name Constructors
/// \{
    /// Construct a surface reconstructor.
    /** \param neighbors the number of neighbors a neighborhood should contain on average.
     *  \param samples the number of points sampled to estimate the mean neighborhood size.
     *  \param sq_radius the squared radius of the
     *  neighborhood size. If this value is negative when
     *  the point set is smoothed or when the surface is computed,
     *  the neighborhood size will be computed automatically.
     */
	Scale_space_surface_reconstructer_3( unsigned int neighbors = 30, unsigned int samples = 200, FT sq_radius = -1 );
/// \}
	~Scale_space_surface_reconstructer_3() { deinit_shape(); }

private:
    void deinit_shape() { if( _shape != 0 ) { delete _shape; _shape = 0; } }
    
    void clear_tree() { _tree.clear(); }
	void clear_surface() { if( has_shape() ) { _shape->clear(); } }
    
	// Once a facet is added to the surface, it is marked as handled.
	bool is_handled( Cell_handle c, unsigned int li ) const;
	inline bool is_handled( const Facet& f ) const { return is_handled( f.first, f.second ); }
	void set_handled( Cell_handle c, unsigned int li );
	inline void set_handled( Facet f ) { set_handled( f.first, f.second ); }

public:
/// \name Accessors
/// \{
    /// Gives the squared radius of the neighborhood ball.
    /** The neighborhood ball is used by
     *  `smooth_scale_space( unsigned int iterations )` to
     *  compute the scale-space and by
     *  `construct_shape()` to construct the shape of
     *  the scale-space.
     *  \return the squared radius of the mean neighborhood,
     *  or -1 if the mean neighborhood has not yet been set.
     */
    FT get_neighborhood_squared_radius() const { return _squared_radius; }

    /// Gives the mean number of neighbors an estimated neighborhood should contain.
    /** This number is only used if the neighborhood radius
     *  has not been set manually.
     *
     *  If the neighborhood ball radius is estimated, it should
     *  on average contain this many neighbors, not counting the
     *  ball center.
     *  \return the number of neighbors a neighborhood ball centered
     *  on a point should contain on average when the radius is estimated.
     */
    unsigned int get_mean_neighbors() const { return _mean_neighbors; }

    /// Gives the number of sample points the neighborhood estimation uses.
    /** This number is only used if the neighborhood radius
     *  has not been set manually.
     *
     *  If the number of samples is larger than the point cloud,
     *  every point is used and the optimal neighborhood radius
     *  is computed exactly in stead of estimated.
     *  \return the number of sample points used for neighborhood estimation.
     */
    unsigned int get_number_neighborhood_samples() const { return _samples; }

    /// Get the shape of the scale space.
    /** If the shape does not exist, it is  constructed first.
     *  \return the shape of the scale space.
     */
    const Shape& get_shape() const;
    
    /// Get the number of triangles of the surface.
    unsigned int get_surface_size() const { return (unsigned int)_surface.size(); }
/// \}

/// \name Mutators
/// \{
    /// Reset the scale-space surface reconstructer.
    void clear() {
		clear_tree();
        clear_surface();
    }
    
    /// Sets the radius of the neighborhood ball.
    /** The neighborhood ball is used by
     *  `smooth_scale_space( unsigned int iterations )` to
     *  compute the scale-space and by
     *  `construct_shape()` to construct the shape of
     *  the scale-space.
     *
     *  If the reconstructor has already constructed a
     *  shape for the scale-space, this may cause the
     *  construction of a new shape.
     *  \param radius the new radius of the mean neighborhood.
     */
    void set_neighborhood_radius( const FT& radius ) {
        _squared_radius = radius * radius;
		if( has_shape() )
            Shape_of_points_3().set_scale( _shape, _squared_radius );
    }
    
    /// Sets the squared radius of the neighborhood ball.
    /** The neighborhood ball is used by
     *  `smooth_scale_space( unsigned int iterations )` to
     *  compute the scale-space and by
     *  `construct_shape()` to construct the shape of
     *  the scale-space.
     *
     *  If the reconstructor has already constructed a
     *  shape for the scale-space, this may cause the
     *  construction of a new shape.
     *  \param radius the new squared radius of the mean neighborhood,
     *  or -1 if the mean neighborhood should be estimated automatically.
     */
    void set_neighborhood_squared_radius( const FT& sq_radius ) {
        _squared_radius = sq_radius;
		if( has_shape() )
            Shape_of_points_3().set_scale( _shape, _squared_radius );
    }
    
    /// Sets the mean number of neighbors an estimated neighborhood should contain.
    /** This number is only used if the neighborhood radius
     *  has not been set manually.
     *
     *  If the neighborhood ball radius is estimated, it should
     *  on average contain this many neighbors, not counting the
     *  ball center.
     */
    void set_mean_neighbors( unsigned int neighbors ) { _mean_neighbors = neighbors; }
    
    /// Sets the number of sample points the neighborhood estimation uses.
    /** This number is only used if the neighborhood radius
     *  has not been set manually.
     *
     *  If the number of samples is larger than the point cloud,
     *  every point is used and the optimal neighborhood radius
     *  is computed exactly in stead of estimated.
     */
    void set_number_neighborhood_samples( unsigned int samples ) { _samples = samples; }

    /// Insert a collection of points.
    /** Note that inserting the points does not automatically construct
     *  the surface.
     *
     *  In order to construct the surface, either run
     *  `construct_surface(unsigned int iterations)` or both
     *  `smooth_scale_space( unsigned int iterations )` and
     *  `collect_surface()`.
     *  \tparam InputIterator an iterator over a collection of points.
     *  The iterator must point to a Point.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa compute_surface( InputIterator start, InputIterator end ).
     */
	template < class InputIterator >
	void insert_points( InputIterator start, InputIterator end ) {
		_tree.insert( start, end );
	}
/// \}

/// \name Query
/// \{
    /// Check whether the neighborhood ball radius has been set.
    /** \return true iff the radius has been set manually or estimated.
     */
    bool has_neighborhood_radius() const {
        return sign( _squared_radius ) == POSITIVE;
    }

    /// Check whether the shape has been constructed.
    /** The shape contains the structure of the point cloud.
     *
     *  Until the shape is constructed, the surface is undefined.
     *  \return true if the shape exists and false otherwise.
     *  \sa construct_shape(InputIterator start, InputIterator end).
     */
	bool has_shape() const { return _shape != 0; }
/// \}

/// \name Iterators
/// \{
    /// Gives an iterator to the first point in the current scale space.
    Const_point_iterator scale_space_begin() const { return _tree.begin(); }
    /// Gives an iterator to the first point in the current scale space.
    Point_iterator scale_space_begin() { return _tree.begin(); }

    /// Gives a past-the-end iterator of the points in the current scale space.
    Const_point_iterator scale_space_end() const { return _tree.begin(); }
    /// Gives a past-the-end iterator of the points in the current scale space.
    Point_iterator scale_space_end() { return _tree.end(); }

    /// Gives an iterator to the first triple in the surface.
    Const_triple_iterator surface_begin() const { return _surface.begin(); }
    /// Gives an iterator to the first triple in the surface.
    Triple_iterator surface_begin() { return _surface.begin(); }
    
    /// Gives a past-the-end iterator of the triples in the surface.
    Const_triple_iterator surface_end() const { return _surface.end(); }
    /// Gives a past-the-end iterator of the triples in the surface.
    Triple_iterator surface_end() { return _surface.end(); }
/// \}

public:

    

    
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
	FT estimate_mean_neighborhood( unsigned int neighbors = 30, unsigned int samples = 200 );

	template < class InputIterator >
	FT estimate_mean_neighborhood(InputIterator start, InputIterator end, unsigned int neighbors = 30, unsigned int samples = 200);



    /// Compute a number of iterations of scale-space transforming.
    /** If earlier iterations have been computed, calling smooth_scale_space()
     *  will add more iterations.
     *
     *  If the mean neighborhood is negative, it will be computed first.
     *  \param iterations the number of iterations to perform.
     */
	void smooth_scale_space(unsigned int iterations = 1);

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

public:

    // make new construct_shape() method for when pts already known...
    void construct_shape() {
		construct_shape( scale_space_begin(), scale_space_end() );
    }

    /*// For if you already have a Delaunay triangulation of the points.
    void construct_shape(Shape_of_points_3::Triangulation& tr ) {
        deinit_shape();
        if( !has_neighborhood_radius() )
            estimate_mean_neighborhood( _mean_neighbors, _samples );
        _shape = Shape_of_points_3()( *tr, r2 );
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
	void construct_shape(InputIterator start, InputIterator end);

private:
    Triple ordered_facet_indices( const Facet& f ) const;

    /// Collect the triangles of one shell of the surface.
    /** A shell is a maximal collection of triangles that are
     *  connected and locally facing towards the same side of the
     *  surface.
     *  \tparam OutputIterator an output iterator for a collection
     *  of triangles. The iterator must point to a Triple type.
     *  \param c a cell touching a triangle of the shell from the
     *  outside of the object.
     *  \param li the index of the vertex of the cell opposite to
     *  the triangle touched.
     *  \param out an iterator to place to output the triangles.
     *  \sa collect_shell(const Facet& f, OutputIterator out).
     *  \sa collect_surface(OutputIterator out).
     */
	void collect_shell( Cell_handle c, unsigned int li );

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
     *  \tparam  the vertex base used in the internal data
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
	void collect_facets( Tag_true );
	void collect_facets( Tag_false );
    void collect_facets() { collect_facets( Shells() ); }

public:

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

    void collect_surface();

	
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
	void construct_surface(InputIterator start, InputIterator end, unsigned int iterations = 4);

    
	void construct_surface( unsigned int iterations = 4 );
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

#include <CGAL/Scale_space_surface_reconstructer_3_impl.h>

#endif // CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTER_H