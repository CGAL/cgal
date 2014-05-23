//A method to construct an interpolating surface.
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

/// Compute the triangular surface mesh interpolating a point set.
/** \ingroup PkgScaleSpaceReconstruction3Classes
 *  The surface mesh indicates the boundary of a sampled object by connecting
 *  the point sample using triangles. The point sample must be dense enough,
 *  such that any region that fits an empty ball with a fixed radius can be
 *  assumed to be ouside the object.
 *
 *  The surface cannot have holes and its triangles are all oriented towards
 *  the outside of the object sampled. This orientation is expressed using the
 *  'right-hand rule' on the ordered vertices of each triangle. Note that a
 *  surface may contain a triangle twice if the surface 'loops around' to
 *  incorporate both facings of the triangles.
 *
 *  The resulting surface need not be manifold and in many cases where the
 *  points sample the surface of an object, the computed surface will contain
 *  both an 'outward-facing' and a similar 'inward-facing' surface, with a thin
 *  volume between them.
 *
 *  If the object is not densely sampled or has disconnected components, a
 *  separate surface is constructed per connected component.
 *
 *  The surface is computed in three steps. First an appropriate neighborhood
 *  size is estimated or set manually. Then, a scale-space of the point set is
 *  constructed. This scale-space contains a smoothed representation of the
 *  point set, which makes the surface reconstruction problem less ill-posed.
 *  Finally, the surface of the points in the scale-space is computed. This
 *  surface in non-self-intersecting, i.e. its triangles cannot pairwise
 *  intersect in their interior.
 *
 *  The connectivity of the surface is stored as triples of indices to the
 *  point set. These steps can either be performed individually, or as a batch
 *  process.
 *
 *  The order of the points in the scale-space is the same as the order of the
 *  original points, meaning that the surface of the scale-space can
 *  interpolate the original point set by applying the indices of the surface
 *  triangles to the original point set. There are no guarantees on the
 *  topology of the resulting surface, e.g. the triangles of the surface may
 *  pairwise intersect in their interior. However, when using appropriate
 *  parameter settings for the number of iterations and neighborhood size the
 *  surface will generally not self-intersect. The appropriate parameter
 *  settings depend on the geometry of the point set.
 *
 *  The shape can be constructed either at a fixed scale, or at a dynamic
 *  scale. The first option is faster when constructing a single shape. It is
 *  undefined whether a shape with fixed scale may have its scale changed, but
 *  if so, this will likely require more computation time than changing a
 *  dynamic scale. In either case, it is possible to change the point set while
 *  maintaining the same scale.
 *
 *  The surface can be stored either as an unordered collection of triangles, 
 *  or as a collection of shells. A shell is a connected component of the 
 *  surface where connected facets are locally oriented towards the same side
 *  of the surface. Shells are separated by a (0,0,0) triple.
 *  
 *  \tparam Kernel the geometric traits class. It should be a model of
 *  DelaunayTriangulationTraits_3.
 *  \tparam Fixed_scale whether the shape is constructed at a fixed scale.
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
    typedef Search_traits_3< Kernel >                   Search_traits;
    typedef Orthogonal_k_neighbor_search< Search_traits >
                                                        Static_search;
    typedef Orthogonal_incremental_neighbor_search< Search_traits >
                                                        Dynamic_search;
    typedef typename Dynamic_search::Tree               Search_tree;

    typedef ::CGAL::Random                              Random;

    // Constructing the surface.
    typedef Shape_of_points_3< Kernel, Fixed_scale >    Shape_of_points_3;

    typedef typename Shape_of_points_3::Shape           Shape;
    typedef typename Shape_of_points_3::Triangulation   Triangulation;

    typedef typename Shape::Vertex_handle               Vertex_handle;
    typedef typename Shape::Cell_handle                 Cell_handle;
    typedef typename Shape::Facet                       Facet;
    
    typedef typename Shape::Vertex_iterator             Vertex_iterator;
    typedef typename Shape::Cell_iterator               Cell_iterator;
    typedef typename Shape::Facet_iterator              Facet_iterator;
    typedef typename Shape::Edge_iterator               Edge_iterator;

    typedef typename Shape::Finite_facets_iterator      Finite_facets_iterator;
    typedef typename Shape::Finite_facets_iterator      Finite_vertices_iterator;

    typedef typename Shape::All_cells_iterator          All_cells_iterator;

    typedef typename Shape::Classification_type         Classification_type;

public:
/// \name Types
/// \{
    typedef Shells                                      Collect_per_shell;      ///< Whether to collect the surface per shell.

	typedef typename Kernel::FT                         FT;                     ///< The number field type.

	typedef typename Kernel::Point_3                    Point;                  ///< The point type.
	typedef typename Kernel::Triangle_3                 Triangle;               ///< The triangle type.

#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                            Point_iterator;         ///< An iterator over the points.
    typedef const unspecified_type                      Const_point_iterator;   ///< A constant iterator over the points.
#else
    typedef typename Search_tree::iterator              Point_iterator;
    typedef typename Search_tree::const_iterator        Const_point_iterator;
#endif // DOXYGEN_RUNNING

    typedef Triple< unsigned int, unsigned int, unsigned int >
                                                        Triple;                 ///< A triple indicating a triangle of the surface.
private:
    typedef std::list< Triple >                         Tripleset;              ///< A collection of triples.

public:
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                            Triple_iterator;        ///< An iterator over the triples.
    typedef const unspecified_type                      Const_triple_iterator;  ///< A constant iterator over the triples.
#else
    typedef Tripleset::iterator                         Triple_iterator;
    typedef Tripleset::const_iterator                   Const_triple_iterator;
#endif // DOXYGEN_RUNNING
/// \}

private:
    class Finite_point_iterator;
    class In_surface_tester;
    typedef Filter_iterator< Facet_iterator, In_surface_tester >
                                                    Surface_facets_iterator;

private:
    Search_tree     _tree;              // To quickly search for nearest neighbors.

	Random          _generator;         // For sampling random points.

    unsigned int    _mean_neighbors;    // The number of nearest neighbors in the mean neighborhood.
    unsigned int    _samples;           // The number of sample points for estimating the neighborhood size.

    FT              _squared_radius;    // The squared neighborhood size.

    // The shape must be a pointer, because the alpha of a Fixed_alpha_shape_3
    // can only be set at construction and its assignment operator is private.
    // We want to be able to set the alpha after constructing the scale-space
    // reconstructer object.
	Shape*          _shape;

    // The surface. If the surface is collected per shell, consecutive triples
    // belong to the same shell and different shells are separated by a (0,0,0)
    // triple.
    Tripleset       _surface;

public:
/// \name Constructors
/// \{
    /// Construct a surface reconstructor that automatically estimates the neighborhood size.
    /** \param neighbors the number of neighbors a point's neighborhood should
     *  contain on average.
     *  \param samples the number of points sampled to estimate the neighborhood size.
     */
	Scale_space_surface_reconstructer_3( unsigned int neighbors, unsigned int samples );

    /// Construct a surface reconstructor with a given neighborhood size.
    /** \param sq_radius the squared radius of the neighborhood. If this value
     *  is negative when the point set is smoothed or when the surface is
     *  computed, the neighborhood size will be computed automatically.
     */
	Scale_space_surface_reconstructer_3( FT sq_radius );
/// \}
	~Scale_space_surface_reconstructer_3() { deinit_shape(); }

private:
    void deinit_shape() { if( _shape != 0 ) { delete _shape; _shape = 0; } }
    
    void clear_tree() { _tree.clear(); }
	void clear_surface() { if( has_shape() ) { _shape->clear(); } }
    
    // SURFACE COLLECTION
	// Once a facet is added to the surface, it is marked as handled.
	bool is_handled( Cell_handle c, unsigned int li ) const;
	inline bool is_handled( const Facet& f ) const { return is_handled( f.first, f.second ); }
	void set_handled( Cell_handle c, unsigned int li );
	inline void set_handled( Facet f ) { set_handled( f.first, f.second ); }
    
    // Get the indices of the points of the facet ordered to point
    // towardds the outside of the shape.
    Triple ordered_facet_indices( const Facet& f ) const;

    //  Collect the triangles of one shell of the surface.
	void collect_shell( Cell_handle c, unsigned int li );

    //  Collect the triangles of one shell of the surface.
	void collect_shell( const Facet& f ) {
		collect_shell( f.first, f.second );
	}

    //  Collect the triangles of the complete surface.
	void collect_facets( Tag_true );
	void collect_facets( Tag_false );
    void collect_facets() { 
        if( !has_neighborhood_radius() )
            estimate_neighborhood_radius();
        collect_facets( Shells() );
    }

public:
/// \name Accessors
/// \{
    /// Gives the squared radius of the neighborhood.
    /** The neighborhood is used by
     *  [advance_scale_space](\ref advance_scale_space) and
     *  [construct_scale_space](\ref construct_scale_space) to compute the
     *  scale-space and by [reconstruct_surface](\ref reconstruct_surface) to
     *  construct the shape of the scale-space.
     *
     *  \return the squared radius of the neighborhood, or -1 if the
     *  neighborhood has not yet been set.
     *
     *  \sa advance_scale_space(unsigned int iterations).
     *  \sa construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations).
     *  \sa reconstruct_surface(unsigned int iterations).
     *  \sa reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     */
    FT get_neighborhood_squared_radius() const { return _squared_radius; }

    /// Gives the mean number of neighbors an estimated neighborhood should contain.
    /** This number is only used if the neighborhood size has not been set
     *  manually.
     *
     *  When the neighborhood size is estimated, it should on average contain
     *  this many neighbors, not counting the neighborhood center.
     *
     *  \return the number of neighbors a neighborhood ball centered on a point
     *  should contain on average when the radius is estimated, not counting
     *  the point itself.
     *
     *  \sa set_mean_neighbors(unsigned int neighbors).
     *  \sa has_neighborhood_radius().
     *  \sa get_neighborhood_samples().
     *  \sa estimate_neighborhood_radius().
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end).
     */
    unsigned int get_mean_neighbors() const { return _mean_neighbors; }

    /// Gives the number of sample points the neighborhood estimation uses.
    /** This number is only used if the neighborhood size has not been set
     *  manually.
     *
     *  If the number of samples is larger than the point cloud, every point is
     *  used and the optimal neighborhood size is computed exactly instead of
     *  estimated.
     *
     *  \return the number of points sampled for neighborhood estimation.
     *
     *  \sa set_neighborhood_samples(unsigned int samples).
     *  \sa has_neighborhood_radius().
     *  \sa get_mean_neighbors().
     *  \sa estimate_neighborhood_radius().
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end).
     */
    unsigned int get_neighborhood_samples() const { return _samples; }
    
    /// Get the number of triangles of the surface.
    unsigned int get_surface_size() const { return (unsigned int)_surface.size(); }
/// \}

private:
    //  Get the shape of the scale space.
    /*  If the shape does not exist, it is  constructed first.
     *  \return the shape of the scale space.
     */
    const Shape& get_shape() const;

public:
/// \name Mutators
/// \{
    /// Reset the scale-space surface reconstructer.
    void clear() {
		clear_tree();
        clear_surface();
    }

    /// Sets the radius of the neighborhood.
    /** The neighborhood is used by
     *  [advance_scale_space](\ref advance_scale_space) and
     *  [construct_scale_space](\ref construct_scale_space) to compute the
     *  scale-space and by [reconstruct_surface](\ref reconstruct_surface) to
     *  construct the shape of the scale-space.
     *
     *  \param radius the radius of the neighborhood.
     *
     *  \sa get_neighborhood_squared_radius().
     *  \sa has_neighborhood_radius().
     *  \sa set_neighborhood_squared_radius(const FT& sq_radius).
     *  \sa advance_scale_space(unsigned int iterations).
     *  \sa construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations).
     *  \sa reconstruct_surface(unsigned int iterations).
     *  \sa reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     */
    void set_neighborhood_radius( const FT& radius ) {
        _squared_radius = radius * radius;
		if( has_shape() )
            Shape_of_points_3().set_scale( _shape, _squared_radius );
    }
    
    /// Sets the squared radius of the neighborhood.
    /** The neighborhood is used by
     *  [advance_scale_space](\ref advance_scale_space) and
     *  [construct_scale_space](\ref construct_scale_space) to compute the
     *  scale-space and by [reconstruct_surface](\ref reconstruct_surface) to
     *  construct the shape of the scale-space.
     *
     *  \param sq_radius the squared radius of the neighborhood. If sq_radius
     *  is negative, the neighborhood size will be estimated automatically.
     *
     *  \sa get_neighborhood_squared_radius().
     *  \sa has_neighborhood_radius().
     *  \sa set_neighborhood_radius(const FT& sq_radius).
     *  \sa advance_scale_space(unsigned int iterations).
     *  \sa construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations).
     *  \sa reconstruct_surface(unsigned int iterations).
     *  \sa reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     */
    void set_neighborhood_squared_radius( const FT& sq_radius ) {
        _squared_radius = sq_radius;
		if( has_shape() )
            Shape_of_points_3().set_scale( _shape, _squared_radius );
    }

    /// Sets the mean number of neighbors an estimated neighborhood should contain.
    /** This number is only used if the neighborhood size has not been set
     *  manually.
     *
     *  When the neighborhood size is estimated, it should on average contain
     *  this many neighbors, not counting the neighborhood center.
     *
     *  \param the number of neighbors a neighborhood ball centered on a point
     *  should contain on average when the radius is estimated, not counting
     *  the point itself.
     *
     *  \sa get_mean_neighbors().
     *  \sa has_neighborhood_radius().
     *  \sa set_neighborhood_samples(unsigned int samples).
     */
    void set_mean_neighbors( unsigned int neighbors ) { _mean_neighbors = neighbors; }
    
    /// Sets the number of sample points the neighborhood estimation uses.
    /** This number is only used if the neighborhood size has not been set
     *  manually.
     *
     *  If the number of samples is larger than the point cloud, every point is
     *  used and the optimal neighborhood size is computed exactly instead of
     *  estimated.
     *
     *  \return the number of points sampled for neighborhood estimation.
     *
     *  \sa get_neighborhood_samples().
     *  \sa has_neighborhood_radius().
     *  \sa set_mean_neighbors(unsigned int neighbors).
     *  \sa estimate_neighborhood_radius().
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end).
     */
    void set_neighborhood_samples( unsigned int samples ) { _samples = samples; }
    
    /// Insert a point into the scale-space.
    /** Note that inserting the point does not automatically construct or
     *  update the surface.
     *
     *  In order to construct the surface, call
     *  reconstruct_surface(unsigned int iterations).
     *
     *  \param p the point to insert.
     *
     *  \sa insert_points(InputIterator begin, InputIterator end).
     */
	void insert_point( const Point& p ) {
		_tree.insert( p );
	}

    /// Insert a collection of points into the scale-space.
    /** Note that inserting the points does not automatically construct or
     *  update the surface.
     *
     *  In order to construct the surface, either call
     *  reconstruct_surface(unsigned int iterations) after inserting the
     *  points, or insert the points using
     *  reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     *
     *  \tparam InputIterator an iterator over a collection of points. The
     *  iterator must point to a Point type.
     *  \param begin an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *
     *  \sa insert_point(const Point& p).
     *  \sa construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations).
     *  \sa reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     */
	template < class InputIterator >
	void insert_points( InputIterator begin, InputIterator end ) {
		_tree.insert( begin, end );
	}
/// \}
    
public:
/// \name Query
/// \{
    /// Check whether the neighborhood ball radius has been set.
    /** The radius can be set manually, or estimated automatically.
     *
     *  \return true iff the radius has been either set manually or estimated.
     *
     *  \sa set_neighborhood_radius().
     *  \sa set_neighborhood_squared_radius().
     *  \sa estimate_neighborhood_radius().
     */
    bool has_neighborhood_radius() const {
        return sign( _squared_radius ) == POSITIVE;
    }
/// \}

private:
    //  Check whether the shape has been constructed.
    /*  The shape contains the structure of the point cloud.
     *
     *  Until the shape is constructed, the surface is undefined.
     *  \return true if the shape exists and false otherwise.
     *  \sa construct_shape(InputIterator begin, InputIterator end).
     */
	bool has_shape() const { return _shape != 0; }
    
public:
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
/// \name Neighborhood Estimation
/// \{
    /// Estimate the neighborhood size.
    /** The neighborhood size is expressed as the radius of the smallest ball
     *  centered on a point such that the ball contains at least a specified
     *  number of points, not counting the point itself.
     *
     *  The neighborhood size is set to the mean of these radii, taken over a
     *  number of point samples.
     *
     *  This method will be called by the scale-space and surface construction
     *  methods if the neighborhood size is not set when they are called.
     *
     *  \return the estimated neighborhood size.
     *
     *  \sa set_mean_neighbors(unsigned int neighbors).
     *  \sa set_neighborhood_samples(unsigned int samples).
     *  \sa estimate_neighborhood_radius(unsigned int neighbors, unsigned int samples).
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end).
     *  \sa advance_scale_space(unsigned int iterations).
     *  \sa reconstruct_surface(unsigned int iterations).
     */
    inline FT estimate_neighborhood_radius() {
        return estimate_neighborhood_radius( get_mean_neighbors(), get_neighborhood_samples() );
    }

    /// Estimate the neighborhood size based on a number of sample points.
    /** The neighborhood size is expressed as the radius of the smallest ball
     *  centered on a point such that the ball contains at least a specified
     *  number of points, not counting the point itself.
     *
     *  The neighborhood size is set to the mean of these radii, taken over a
     *  number of point samples.
     *
     *  \param neighbors the number of neighbors a point's neighborhood
     *  should contain on average, not counting the point itself.
     *  \param samples the number of points sampled to estimate the
     *  neighborhood size.
     *  \return the estimated neighborhood size.
     *
     *  \sa estimate_neighborhood_radius().
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples).
     */
	FT estimate_neighborhood_radius( unsigned int neighbors, unsigned int samples );
    
    /// Estimate the neighborhood size of a collection of points.
    /** The neighborhood size is expressed as the radius of the smallest ball
     *  centered on a point such that the ball contains at least a specified
     *  number of points, not counting the point itself.
     *
     *  The neighborhood size is set to the mean of these radii, taken over a
     *  number of point samples.
     *  
     *  This method is equivalent to running
     *  [insert_points(begin, end)](\ref insert_points) followed by
     *  estimate_neighborhood_radius().
     *
     *  This method will be called by the scale-space and surface construction
     *  methods if the neighborhood size is not set when they are called.
     *
     *  \tparam InputIterator an iterator over a collection of points. The
     *  iterator must point to a Point type.
     *  \param begin an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \return the estimated neighborhood size.
     *
     *  \sa set_mean_neighbors(unsigned int neighbors).
     *  \sa set_neighborhood_samples(unsigned int samples).
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples).
     *  \sa estimate_neighborhood_radius().
     *  \sa insert_points(InputIterator begin, InputIterator end).
     *  \sa construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations).
     *  \sa reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     */
	template < class InputIterator >
    FT estimate_neighborhood_radius( InputIterator begin, InputIterator end ) {
        return estimate_neighborhood_radius( begin, end, get_mean_neighbors(), get_neighborhood_samples() );
    }
    
    /// Estimate the neighborhood size of a collection of points based on a number of sample points.
    /** The neighborhood size is expressed as the radius of the smallest ball
     *  centered on a point such that the ball contains at least a specified
     *  number of points, not counting the point itself.
     *
     *  The neighborhood size is set to the mean of these radii, taken over a
     *  number of point samples.
     *  
     *  This method is equivalent to running
     *  [insert_points(begin, end)](\ref insert_points) followed by
     *  [estimate_neighborhood_radius(neighbors, samples)](\ref estimate_neighborhood_radius).
     *
     *  \tparam InputIterator an iterator over a collection of points. The
     *  iterator must point to a Point type.
     *  \param begin an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \param neighbors the number of neighbors a point's neighborhood
     *  should contain on average, not counting the point itself.
     *  \param samples the number of points sampled to estimate the
     *  neighborhood size.
     *  \return the estimated neighborhood size.
     *
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end).
     *  \sa estimate_neighborhood_radius(unsigned int neighbors, unsigned int samples).
     */
	template < class InputIterator >
	FT estimate_neighborhood_radius( InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples );
/// \}
    
/// \name Scale-space
/// \{
    /// Advance the scale-space by a number of iterations.
    /** Each iteration the scale-space is advanced, a higher scale-space is
     *  computed. At a higher scale, the scale-space contains a smoother
     *  representation of the point set.
     *
     *  In case the scale-space is not at the scale of the original point set,
     *  calling [advance_scale_space(iterations)](\ref advance_scale_space)
     *  with `iterations > 0` will advance the scale-space further.
     *
     *  If the neighborhood has not been set before, it is automatically
     *  estimated.
     *
     *  \param iterations the number of iterations to perform. If 0, nothing
     *  happens.
     *
     *  \sa construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations).
     *  \sa estimate_neighborhood_radius().
     *  \sa reconstruct_surface(unsigned int iterations).
     */
	void advance_scale_space( unsigned int iterations = 1 );
    
    /// Construct a scale-space of a collection of points.
    /** If the scale-space is advanced, the scale-space is computed at a higher
     *  scale. At a higher scale, the scale-space contains a smoother
     *  representation of the point set.
     *
     *  If the neighborhood has not been set before, it is automatically
     *  estimated.
     *
     *  This method is equivalent to running
     *  [insert_points(begin, end)](\ref insert_points) followed by
     *  [advance_scale_space(iterations)](\ref advance_scale_space).
     *
     *  \tparam InputIterator an iterator over a collection of points. The
     *  iterator must point to a Point type.
     *  \param begin an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \param iterations the number of iterations to perform. If 0, nothing
     *  happens.
     *
     *  \sa insert_points(InputIterator begin, InputIterator end).
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end).
     *  \sa advance_scale_space(unsigned int iterations).
     *  \sa reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     */
	template < class InputIterator >
	void construct_scale_space( InputIterator begin, InputIterator end, unsigned int iterations = 1 ) {
        clear_tree();
		insert_points( begin, end );
		advance_scale_space( iterations );
	}
/// \}
private:
    // Construct the scale-space from a triangulation.
    void construct_scale_space( Triangulation& tr ) {
        insert_points( tr.finite_vertices_begin(), tr.finite_vertices_end() );
    }
    
private:
/// \name Shape
/// \{
    /// Construct the shape of the scale-space.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  If the neighborhood has not been set before, it is automatically
     *  estimated.
     */
    void construct_shape() {
		construct_shape( scale_space_begin(), scale_space_end() );
    }
    
    /// Construct the shape from an existing triangulation.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  Note that this does not set the current scale-space.
     *  To set this as well, use `construct_scale_space(Triangulation& tr)`.
     *
     *  If the neighborhood has not been set before, it is automatically
     *  estimated.
     *  \param tr the triangulation to construct the shape of.
     */
    void construct_shape(Triangulation& tr ) {
        deinit_shape();
        if( !has_neighborhood_radius() )
            estimate_neighborhood_radius();
        _shape = Shape_of_points_3()( *tr, _squared_radius );
	}

    /// Construct the shape from a collection of points.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  Note that this does not set the current scale-space.
     *  To set this as well, use `insert_points( InputIterator begin, InputIterator end )`.
     *
     *  If the neighborhood has not been set before, it is automatically
     *  estimated.
     *  \tparam InputIterator an iterator over the point sample.
     *  The iterator must point to a Point type.
     *  \param begin an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa is_constructed().
     */
	template < class InputIterator >
	void construct_shape( InputIterator begin, InputIterator end );
    
    // Collect the surface mesh from the shape.
    // If the sahep does not yet exist, it is constructed.
    void collect_surface();
/// \}

public:
/// \name Surfaces
/// \{
    /// Construct a triangle mesh from the scale-space.
    /** The order of the points in the scale-space is the same as the order of
     *  the original points, meaning that the surface of the scale-space can
     *  interpolate the original point set by applying the indices of the
     *  surface triangles to the original point set.
     *
     *  After construction, the triangles of the surface can be iterated using
     *  surface_begin() and surface_end().
     *
     *  If the neighborhood has not been set before, it is automatically
     *  estimated.
     *
     *  \param iterations the number of scale-space advancement iterations to
     *  apply. If 0, the current scale-space is used.
     *
     *  \sa reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations).
     *  \sa estimate_neighborhood_radius().
     *  \sa advance_scale_space(unsigned int iterations).
     */
	void reconstruct_surface( unsigned int iterations = 0 );
    
    /// Construct a surface mesh from the scale-space of a collection of points.
    /**  The order of the points in the scale-space is the same as the order of
     *  the original points, meaning that the surface of the scale-space can
     *  interpolate the original point set by applying the indices of the
     *  surface triangles to the original point set.
     *
     *  After construction, the triangles of the surface can be iterated using
     *  surface_begin() and surface_end().
     *
     *  If the neighborhood has not been set before, it is automatically
     *  estimated.
     *
     *  \tparam InputIterator an iterator over the point collection.
     *  The iterator must point to a Point type.
     *  \param begin an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \param iterations the number of scale-space advancement iterations to
     *  apply. If 0, the current scale-space is used.
     *  
     *  \sa reconstruct_surface(unsigned int iterations).
     *  \sa insert_points(InputIterator begin, InputIterator end).
     *  \sa estimate_neighborhood_radius(InputIterator begin, InputIterator end).
     *  \sa construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations).
     */
	template < class InputIterator >
	void reconstruct_surface( InputIterator begin, InputIterator end, unsigned int iterations = 0 );
/// \}
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