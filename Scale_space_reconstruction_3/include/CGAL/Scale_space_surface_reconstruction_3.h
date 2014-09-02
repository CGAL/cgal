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

#ifndef CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_H
#define CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_H

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/utility.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Random.h>

#include <CGAL/Scale_space_reconstruction_3/internal/check3264.h>
#include <CGAL/Scale_space_reconstruction_3/Shape_construction_3.h>

#include <CGAL/Scale_space_reconstruction_3/Weighted_PCA_projection_3.h>

#include <boost/mpl/and.hpp>

namespace CGAL {

/// computes a triangulated surface mesh interpolating a point set.
/** \ingroup PkgScaleSpaceReconstruction3Classes
 *  The scale-space surface reconstruction method is twofold. Firstly, a
 *  scale-space of the point set is constructed. This scale-space contains a
 *  smoothed representation of the point set, which makes the surface
 *  reconstruction problem less ill-posed. Then, a triangulated surface mesh of
 *  the points in the scale-space is computed. This mesh is represented as
 *  triples of indices of the point set. Neighboring triples have consistent
 *  orientation, as expressed using the '<em>right-hand rule</em>' on the
 *  ordered vertices of the triple.
 * 
 *  The points maintain their original order in the scale-space. This means it
 *  is straightforward to revert the scale-space on the triangulated surface
 *  mesh. In essence, the method can reconstruct a smoothed surface or a
 *  surface interpolating the original points. The only change is whether to
 *  apply the index triples to the scale-space or to the original point set.
 *
 *  When applied to the scale-space, the surface mesh is non-self-intersecting.
 *  The interior of the triangles cannot pairwise intersect in a line segment.
 *  However, the surface does not need to be manifold. An edge may be incident
 *  to more than two triangles and triangles may overlap exactly if the regions
 *  on both sides of the triangle are not part of the object. Note that we
 *  count overlapping triangles with opposite orientations as separate
 *  triangles. In many cases where the points sample the surface of an object,
 *  the computed surface will contain both an '<em>outward-facing</em>' and a
 *  similar '<em>inward-facing</em>' surface, with a thin volume between them.
 *
 *  The surface mesh will not have edges incident to only one triangle or
 *  holes, loops of such edges, and the triangles are all oriented towards the
 *  outside of the object sampled. If the point sample has '<em>holes</em>', it
 *  is likely that the surface mesh will contain overlapping triangles with
 *  opposite orientation touching this hole.
 *
 *  If the object is not densely sampled or has disconnected components, the
 *  surface may have several disconnected components. The surface may be
 *  presented as an unordered collection of triangles, or as a collection
 *  sorted per \em shell. A shell is a collection of connected triangles that
 *  are locally oriented towards the same side of the surface.
 *
 *  When applied to the original points, we are unable to guarantee the same
 *  topology of the surface. The triangles of this surface may pairwise
 *  intersect in their interior and the surface could have boundary edges.
 *  However, when using appropriate parameter settings for the number of
 *  iterations and neighborhood size the surface will generally not
 *  self-intersect. The appropriate parameter settings depend on the geometry
 *  of the point set and generally need to be fine-tuned per data set.
 *
 *  Both the smoothing operator and the mesh reconstruction assume that points
 *  near each other belong to the same part of the object. This is expressed in
 *  the notion of balls with a fixed size, the neighborhood radius. If such a
 *  ball contains multiple points, these points are near each other and will
 *  influence each other while advancing the scale-space. If such a ball is
 *  empty, it lies outside the object. Note that \em outside is based on
 *  regions empty of points, not on whether a volume is enclosed by the
 *  surface.
 *
 *  The scale-space is constructed by projecting each point to the
 *  '<em>density</em>'-weighted principal component analysis (PCA) of the local
 *  (\f$ \delta \f$-distance) neighborhood. If the point set was sampled from a
 *  surface for which any high-frequency deformation and sampling noise is
 *  smaller than the neighborhood size, the scale is coarse enough for mesh
 *  reconstruction after a few iterations of advancing the scale-space.
 *
 *  The mesh reconstruction of the scale-space is performed by filtering a 3D
 *  \f$ \alpha \f$-shape. The result is returned as a collection of triples on
 *  the indices of the points of the surface triangles. Recall that this
 *  collection may be sorted per shell, where a shell is a collection of
 *  connected triangles that are locally oriented towards the same side of the
 *  surface.
 *
 *  The reconstruction method requires a neighborhood radius, related to the
 *  resolution of the data. This parameter can be estimated through statistical
 *  analysis. Specifically, we use a kD-tree to estimate the mean distance to
 *  the n-th nearest neighbor and we use this distance as an approximator for
 *  the resolution.
 *
 *  The method generally works well as long as this neighborhood radius is not
 *  too small and the number of scale-space advancement iterations necessary to
 *  reduce the high-frequency signals does not disturb the topology of the data
 *  so much that applying the surface connectivity to the original point set
 *  produces too many self-intersections. As a general rule of thumb, a
 *  neighborhood containg 30 points on average provides a good estimate for the
 *  radius and 4 iterations of smoothing proved a nice scale-space.
 *
 *  This class stores several of the (intermediate) results. This makes it
 *  easier and more efficient to adjust the parameter settings based on
 *  preliminary results, or to further advance the scale-space to improve the
 *  results. The class stores the current scale-space and the reconstructed
 *  surface, possibly with iterators over the shells.
 *
 *  The class also stores the parameters for estimating the optimal
 *  neighborhood radius and either the lastest estimate or the manually set
 *  radius. This way, the radius can be estimated (again) whenever necessary.
 *  Also note that both advancing the scale-space and reconstructing the
 *  surface use this radius. By changing or re-estimating the radius between
 *  these operations, they can use separate parameter settings.
 *
 *  The shape can be constructed either at a fixed scale, or at a dynamic
 *  scale. When constructing the surface for exactly one neighborhood radius,
 *  it is faster to set `FixedScale` to `Tag_true`. If the correct neighborhood
 *  radius should be changed or estimated multiple times, it is faster to set
 *  `FixedScale` to `Tag_false`.
 *
 *  It is undefined whether a shape with fixed scale may have its scale
 *  changed, but if so, this will likely require more computation time than
 *  changing a dynamic scale. In either case, it is possible to change the
 *  point set while maintaining the same scale.
 *
 *  The surface can be stored either as an unordered collection of triangles, 
 *  or as a collection of shells. A shell is a connected component of the 
 *  surface where connected facets are locally oriented towards the same side
 *  of the surface.
 *  
 *  \tparam DelaunayTriangulationTraits_3 is the geometric traits class.
 *  Generally, `Exact_predicates_inexact_constructions_kernel` is preferred.
 *  \tparam FixedScale determines whether the shape is constructed at a fixed
 *  scale. It must be a `Boolean_tag` type. The default value is `Tag_true`.
 *  \tparam OrderShells determines whether to collect the surface per shell. It
 *  must be a `Boolean_tag` type. The default value is `Tag_true`.
 *  \tparam WeightedPCAProjection_3 is the type of weighted PCA to use. The
 *  default value is `Weighted_PCA_projection_3<DelaunayTriangulationTraits_3>`.
 */
#ifdef DOXYGEN_RUNNING
template < class DelaunayTriangulationTraits_3, class FixedScale, class OrderShells, class WeightedPCAProjection_3 >
#else
template < class Gt, class FS = Tag_true, class OS = Tag_true, class Ct = Parallel_tag, class WPCA = Weighted_PCA_projection_3< Gt > >
#endif
class Scale_space_surface_reconstruction_3 {
public:
    typedef FS                                          FixedScale;
    typedef OS                                          OrderShells;
    typedef Ct                                          Concurrency_tag;

private:
    // Searching for neighbors.
    typedef Search_traits_3< Gt >                       Search_traits;
    typedef Orthogonal_k_neighbor_search< Search_traits >
                                                        Static_search;
    typedef Orthogonal_incremental_neighbor_search< Search_traits >
                                                        Dynamic_search;
    typedef typename Dynamic_search::Tree               Search_tree;

    typedef CGAL::Random                                Random;

    // Constructing the surface.
    typedef CGAL::Shape_construction_3< Gt, FS >        Shape_construction_3;

    typedef typename Shape_construction_3::Shape           Shape;
    typedef typename Shape_construction_3::Triangulation   Triangulation;

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
#ifdef DOXYGEN_RUNNING
    typedef typename DelaunayTriangulationTraits_3::FT          FT;             ///< defines the number field type.

	typedef typename DelaunayTriangulationTraits_3::Point_3     Point;          ///< defines the point type.
	typedef typename DelaunayTriangulationTraits_3::Triangle_3  Triangle;       ///< defines the triangle type.

    typedef unspecified_type                            Point_iterator;         ///< defines an iterator over the points.
    typedef const unspecified_type                      Const_point_iterator;   ///< defines a constant iterator over the points.
#else // DOXYGEN_RUNNING
    typedef typename Gt::FT                             FT;

	typedef typename Gt::Point_3                        Point;
	typedef typename Gt::Triangle_3                     Triangle;

    typedef typename Search_tree::iterator              Point_iterator;
    typedef typename Search_tree::const_iterator        Const_point_iterator;
#endif // DOXYGEN_RUNNING

    typedef CGAL::cpp11::array< unsigned int, 3 >       Triple;                 ///< defines a triple of point indices indicating a triangle of the surface.
private:
    typedef std::list< Triple >                         Tripleset;              ///< defines a collection of triples.
    // Note that this is a list for two reasons: iterator validity for the shell iterators, and memory requirements for the expected huge collections.

public:
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                            Triple_iterator;        ///< defines an iterator over the triples.
    typedef const unspecified_type                      Const_triple_iterator;  ///< defines a constant iterator over the triples.
#else // DOXYGEN_RUNNING
    typedef Tripleset::iterator                         Triple_iterator;
    typedef Tripleset::const_iterator                   Const_triple_iterator;
#endif // DOXYGEN_RUNNING

/// \}

private:
    typedef std::vector< Triple_iterator >              TripleIterSet;
    
private:
    class Finite_point_iterator;
    class In_surface_tester;
    typedef Filter_iterator< Facet_iterator, In_surface_tester >
                                                    Surface_facets_iterator;

    // Parallel processing functors.
    class ComputeNN;
    class AdvanceSS;

private:
    Search_tree     _tree;              // To quickly search for nearest neighbors.

	Random          _generator;         // For sampling random points.

    unsigned int    _mean_neighbors;    // The number of nearest neighbors in the mean neighborhood.
    unsigned int    _samples;           // The number of sample points for estimating the neighborhood radius.

    FT              _squared_radius;    // The squared neighborhood radius.

    // The shape must be a pointer, because the alpha of a Fixed_alpha_shape_3
    // can only be set at construction and its assignment operator is private.
    // We want to be able to set the alpha after constructing the scale-space
    // reconstructer object.
	Shape*          _shape;

    // The surface. If the surface is collected per shell, the triples of the
    // same shell are stored consecutively.
    Tripleset       _surface;

    // The shells can be accessed through iterators to the surface.
    TripleIterSet   _shells;

public:
/// \name Constructors
/// \{
    /// constructs a surface reconstructor that will automatically estimate the neighborhood radius.
    /** \param neighbors is the number of neighbors a point's neighborhood should
     *  contain on average.
     *  \param samples is the number of points sampled to estimate the
     *  neighborhood radius.
     */
	Scale_space_surface_reconstruction_3( unsigned int neighbors, unsigned int samples );

    /// constructs a surface reconstructor with a given neighborhood radius.
    /** \param sq_radius is stored as the squared radius of the neighborhood.
     *
     *  \note If the neighborhood squared radius is negative when the point set
     *  is smoothed or when the surface is computed, the neighborhood radius
     *  will be computed automatically.
     */
	Scale_space_surface_reconstruction_3( FT sq_radius );

/// \}
	~Scale_space_surface_reconstruction_3() { deinit_shape(); }

private:
    void deinit_shape() { if( _shape != 0 ) { delete _shape; _shape = 0; } }

    void clear_tree() { _tree.clear(); }
	void clear_surface() { _surface.clear(); deinit_shape(); }
    
    // SURFACE COLLECTION
	// Once a facet is added to the surface, it is marked as handled.
	bool is_handled( Cell_handle c, unsigned int li ) const;
	inline bool is_handled( const Facet& f ) const { return is_handled( f.first, f.second ); }
	void mark_handled( Cell_handle c, unsigned int li );
	inline void mark_handled( Facet f ) { mark_handled( f.first, f.second ); }
    
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
        collect_facets( OS() );
    }

private:
    //  Get the shape of the scale space.
    /*  If the shape does not exist, it is  constructed first.
     *  \return the shape of the scale space.
     */
    const Shape& shape() const;

public:
/// \name Point Set Manipulation
/// \{
    /// inserts a point into the scale-space.
    /** \param p is the point to insert.
     *
     *  \note Inserting the point does not automatically construct or
     *  update the surface.
     *
     *  \note In order to construct the surface, call
     *  `reconstruct_surface(unsigned int iterations)`.
     *
     *  \warning Inserting a new point may invalidate the neighborhood radius
     *  if it was previously estimated.
     *
     *  \sa `add_points(InputIterator begin, InputIterator end)`.
     */
	void add_point( const Point& p ) {
		_tree.insert( p );
	}

    /// inserts a collection of points into the scale-space.
    /** \tparam InputIterator is an iterator over the point collection.
     *  The iterator must point to a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *
     *  \note Inserting the points does not automatically construct or
     *  update the surface.
     *
     *  \note In order to construct the surface, either call
     *  `reconstruct_surface(unsigned int iterations)` after inserting the
     *  points, or insert the points using
     *  `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     *
     *  \warning Inserting new points may invalidate the neighborhood radius if
     *  it was previously estimated.
     *
     *  \sa `add_point(const Point& p)`.
     *  \sa `construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations)`.
     *  \sa `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	void add_points( InputIterator begin, InputIterator end ) {
#else // DOXYGEN_RUNNING
	void add_points( InputIterator begin, InputIterator end,
                     typename boost::enable_if<
                        boost::is_convertible< typename std::iterator_traits<InputIterator>::value_type,
                                               Point > >::type* = NULL ) {
#endif // DOXYGEN_RUNNING
		_tree.insert( begin, end );
	}
    
    /// clears the scale-space surface reconstruction data.
    /** This includes discarding the surface, the scale-space, and any
     *  estimation of the neighborhood radius.
     */
    void clear() {
		clear_tree();
        clear_surface();
        _squared_radius = -1;
    }

/// \}
    
private:
    //  checks whether the shape has been constructed.
    /*  The shape contains the structure of the point cloud.
     *
     *  Until the shape is constructed, the surface is undefined.
     *  \return true if the shape exists and false otherwise.
     *  \sa construct_shape(InputIterator begin, InputIterator end).
     */
	bool has_shape() const { return _shape != 0; }

public:
/// \name Neighborhood Size Estimation
/// \{
    /// checks whether the neighborhood radius has been set.
    /** The radius can be set manually, or estimated automatically.
     *
     *  \return `true` iff the radius has been either set manually or estimated.
     *
     *  \sa `set_neighborhood_squared_radius()`.
     *  \sa `estimate_neighborhood_radius()`.
     */
    bool has_neighborhood_radius() const {
        return sign( _squared_radius ) == POSITIVE;
    }

    /// gives the squared radius of the neighborhood.
    /** The neighborhood radius is used by
     *  <code>[advance_scale_space](\ref advance_scale_space)</code> and
     *  <code>[construct_scale_space](\ref construct_scale_space)</code> to
     *  compute the scale-space and by
     *  <code>[reconstruct_surface](\ref reconstruct_surface)</code> to
     *  construct the shape of the scale-space.
     *
     *  \return the squared radius of the neighborhood, or -1 if the
     *  neighborhood radius has not yet been set.
     *
     *  \sa `advance_scale_space(unsigned int iterations)`.
     *  \sa `construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations)`.
     *  \sa `reconstruct_surface(unsigned int iterations)`.
     *  \sa `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
    FT neighborhood_squared_radius() const { return _squared_radius; }

    /// gives the mean number of neighbors an estimated neighborhood should contain.
    /** This number is only used if the neighborhood radius has not been set
     *  manually.
     *
     *  When the neighborhood radius is estimated, it should on average contain
     *  this many neighbors, not counting the neighborhood center.
     *
     *  \return the number of neighbors a neighborhood ball centered on a point
     *  should contain on average when the radius is estimated, not counting
     *  the point itself.
     *
     *  \sa `set_mean_number_of_neighbors(unsigned int neighbors)`.
     *  \sa `has_neighborhood_radius()`.
     *  \sa `neighborhood_sample_size()`.
     *  \sa `estimate_neighborhood_radius()`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end)`.
     */
    unsigned int mean_number_of_neighbors() const { return _mean_neighbors; }

    /// gives the number of sample points the neighborhood estimation uses.
    /** This number is only used if the neighborhood radius has not been set
     *  manually.
     *
     *  If the number of samples is larger than the point cloud, every point is
     *  used and the optimal neighborhood radius is computed exactly instead of
     *  estimated.
     *
     *  \return the number of points sampled for neighborhood estimation.
     *
     *  \sa `set_neighborhood_sample_size(unsigned int samples)`.
     *  \sa `has_neighborhood_radius()`.
     *  \sa `mean_number_of_neighbors()`.
     *  \sa `estimate_neighborhood_radius()`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end)`.
     */
    unsigned int neighborhood_sample_size() const { return _samples; }
    
    /// sets the squared radius of the neighborhood.
    /** The neighborhood radius is used by
     *  <code>[advance_scale_space](\ref advance_scale_space)</code> and
     *  <code>[construct_scale_space](\ref construct_scale_space)</code> to
     *  compute the scale-space and by
     *  <code>[reconstruct_surface](\ref reconstruct_surface)</code> to
     *  construct the shape of the scale-space.
     *
     *  \param sq_radius is stored as the squared radius of the neighborhood.
     *
     *  \note If the neighborhood squared radius is negative when the point set
     *  is smoothed or when the surface is computed, the neighborhood radius
     *  will be computed automatically.
     *
     *  \warning If the surface was already constructed, changing the
     *  neighborhood radius will automatically adjust the surface.
     *
     *  \sa `neighborhood_squared_radius()`.
     *  \sa `has_neighborhood_radius()`.
     *  \sa `advance_scale_space(unsigned int iterations)`.
     *  \sa `construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations)`.
     *  \sa `reconstruct_surface(unsigned int iterations)`.
     *  \sa `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
    void set_neighborhood_squared_radius( const FT& sq_radius ) {
        _squared_radius = sq_radius;
		if( has_neighborhood_radius() && has_shape() )
            Shape_construction_3().change_scale( _shape, _squared_radius );
    }

    /// sets the mean number of neighbors an estimated neighborhood should contain.
    /** This number is only used if the neighborhood radius has not been set
     *  manually.
     *
     *  When the neighborhood radius is estimated, it should on average contain
     *  this many neighbors, not counting the neighborhood center.
     *
     *  \param neighbors is the number of neighbors a neighborhood ball centered on a point
     *  should contain on average when the radius is estimated, not counting
     *  the point itself.
     *
     *  \note This does not start the estimation process.
     *
     *  \sa `mean_number_of_neighbors()`.
     *  \sa `has_neighborhood_radius()`.
     *  \sa `set_neighborhood_sample_size(unsigned int samples)`.
     */
    void set_mean_number_of_neighbors( unsigned int neighbors ) { _mean_neighbors = neighbors; }
    
    /// sets the number of sample points the neighborhood estimation uses.
    /** This number is only used if the neighborhood radius has not been set
     *  manually.
     *
     *  If the number of samples is larger than the point cloud, every point is
     *  used and the optimal neighborhood radius is computed exactly instead of
     *  estimated.
     *
     *  \param samples is the number of points to sample for neighborhood
     *  estimation.
     *
     *  \note This does not start the estimation process.
     *
     *  \sa `neighborhood_sample_size()`.
     *  \sa `has_neighborhood_radius()`.
     *  \sa `set_mean_number_of_neighbors(unsigned int neighbors)`.
     *  \sa `estimate_neighborhood_radius()`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end)`.
     */
    void set_neighborhood_sample_size( unsigned int samples ) { _samples = samples; }

    /// estimates the neighborhood radius.
    /** This method is equivalent to running
     *  <code>[estimate_neighborhood_radius( mean_number_of_neighbors(), neighborhood_sample_size() )](\ref estimate_neighborhood_radius)</code>.
     *
     *  This method will be called by the scale-space and surface construction
     *  methods if the neighborhood radius is not set when they are called.
     *
     *  \return the estimated neighborhood radius.
     *
     *  \note This method processes the current point scale-space. The points
     *  in the scale-space can be set with <code>[add_points(begin, end)](\ref add_points)</code>.
     *
     *  \warning If the surface was already constructed, estimating the
     *  neighborhood radius will automatically adjust the surface.
     *
     *  \sa `set_mean_number_of_neighbors(unsigned int neighbors)`.
     *  \sa `set_neighborhood_sample_size(unsigned int samples)`.
     *  \sa `estimate_neighborhood_radius(unsigned int neighbors, unsigned int samples)`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end)`.
     *  \sa `advance_scale_space(unsigned int iterations)`.
     *  \sa `reconstruct_surface(unsigned int iterations)`.
     */
    inline FT estimate_neighborhood_radius() {
        return estimate_neighborhood_radius( mean_number_of_neighbors(), neighborhood_sample_size() );
    }

    /// estimates the neighborhood radius based on a number of sample points.
    /** The neighborhood radius is expressed as the radius of the smallest ball
     *  centered on a point such that the ball contains at least a specified
     *  number of points, not counting the point itself.
     *
     *  The neighborhood radius is set to the mean of these radii, taken over a
     *  number of point samples.
     *
     *  \param neighbors is the number of neighbors a point's neighborhood
     *  should contain on average, not counting the point itself.
     *  \param samples is the number of points sampled to estimate the
     *  neighborhood radius.
     *  \return the estimated neighborhood radius.
     *
     *  \note This method processes the current point scale-space. The points
     *  in the scale-space can be set with <code>[add_points(begin, end)](\ref add_points)</code>.
     *
     *  \warning If the surface was already constructed, estimating the
     *  neighborhood radius will automatically adjust the surface.
     *
     *  \sa `estimate_neighborhood_radius()`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples)`.
     */
	FT estimate_neighborhood_radius( unsigned int neighbors, unsigned int samples );
    
    /// estimates the neighborhood radius of a collection of points.
    /** This method is equivalent to running
     *  `clear()` followed by
     *  <code>[add_points(begin, end)](\ref add_points)</code> and
     *  finally <code>[estimate_neighborhood_radius( mean_number_of_neighbors(), neighborhood_sample_size() )](\ref estimate_neighborhood_radius)</code>.
     *
     *  This method will be called by the scale-space and surface construction
     *  methods if the neighborhood radius is not set when they are called.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The iterator must point to a `Point`.
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \return the estimated neighborhood radius.
     *
     *  \sa `set_mean_number_of_neighbors(unsigned int neighbors)`.
     *  \sa `set_neighborhood_sample_size(unsigned int samples)`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples)`.
     *  \sa `estimate_neighborhood_radius()`.
     *  \sa `add_points(InputIterator begin, InputIterator end)`.
     *  \sa `construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations)`.
     *  \sa `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
    FT estimate_neighborhood_radius( InputIterator begin, InputIterator end ) {
#else // DOXYGEN_RUNNING
    FT estimate_neighborhood_radius( InputIterator begin, InputIterator end,
                     typename boost::enable_if<
                        boost::is_convertible< typename std::iterator_traits<InputIterator>::value_type,
                                               Point > >::type* = NULL ) {
#endif // DOXYGEN_RUNNING
        return estimate_neighborhood_radius( begin, end, mean_number_of_neighbors(), neighborhood_sample_size() );
    }
    
    /// estimates the neighborhood radius of a collection of points based on a number of sample points.
    /** The neighborhood radius is expressed as the radius of the smallest ball
     *  centered on a point such that the ball contains at least a specified
     *  number of points, not counting the point itself.
     *
     *  The neighborhood radius is set to the mean of these radii, taken over a
     *  number of point samples.
     *  
     *  This method is equivalent to running
     *  `clear()` followed by
     *  <code>[add_points(begin, end)](\ref add_points)</code> and finally
     *  <code>[estimate_neighborhood_radius(neighbors, samples)](\ref estimate_neighborhood_radius)</code>.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The iterator must point to a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \param neighbors is the number of neighbors a point's neighborhood
     *  should contain on average, not counting the point itself.
     *  \param samples is the number of points sampled to estimate the
     *  neighborhood radius.
     *  \return the estimated neighborhood radius.
     *
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end)`.
     *  \sa `estimate_neighborhood_radius(unsigned int neighbors, unsigned int samples)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	FT estimate_neighborhood_radius( InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples );
#else // DOXYGEN_RUNNING
	FT estimate_neighborhood_radius( InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples,
                                     typename boost::enable_if<
                                        boost::is_convertible< typename std::iterator_traits<InputIterator>::value_type,
                                                               Point > >::type* = NULL );
#endif // DOXYGEN_RUNNING

/// \}
    
/// \name Scale-space Manipulation
/// \{
    /// advances the scale-space by a number of iterations.
    /** Each iteration the scale-space is advanced, a higher scale-space is
     *  computed. At a higher scale, the scale-space contains a smoother
     *  representation of the point set.
     *
     *  In case the scale-space is not at the scale of the original point set,
     *  calling <code>[advance_scale_space(iterations)](\ref advance_scale_space)</code>
     *  with `iterations > 0` will advance the scale-space further.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_radius()`.
     *
     *  \param iterations is the number of iterations to perform. If
     *  `iterations` is 0, nothing happens.
     *
     *  \note This method processes the current point scale-space. The points
     *  in the scale-space can be set with <code>[add_points(begin, end)](\ref add_points)</code>.
     *
     *  \note If the surface was already constructed, advancing the scale-space
     *  will not automatically adjust the surface.
     *
     *  \sa `construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations)`.
     *  \sa `estimate_neighborhood_radius()`.
     *  \sa `reconstruct_surface(unsigned int iterations)`.
     */
	void advance_scale_space( unsigned int iterations = 1 );
    
    /// constructs a scale-space of a collection of points.
    /** If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_radius()`.
     *
     *  This method is equivalent to running
     *  `clear()` followed by
     *  <code>[add_points(begin, end)](\ref add_points)</code> and finally
     *  <code>[advance_scale_space(iterations)](\ref advance_scale_space)</code>.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The iterator must point to a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \param iterations is the number of iterations to perform. If
     *  `iterations` is 0, nothing happens.
     *
     *  \sa `add_points(InputIterator begin, InputIterator end)`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end)`.
     *  \sa `advance_scale_space(unsigned int iterations)`.
     *  \sa `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	void construct_scale_space( InputIterator begin, InputIterator end, unsigned int iterations = 1 ) {
#else // DOXYGEN_RUNNING
	void construct_scale_space( InputIterator begin, InputIterator end, unsigned int iterations = 1,
                                     typename boost::enable_if<
                                        boost::is_convertible< typename std::iterator_traits<InputIterator>::value_type,
                                                               Point > >::type* = NULL ) {
#endif // DOXYGEN_RUNNING
        clear();
		add_points( begin, end );
		advance_scale_space( iterations );
	}

/// \}
private:
    // constructs the scale-space from a triangulation.
    void construct_scale_space( Triangulation& tr ) {
        add_points( tr.finite_vertices_begin(), tr.finite_vertices_end() );
    }

    // tries to perform a functor in parallel.
    template< class F > void try_parallel( const F& func, size_t begin, size_t end, Sequential_tag ) const;
    template< class F > void try_parallel( const F& func, size_t begin, size_t end, Parallel_tag ) const;
    template< class F > void try_parallel( const F& func, size_t begin, size_t end ) const { try_parallel( func, begin, end, Ct() );}
    
private:
/// \name Shape
/// \{
    /// constructs the shape of the scale-space.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated.
     */
    void construct_shape() {
		construct_shape( scale_space_begin(), scale_space_end() );
    }
    
    /// constructs the shape from an existing triangulation.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  \param tr is the triangulation to construct the shape of.
     *
     *  \note This does not set the current scale-space.
     *  To set this as well, use `construct_scale_space(Triangulation& tr)`.
     *
     *  \note If the neighborhood radius has not been set before, it is automatically
     *  estimated.
     */
    void construct_shape(Triangulation& tr ) {
        deinit_shape();
        if( !has_neighborhood_radius() )
            estimate_neighborhood_radius();
        _shape = Shape_construction_3()( *tr, _squared_radius );
	}

    /// constructs the shape from a collection of points.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  \tparam InputIterator is an iterator over the point sample.
     *  The iterator must point to a `Point`.
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *
     *  \note This does not set the current scale-space.
     *  To set this as well, use `add_points( InputIterator begin, InputIterator end )`.
     *
     *  \note If the neighborhood radius has not been set before, it is automatically
     *  estimated.
     *
     *  \sa `is_constructed()`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	void construct_shape( InputIterator begin, InputIterator end );
#else // DOXYGEN_RUNNING
	void construct_shape( InputIterator begin, InputIterator end,
                                     typename boost::enable_if<
                                        boost::is_convertible< typename std::iterator_traits<InputIterator>::value_type,
                                                               Point > >::type* = NULL );
#endif // DOXYGEN_RUNNING
    
    // collects the surface mesh from the shape.
    // If the sahep does not yet exist, it is constructed.
    void collect_surface();

/// \}

public:
/// \name Surface Reconstruction
/// \{
    /// gives the number of triangles of the surface.
    std::size_t number_of_triangles() const { return _surface.size(); }
    
    /// gives the number of shells of the surface.
    std::size_t number_of_shells() const {
        CGAL_assertion( OS::value == true );
        return _shells.size();
    }
    
    /// constructs a triangle mesh from the scale-space.
    /** The order of the points in the scale-space is the same as the order of
     *  the original points, meaning that the surface of the scale-space can
     *  interpolate the original point set by applying the indices of the
     *  surface triangles to the original point set.
     *
     *  After construction, the triangles of the surface can be iterated using
     *  `surface_begin()` and `surface_end()`.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_radius()`.
     *
     *  \param iterations is the number of scale-space advancement iterations to
     *  apply. If `iterations` is 0, the current scale-space is used.
     *
     *  \note This method processes the current point scale-space. The points
     *  in the scale-space can be set with <code>[add_points(begin, end)](\ref add_points)</code>.
     *
     *  \sa `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     *  \sa `estimate_neighborhood_radius()`.
     *  \sa `advance_scale_space(unsigned int iterations)`.
     */
	void reconstruct_surface( unsigned int iterations = 0 );
    
    /// constructs a surface mesh from the scale-space of a collection of points.
    /** This method is equivalent to running
     *  `clear()` followed by
     *  <code>[add_points(begin, end)](\ref add_points)</code> and finally
     *  <code>[reconstruct_surface(iterations)](\ref reconstruct_surface)</code>.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_radius()`.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The iterator must point to a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \param iterations is the number of scale-space advancement iterations to
     *  apply. If `iterations` is 0, the current scale-space is used.
     *  
     *  \sa `reconstruct_surface(unsigned int iterations)`.
     *  \sa `add_points(InputIterator begin, InputIterator end)`.
     *  \sa `estimate_neighborhood_radius(InputIterator begin, InputIterator end)`.
     *  \sa `construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	void reconstruct_surface( InputIterator begin, InputIterator end, unsigned int iterations = 0 );
#else // DOXYGEN_RUNNING
	void reconstruct_surface( InputIterator begin, InputIterator end, unsigned int iterations = 0,
                                     typename boost::enable_if<
                                        boost::is_convertible< typename std::iterator_traits<InputIterator>::value_type,
                                                               Point > >::type* = NULL );
#endif // DOXYGEN_RUNNING

/// \}

public:
/// \name Iterators
/// \{
    /// gives an iterator to the first point in the current scale space.
    Const_point_iterator scale_space_begin() const { return _tree.begin(); }
    /// gives an iterator to the first point in the current scale space.
    /** \warning Changes to the scale-space do not cause an automatic update to
     *  the surface.
     */
    Point_iterator scale_space_begin() { return _tree.begin(); }

    /// gives a past-the-end iterator of the points in the current scale space.
    Const_point_iterator scale_space_end() const { return _tree.begin(); }
    /// gives a past-the-end iterator of the points in the current scale space.
    /** \warning Changes to the scale-space do not cause an automatic update to
     *  the surface.
     */
    Point_iterator scale_space_end() { return _tree.end(); }

    /// gives an iterator to the first triple in the surface.
    Const_triple_iterator surface_begin() const { return _surface.begin(); }
    /// gives an iterator to the first triple in the surface.
    /** \warning Changes to the surface may change its topology.
     */
    Triple_iterator surface_begin() { return _surface.begin(); }
    
    /// gives a past-the-end iterator of the triples in the surface.
    Const_triple_iterator surface_end() const { return _surface.end(); }
    /// gives a past-the-end iterator of the triples in the surface.
    /** \warning Changes to the surface may change its topology.
     */
    Triple_iterator surface_end() { return _surface.end(); }

    /// gives an iterator to the first triple in a given shell.
    /** \param shell is the index of the shell to access.
     *
     *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
     */
    Const_triple_iterator shell_begin( std::size_t shell ) const;
    /// gives an iterator to the first triple in a given shell.
    /** \param shell is the index of the shell to access.
     *
     *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
     *
     *  \warning Changes to a shell may invalidate the topology of the surface.
     */
    Triple_iterator shell_begin( std::size_t shell );

    /// gives a past-the-end iterator of the triples in a given shell.
    /** \param shell is the index of the shell to access.
     *
     *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
     */
    Const_triple_iterator shell_end( std::size_t shell ) const;

    /// gives a past-the-end iterator of the triples in a given shell.
    /** \param shell is the index of the shell to access.
     *
     *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
     *
     *  \warning Changes to a shell may invalidate the topology of the surface.
     */
    Triple_iterator shell_end( std::size_t shell );

/// \}
}; // class Scale_space_surface_reconstruction_3

} // namespace CGAL

template< typename T >
std::ostream&
operator<<( std::ostream& os, const CGAL::cpp11::array< T, 3 >& t ) {
    return os << get<0>(t) << " " << get<1>(t) << " " << get<2>(t);
}

template< typename T >
std::istream&
operator>>( std::istream& is, CGAL::cpp11::array< T, 3 >& t ) {
    return is >> get<0>(t) >> get<1>(t) >> get<2>(t);
}

#include <CGAL/Scale_space_surface_reconstruction_3_impl.h>

#endif // CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_H
