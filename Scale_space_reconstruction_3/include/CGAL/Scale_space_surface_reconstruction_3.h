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
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <CGAL/utility.h>
#include <CGAL/is_iterator.h>
#include <CGAL/Default.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Random.h>

#include <CGAL/Scale_space_reconstruction_3/Shape_construction_3.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Scale_space_reconstruction_3/Weighted_PCA_approximation_3.h>
#endif // CGAL_EIGEN3_ENABLED

#include <boost/mpl/and.hpp>

namespace CGAL {


  struct Null_function {
    typedef int result_type;
    template <typename T>
    int operator()(const T&) const
    {
      return 0;
    }

  }; 

/// computes a triangulated surface mesh interpolating a point set.
/** \ingroup PkgScaleSpaceReconstruction3Classes
 *  This class stores several of the (intermediate) results. This makes it
 *  easier and more efficient to adjust the parameter settings based on
 *  preliminary results, or to further increase the scale to improve the
 *  results. The class stores the point set at the current scale and the
 *  reconstructed surface, possibly with iterators over the shells.
 *
 *  The class also stores the parameters for estimating the optimal
 *  neighborhood radius and either the lastest estimate or the manually set
 *  radius. This way, the radius can be estimated (again) whenever necessary.
 *  Also note that both increasing the scale and reconstructing the
 *  surface use this radius. By changing or re-estimating the radius between
 *  these operations, they can use separate parameter settings.
 *
 *  The surface can be constructed either for a fixed neighborhood radius, or
 *  for a dynamic radius. When constructing the surface for exactly one
 *  neighborhood radius, it is faster to set `FS` to `Tag_true`. If
 *  the correct neighborhood radius should be changed or estimated multiple
 *  times, it is faster to set `FS` to `Tag_false`.
 *
 *  It is undefined whether a surface with fixed radius may have its radius
 *  changed, but if so, this will likely require more computation time than
 *  changing the radius of a dynamic surface. In either case, it is possible to
 *  change the point set while maintaining the same radius.
 *
 *  The surface can be stored either as an unordered collection of triangles, 
 *  or as a collection ordered by shells. A shell is a maximally connected
 *  component of the surface where connected facets are locally oriented
 *  towards the same side of the surface.
 *  
 *  \tparam Gt is the geometric traits class. It must be a model of
 *  `DelaunayTriangulationTraits_3`. It must have a `RealEmbeddable` field
 *  number type. Generally, `Exact_predicates_inexact_constructions_kernel` is
 *  preferred.
 *  \tparam FS determines whether the surface is expected to be constructed
 *  for a fixed neighborhood radius. It must be a `Boolean_tag` type. The default value is
 *  `Tag_true`. Note that the value of this parameter does not change the result but
 *  only has an impact on the run-time.
 *  \tparam Sh determines whether to collect the surface per shell. It
 *  must be a `Boolean_tag` type. The default value is `Tag_true`.
 *  \tparam wA must be a model of `WeightedPCAProjection_3` and determines how
 *  to approximate a weighted point set. If \ref thirdpartyEigen 3.1.2 (or
 *  greater) is available and CGAL_EIGEN3_ENABLED is defined, then
 *  `Weighted_PCA_approximation_3<DelaunayTriangulationTraits_3>` is used by default.
 *  \tparam Ct indicates whether to use concurrent processing. It must be
 *  either `Sequential_tag` or `Parallel_tag` (the default value).
 */
template < class Gt, class FS = Tag_true, class Sh = Tag_true, class wA = Default, class Ct = Parallel_tag >
class Scale_space_surface_reconstruction_3 {
    typedef typename Default::Get< wA,
#ifdef CGAL_EIGEN3_ENABLED
                                   Weighted_PCA_approximation_3<Gt>
#else // CGAL_EIGEN3_ENABLED
                                   void
#endif // CGAL_EIGEN3_ENABLED
                                 >::type                Approximation;

public:
	typedef typename Gt::Point_3                        Point;          ///< defines the point type.
typedef boost::tuple<Point,int>                           Point_and_int;

private:
    // Searching for neighbors.
    typedef Search_traits_3< Gt >                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
  CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
  Traits_base>                                              Search_traits;
    typedef Orthogonal_k_neighbor_search< Search_traits >
                                                        Static_search;
    typedef Orthogonal_incremental_neighbor_search< Search_traits >
                                                        Dynamic_search;
    typedef typename Dynamic_search::Tree               Search_tree;
    typedef Fuzzy_sphere< Search_traits >               Sphere;
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
     typedef typename Gt::FT                             FT;             ///< defines the field number type.
  //	typedef typename Gt::Point_3                        Point;          ///< defines the point type.


#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                            Point_iterator;         ///< defines an iterator over the points.
    typedef const unspecified_type                      Point_const_iterator;   ///< defines a constant iterator over the points.
#else // DOXYGEN_RUNNING

  // typedef typename Search_tree::iterator              Point_iterator;
  // typedef typename Search_tree::const_iterator        Point_const_iterator;
  typedef typename std::vector<Point>::iterator                          Point_iterator;
  typedef typename std::vector<Point>::const_iterator                          Point_const_iterator;
#endif // DOXYGEN_RUNNING

    typedef CGAL::cpp11::array< unsigned int, 3 >       Triple;                 ///< defines a triple of point indices indicating a triangle of the surface.
private:
    typedef std::list< Triple >                         Tripleset;              ///< defines a collection of triples.
    // Note that this is a list for two reasons: iterator validity for the shell iterators, and memory requirements for the expected huge collections.

public:
#ifdef DOXYGEN_RUNNING
    typedef unspecified_type                            Triple_iterator;        ///< defines an iterator over the triples.
    typedef const unspecified_type                      Triple_const_iterator;  ///< defines a constant iterator over the triples.
#else // DOXYGEN_RUNNING
    typedef Tripleset::iterator                         Triple_iterator;
    typedef Tripleset::const_iterator                   Triple_const_iterator;
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

    typedef std::vector<Point> Pointset;
    Pointset _points;

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
    /** \param sq_radius is the squared radius of the neighborhood.
     *
     *  \pre `sq_radius` is not negative.
     */
	Scale_space_surface_reconstruction_3( FT sq_radius );

/// \}
	~Scale_space_surface_reconstruction_3() { deinit_shape(); }

private:
    void deinit_shape() { if( _shape != 0 ) { delete _shape; _shape = 0; } }

    void clear_tree() { _tree.clear(); }
	void clear_surface() { _shells.clear(); _surface.clear(); deinit_shape(); }
    
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
        if( !has_neighborhood_squared_radius() )
            estimate_neighborhood_squared_radius();
        collect_facets( Sh() );
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
    /// inserts a collection of points into the scale-space at the current scale.
    /** \tparam InputIterator is an iterator over the point collection.
     *  The value type of the iterator must be a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *
     *  \note Inserting the points does not automatically construct or
     *  update the surface.
     *
     *  \note In order to construct the surface, call
     *  `reconstruct_surface()` after inserting the points.
     *
     *  \warning Inserting new points may invalidate the neighborhood radius if
     *  it was previously estimated.
     *
     *  \sa `insert(const Point& p)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	void insert( InputIterator begin, InputIterator end ) {
#else // DOXYGEN_RUNNING
	void insert( InputIterator begin, InputIterator end,
                     typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* = NULL ) {
#endif // DOXYGEN_RUNNING
                Null_function zero;
		_tree.insert(
                             boost::make_zip_iterator(boost::make_tuple( begin, boost::make_transform_iterator(begin, zero))),
                             boost::make_zip_iterator(boost::make_tuple( end , boost::make_transform_iterator(end, zero))));
                _points.insert(_points.end(), begin, end);
	}
    
    /// inserts a point into the scale-space at the current scale.
    /** \param p is the point to insert.
     *
     *  \note Inserting the point does not automatically construct or
     *  update the surface.
     *
     *  \note In order to construct the surface, call
     *  `#reconstruct_surface()`.
     *
     *  \warning Inserting a new point may invalidate the neighborhood radius
     *  if it was previously estimated.
     *
     *  \sa `insert(InputIterator begin, InputIterator end)`.
     */
	void insert( const Point& p ) {
          _tree.insert( boost::make_tuple(p,0) );
                _points.push_back(p);
	}
    
    /// clears the stored scale-space surface reconstruction data.
    /** This includes discarding the surface, the scale-space and all its
     *  points, and any estimation of the neighborhood radius. 
     *
     *  Methods called after this point may have to re-estimate the
     *  neighborhood radius. This method does not discard the parameters for
     *  estimating this radius (the mean number of neighbors and the sample
     *  size).
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
    /// estimates the neighborhood radius.
    /** This method is equivalent to running
     *  <code>[estimate_neighborhood_squared_radius( mean_number_of_neighbors(), neighborhood_sample_size() )](\ref estimate_neighborhood_squared_radius)</code>.
     *
     *  This method can be called by the scale-space and surface construction
     *  methods if the neighborhood radius is not set when they are called.
     *
     *  \return the estimated neighborhood radius.
     *
     *  \note This method processes the point set at the current scale. The
     *  points can be set with <code>[insert(begin, end)](\ref insert)</code>.
     *
     *  \warning If the surface was already constructed, estimating the
     *  neighborhood radius will automatically adjust the surface.
     *
     *  \sa `set_mean_number_of_neighbors(unsigned int neighbors)`.
     *  \sa `set_neighborhood_sample_size(unsigned int samples)`.
     *  \sa `estimate_neighborhood_squared_radius(unsigned int neighbors, unsigned int samples)`.
     *  \sa `increase_scale(unsigned int iterations)`.
     *  \sa `reconstruct_surface()`.
     */
    inline FT estimate_neighborhood_squared_radius() {
        return estimate_neighborhood_squared_radius( mean_number_of_neighbors(), neighborhood_sample_size() );
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
     *  \note This method processes the point set at the current scale. The
     *  points can be set with <code>[insert(begin, end)](\ref insert)</code>.
     *
     *  \warning If the surface was already constructed, estimating the
     *  neighborhood radius will automatically adjust the surface.
     *
     *  \sa `estimate_neighborhood_squared_radius()`.
     */
	FT estimate_neighborhood_squared_radius( unsigned int neighbors, unsigned int samples );


    /// sets the squared radius of the neighborhood.
    /** The neighborhood radius is used by
     *  `#[increase_scale()` to
     *  compute the point set at the desired scale and by
     *  `reconstruct_surface()` to
     *  construct a surface from the point set at the current scale.
     *
     *  \param sq_radius is the squared radius of the neighborhood.
     *
     *  \note If the neighborhood squared radius is negative when the point set
     *  is smoothed or when the surface is computed, the neighborhood radius
     *  will be computed automatically.
     *
     *  \warning If the surface was already constructed, changing the
     *  neighborhood radius will automatically adjust the surface.
     *
     *  \sa `neighborhood_squared_radius()`.
     *  \sa `has_neighborhood_squared_radius()`.
     *  \sa `increase_scale(unsigned int iterations)`.
     *  \sa `reconstruct_surface()`.
     */
    void set_neighborhood_squared_radius( const FT& sq_radius ) {
        _squared_radius = sq_radius;
		if( has_neighborhood_squared_radius() && has_shape() )
            Shape_construction_3().change_scale( _shape, _squared_radius );
    }

    /// gives the squared radius of the neighborhood.
    /** The neighborhood radius is used by
     *  `#increase_scale()` to
     *  compute the point set at the desired scale and by
     *  `#reconstruct_surface()` to
     *  construct a surface from the point set at the current scale.
     *
     *  \return the squared radius of the neighborhood, or -1 if the
     *  neighborhood radius has not yet been set.
     *
     *  \sa `increase_scale(unsigned int iterations)`.
     *  \sa `reconstruct_surface()`.
     */
    FT neighborhood_squared_radius() const { return _squared_radius; }

    /// checks whether the neighborhood radius has been set.
    /** The radius can be set manually, or estimated automatically.
     *
     *  \return `true` iff the radius has been either set manually or estimated.
     *
     *  \sa `set_neighborhood_squared_radius()`.
     *  \sa `estimate_neighborhood_squared_radius()`.
     */
    bool has_neighborhood_squared_radius() const {
        return sign( _squared_radius ) == POSITIVE;
    }
    
    /// \cond internal_doc
    /// estimates the neighborhood radius of a collection of points.
    /** This method is equivalent to running
     *  `clear()` followed by
     *  <code>[insert(begin, end)](\ref insert)</code> and
     *  finally <code>[estimate_neighborhood_squared_radius( mean_number_of_neighbors(), neighborhood_sample_size() )](\ref estimate_neighborhood_squared_radius)</code>.
     *
     *  This method can be called by the scale-space and surface construction
     *  methods if the neighborhood radius is not set when they are called.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The value type of the iterator must be a `Point`.
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \return the estimated neighborhood radius.
     *
     *  \sa `set_mean_number_of_neighbors(unsigned int neighbors)`.
     *  \sa `set_neighborhood_sample_size(unsigned int samples)`.
     *  \sa `estimate_neighborhood_squared_radius()`.
     *  \sa `insert(InputIterator begin, InputIterator end)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
    FT estimate_neighborhood_squared_radius( InputIterator begin, InputIterator end ) {
#else // DOXYGEN_RUNNING
    FT estimate_neighborhood_squared_radius( InputIterator begin, InputIterator end,
                     typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* = NULL ) {
#endif // DOXYGEN_RUNNING
        return estimate_neighborhood_squared_radius( begin, end, mean_number_of_neighbors(), neighborhood_sample_size() );
    }
    /// \endcond

    /// \cond internal_doc
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
     *  <code>[insert(begin, end)](\ref insert)</code> and finally
     *  <code>[estimate_neighborhood_squared_radius(neighbors, samples)](\ref estimate_neighborhood_squared_radius)</code>.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The value type of the iterator must be a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \param neighbors is the number of neighbors a point's neighborhood
     *  should contain on average, not counting the point itself.
     *  \param samples is the number of points sampled to estimate the
     *  neighborhood radius.
     *  \return the estimated neighborhood radius.
     *
     *  \sa `estimate_neighborhood_squared_radius(unsigned int neighbors, unsigned int samples)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	FT estimate_neighborhood_squared_radius( InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples );
#else // DOXYGEN_RUNNING
	FT estimate_neighborhood_squared_radius( InputIterator begin, InputIterator end, unsigned int neighbors, unsigned int samples,
                                         typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* = NULL);
#endif // DOXYGEN_RUNNING
    /// \endcond
/// \}

/// \name Neighborhood Size Estimation Parameters
/// \{
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
     *  \sa `has_neighborhood_squared_radius()`.
     *  \sa `neighborhood_sample_size()`.
     *  \sa `estimate_neighborhood_squared_radius()`.
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
     *  \sa `has_neighborhood_squared_radius()`.
     *  \sa `mean_number_of_neighbors()`.
     *  \sa `estimate_neighborhood_squared_radius()`.
     */
    unsigned int neighborhood_sample_size() const { return _samples; }  

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
     *  \sa `has_neighborhood_squared_radius()`.
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
     *  \sa `has_neighborhood_squared_radius()`.
     *  \sa `set_mean_number_of_neighbors(unsigned int neighbors)`.
     *  \sa `estimate_neighborhood_squared_radius()`.
     */
    void set_neighborhood_sample_size( unsigned int samples ) { _samples = samples; }
/// \}

/// \name Scale-Space Manipulation
/// \{
    /// increases the scale by a number of iterations.
    /** Each iteration the scale is increased, the points set at a higher scale
     *  is computed. At a higher scale, the points set is smoother.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_squared_radius()`.
     *
     *  \param iterations is the number of iterations to perform. If
     *  `iterations` is 0, nothing happens.
     *
     *  \note This method processes the point set at the current scale. The
     *  points can be set with <code>[insert(begin, end)](\ref insert)</code>.
     *
     *  \note If the surface was already constructed, increasing the scale
     *  will not automatically adjust the surface.
     *
     *  \sa `estimate_neighborhood_squared_radius()`.
     *  \sa `reconstruct_surface()`.
     */
	void increase_scale( unsigned int iterations = 1 );

    /// \cond internal_doc
    /// constructs a scale-space of a collection of points.
    /** If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_squared_radius()`.
     *
     *  This method is equivalent to running
     *  `clear()` followed by
     *  <code>[insert(begin, end)](\ref insert)</code> and finally
     *  <code>[increase_scale(iterations)](\ref increase_scale)</code>.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The value type of the iterator must be a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \param iterations is the number of iterations to perform. If
     *  `iterations` is 0, nothing happens.
     *
     *  \sa `insert(InputIterator begin, InputIterator end)`.
     *  \sa `estimate_neighborhood_squared_radius(InputIterator begin, InputIterator end)`.
     *  \sa `increase_scale(unsigned int iterations)`.
     *  \sa `reconstruct_surface(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
	template < class InputIterator >
	void construct_scale_space( InputIterator begin, InputIterator end, unsigned int iterations = 1,
                                    typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* = NULL ) {
        clear();
		insert( begin, end );
		increase_scale( iterations );
	}
    ///\endcond

/// \}
private:
    // constructs the scale-space from a triangulation.
    void construct_scale_space( Triangulation& tr ) {
        insert( tr.finite_vertices_begin(), tr.finite_vertices_end() );
    }

    // tries to perform a functor in parallel.
    template< class F > void try_parallel( const F& func, std::size_t begin, std::size_t end, Sequential_tag ) const;
    template< class F > void try_parallel( const F& func, std::size_t begin, std::size_t end, Parallel_tag ) const;
    template< class F > void try_parallel( const F& func, std::size_t begin, std::size_t end ) const { try_parallel( func, begin, end, Ct() );}
    
private:
/// \name Shape
/// \{
    /// constructs the shape of the points at a fixed scale.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated.
     */
    void construct_shape() {
		construct_shape( points_begin(), points_end() );
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
        if( !has_neighborhood_squared_radius() )
            estimate_neighborhood_squared_radius();
        _shape = Shape_construction_3()( *tr, _squared_radius );
	}

    /// constructs the shape from a collection of points.
    /** The shape contains geometric and connectivity information
     *  of the scale space.
     *
     *  \tparam InputIterator is an iterator over the point sample.
     *  The value type of the iterator must be a `Point`.
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *
     *  \note This does not set the current scale-space.
     *  To set this as well, use `insert( InputIterator begin, InputIterator end )`.
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
                              typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* = NULL );
#endif // DOXYGEN_RUNNING
    
    // collects the surface mesh from the shape.
    // If the sahep does not yet exist, it is constructed.
    void collect_surface();

/// \}

public:
/// \name Surface Reconstruction
/// \{
    /// constructs a triangle mesh from the point set at a fixed scale.
    /** The order of the points at the current scale is the same as the order
     *  at the original scale, meaning that the surface can interpolate the
     *  point set at the original scale by applying the indices of the surface
     *  triangles to the original point set.
     *
     *  After construction, the triangles of the surface can be iterated using
     *  `surface_begin()` and `surface_end()`.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_squared_radius()`.
     *
     *  \note This method processes the point set at the current scale. The
     *  points can be set with <code>[insert(begin, end)](\ref insert)</code>.
     *
     *  \sa `estimate_neighborhood_squared_radius()`.
     *  \sa `increase_scale(unsigned int iterations)`.
     */
        void reconstruct_surface()
        {
          reconstruct_surface(0);
        }

    /// gives the number of triangles of the surface.
    std::size_t number_of_triangles() const { return _surface.size(); }
    
    /// gives the number of points of the surface.
    std::size_t number_of_points() const { return _tree.size(); }


    /// gives the number of shells of the surface.
    std::size_t number_of_shells() const {
        CGAL_assertion( Sh::value == true );
        return _shells.size();
    }

    /// \cond internal_doc
        void reconstruct_surface( unsigned int iterations);
    /// \endcond

    /// \cond internal_doc
    /// constructs a surface mesh from a collection of points at a fixed scale.
    /** This method is equivalent to running
     *  `clear()` followed by
     *  <code>[insert(begin, end)](\ref insert)</code> and finally
     *  <code>[reconstruct_surface(iterations)](\ref reconstruct_surface)</code>.
     *
     *  If the neighborhood radius has not been set before, it is automatically
     *  estimated using `estimate_neighborhood_squared_radius()`.
     *
     *  \tparam InputIterator is an iterator over the point collection.
     *  The value type of the iterator must be a `Point`.
     *
     *  \param begin is an iterator to the first point of the collection.
     *  \param end is a past-the-end iterator for the point collection.
     *  \param iterations is the number of scale increase iterations to apply.
     *  If `iterations` is 0, the point set at the current scale is used.
     *  
     *  \sa `reconstruct_surface(unsigned int iterations)`.
     *  \sa `insert(InputIterator begin, InputIterator end)`.
     *  \sa `estimate_neighborhood_squared_radius(InputIterator begin, InputIterator end)`.
     *  \sa `construct_scale_space(InputIterator begin, InputIterator end, unsigned int iterations)`.
     */
	template < class InputIterator >
#ifdef DOXYGEN_RUNNING
	void reconstruct_surface( InputIterator begin, InputIterator end, unsigned int iterations = 0 );
#else // DOXYGEN_RUNNING
	void reconstruct_surface( InputIterator begin, InputIterator end, unsigned int iterations = 0,
                                  typename boost::enable_if< CGAL::is_iterator<InputIterator> >::type* = NULL );
#endif // DOXYGEN_RUNNING
    /// \endcond
/// \}

public:
/// \name Iterators
/// \{
    /// gives an iterator to the first point at the current scale.
    Point_const_iterator points_begin() const { return _points.begin(); }
    /// gives an iterator to the first point at the current scale.
    /** \warning Changes to the scale-space do not cause an automatic update to
     *  the surface.
     */
    Point_iterator points_begin() { return _points.begin(); }

    /// gives a past-the-end iterator of the points at the current scale.
    Point_const_iterator points_end() const { return _points.end(); }
    /// gives a past-the-end iterator of the points at the current scale.
    /** \warning Changes to the scale-space do not cause an automatic update to
     *  the surface.
     */
    Point_iterator points_end() { return _points.end(); }

    /// gives an iterator to the first triple in the surface.
    Triple_const_iterator surface_begin() const { return _surface.begin(); }
    /// gives an iterator to the first triple in the surface.
    /** \warning Changes to the surface may change its topology.
     */
    Triple_iterator surface_begin() { return _surface.begin(); }
    
    /// gives a past-the-end iterator of the triples in the surface.
    Triple_const_iterator surface_end() const { return _surface.end(); }
    /// gives a past-the-end iterator of the triples in the surface.
    /** \warning Changes to the surface may change its topology.
     */
    Triple_iterator surface_end() { return _surface.end(); }

    /// gives an iterator to the first triple in a given shell.
    /** \param shell is the index of the shell to access.
     *
     *  \pre `shell` is in the range [ 0, `number_of_shells()` ).
     */
    Triple_const_iterator shell_begin( std::size_t shell ) const;
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
    Triple_const_iterator shell_end( std::size_t shell ) const;

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
    return os << t[0] << " " << t[1] << " " << t[2];
}

template< typename T >
std::istream&
operator>>( std::istream& is, CGAL::cpp11::array< T, 3 >& t ) {
    return is >> get<0>(t) >> get<1>(t) >> get<2>(t);
}

#include <CGAL/Scale_space_reconstruction_3/Scale_space_surface_reconstruction_3_impl.h>

#endif // CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_H
