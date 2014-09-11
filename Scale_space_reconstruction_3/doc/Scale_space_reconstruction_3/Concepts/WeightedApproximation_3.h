//Orthogonal projection of a point onto a weighted PCA plane.
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

    
/// A concept for computing a approximation of a weighted point set.
/** \ingroup PkgScaleSpaceReconstruction3Concepts
 *  \cgalConcept
 *  This approximation can be used to fit other points to the point set.
 *
 *  \cgalHasModel `CGAL::Weighted_PCA_approximation_3`
 */
class WeightedApproximation_3 {
public:
/// \name Types
/// \{
	typedef unspecified_type    FT;     ///< defines the field number type.
	typedef unspecified_type    Point;  ///< defines the point type.

/// \}
    
public:
/// \name Constructors
/// \{
    /// constructs an approximation of a undefined point set.
    /** The point set holds a fixed number of points with undefined coordinates.
     *
     *  \param size is the size of the points set.
     *
     *  \note this does not compute the approximation.
     */
    WeightedApproximation_3( unsigned int size );
    
/// \}

    //  computes the approximation of a point set.
    /*  Similar to constructing a approximation and calling
     *  <code>[set_points(points_begin, points_end, weights_begin)](\ref WeightedApproximation_3::set_points )</code>
     *
     *  \tparam PointIterator is an input iterator over the point collection.
     *  The value type of the iterator must be a `Point`.
     *  \tparam WeightIterator is an input iterator over the collection of
     *  weights. The value type of the iterator must be an `FT`.
     *
     *  \param points_begin is an iterator to the first point.
     *  \param points_end is a past-the-end oterator for the points.
     *  \param weights_begin is an iterator to the weight of the first point.
     *
     *  \pre The number of weights must be at least equal to the number of
     *  points.
     */
    template < typename PointIterator, typename WeightIterator >
    WeightedApproximation_3( PointIterator points_begin, PointIterator points_end, WeightIterator weights_begin );

public:
/// \name Point Set.
/// \{
    /// changes a weighted point in the set.
    /** This invalidates the approximation. `compute()` should be called
     *  after all points have been set.
     *  \pre i must be smaller than the total size of the point set.
     */
    void set_point( unsigned int i, const Point& p, const FT& w );
    
    /// gives the size of the weighted point set.
    std::size_t size() const;

/// \}

    //  changes the weighted point set.
    /*  After these points are set, the approximation is immediately computed.
     *
     *  \tparam PointIterator is an input iterator over the point collection.
     *  The value type of the iterator must be a `Point`.
     *  \tparam WeightIterator is an input iterator over the collection of
     *  weights. The value type of the iterator must be an `FT`.
     *
     *  \param points_begin is an iterator to the first point.
     *  \param points_end is a past-the-end iterator for the points.
     *  \param weights_begin is an iterator to the weight of the first point.
     *
     *  \return whether the approximation converges. If the approximation does
     *  not converge this may indicate that the point set is too small, or the
     *  affine hull of the points cannot contain the approximation.
     *
     *  \pre The points must fit in the approximation.
     *  \pre The number of weights must be at least equal to the number of
     *  points.
     */
    template < typename PointIterator, typename WeightIterator >
    bool set_points( PointIterator points_begin, PointIterator points_end,
                     WeightIterator weights_begin );

public:
/// \name Approximation
/// \{
    
    /// computes the approximation.
    /** \return whether the approximation converges. If the approximation does
     *  not converge this may indicate that the point set is too small, or the
     *  affine hull of the points cannot contain the approximation.
     */
    bool compute();

    /// checks whether the approximation has been computed successfully.
    bool is_computed() const;
/// \}

public:
/// \name Fitting
/// \{
    
    /// fits a point to the approximation.
    /** \param p is the point to fit.
     *
     *  \return the point on the approximation closest to `p`.
     *
     *  \pre The approximation must have been computed.
     */
    Point fit( const Point& p );
/// \}
}; // class WeightedApproximation_3
