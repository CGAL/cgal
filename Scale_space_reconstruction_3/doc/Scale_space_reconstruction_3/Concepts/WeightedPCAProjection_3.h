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

    
/// A concept for projecting a point orthogonally onto a weighted least-squares planar approximation of a point set.
/** \ingroup PkgScaleSpaceReconstruction3Concepts
 *  \cgalConcept
 *
 *  The weighted least-squares planar approximation intersects the barycenter
 *  of the points and is orthogonal to the eigenvector corresponding to the
 *  smallest eigenvalue.
 *
 *  \cgalHasModel `CGAL::Weighted_PCA_projection_3`
 */
class WeightedPCAProjection_3 {
public:
/// \name Types
/// \{
	typedef unspecified_type    FT;     ///< defines the number field type.
	typedef unspecified_type    Point;  ///< defines the point type.

/// \}

public:
/// \name Point Insertions.
/// \{
    /// sets a weighted point in the collection.
    /** \pre i must be smaller than the total size of the point set.
     */
    void set_point( unsigned int i, const Point& p, const FT& w );
    
    /// sets the collection of weighted points.
    /** After these points are set, the weighted least-squares planar
     *  approximation is immediately constructed.
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
     *  affine hull of the points is smaller than 2D.
     *
     *  \pre The number of weights must be at least equal to the number of
     *  points.
     */
    template < typename PointIterator, typename WeightIterator >
    bool set_points( PointIterator points_begin, PointIterator points_end,
                     WeightIterator weights_begin );

/// \}

public:
/// \name Weighted PCA Approximation
/// \{
    
    /// computes the weighted least-squares planar approximation of the point set.
    /** \return whether the approximation converges. If the approximation does
     *  not converge this may indicate that the point set is too small, or the
     *  affine hull of the points is smaller than 2D.
     */
    bool approximate();

    /// checks whether the weighted least-squares planar approximation has been computed.
    /** \return `true` iff the approximation has been computed successfully.
     */
    bool is_approximated() const;
/// \}

public:
    /// projects a point onto the weighted least-squares planar approximation.
    /** \param p is the point to project.
     *
     *  \return the orthogonal projection of `p` onto the weighted
     *  least-squares planar approximating of the point set.
     *
     *  \pre The approximating plane must already be computed.
     */
    Point project( const Point& p );
}; // class WeightedPCAProjection_3
