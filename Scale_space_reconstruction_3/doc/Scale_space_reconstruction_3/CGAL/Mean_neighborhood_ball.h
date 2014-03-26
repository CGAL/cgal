//A method that computes the size of a mean neighborhood.
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


#ifndef MEAN_NEIGHBORHOOD_BALL
#define MEAN_NEIGHBORHOOD_BALL

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Random.h>

/// Compute the mean neighborhood size.
/** Each neighborhood size is expressed as the radius of the  smallest
 *  ball centered on a sample point such that the ball contains at least
 *  a specified number of points.
 *
 *  The mean neighborhood size is then the mean of these radii, taken over
 *  a number of point samples.
 *  \tparam Kernel the geometric traits class. This class specifies, amongst
 *  other things, the number types and predicates used.
 */
template < typename Kernel >
class Mean_neighborhood_ball {
	typedef typename Kernel::FT											Scalar;             ///< The number type.

	typedef typename Kernel::Point_3									Point;              ///< The point type.

	typedef Scalar														result_type;        ///< The type for expressing the mean neigborhood size.

public:
    /// The constructor.
    /** \param nn the number of nearest neighbors that the neigborhood should contain.
     *  \param smp the number of smaple points used to estimate the neighborhood size.
     */
	Mean_neighborhood_ball(unsigned int nn = 30, unsigned int smp = 200);

    /// Accessor for the number of neighbors required in the neighborhood.
    /** \return the number of neighbors the neighborhood must contain.
     */
    unsigned int get_neighbors() const;

    /// Accessor for the number of sample points used to estimate the neighborhood size.
    /** \return the number of sample points used to estimate the neighborhood size.
     */
    unsigned int get_samples() const;

    /// Mutator for the number of neighbors required in the neighborhood.
    /** \param nn the number of neighbors required in the neighborhood.
     */
    void set_neighbors( unsigned int nn );

    /// Mutator for the number of sample points used to estimate the neighborhood size.
    /** \param smp the number of sample points used to estimate the neighborhood size.
     */
    void set_samples( unsigned int smp );
	
    /// Construct a search tree from a collection of points.
    /** This tree must be constructed before handling sample points.
     *  \tparam InputIterator an iterator over a collection of points.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa compute().
     *  \sa operator()(InputIterator start, InputIterator end).
     */
	template < class InputIterator >
	void constructTree(InputIterator start, InputIterator end);

    /// Incorporate the neighborhood of a sample point into the estimate.
    /** Note that a total of get_samples() points should be handled.
     *  \param p the point to incorporate in the estimate.
     *  \sa compute().
     */
	void handlePoint(const Point& p);

    /// Estimate the mean neighborhood size based on a number of sample points.
    /** It is assumed that the seach tree has already been constructed.
     *  \return the estimated mean neighborhood size.
     *  \sa handlePoint(const Point& p).
     *  \sa operator()(InputIterator start, InputIterator end).
     */
	result_type compute();

    /// Estimate the mean neighborhood size of a collection of points.
    /** This method is equivalent to running [constructTree(start, end)](\ref constructTree)
     *  followed by compute().
     *  \tparam InputIterator an iterator over a collection of points.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa constructTree(InputIterator start, InputIterator end).
     *  \sa compute().
     */
	template < class InputIterator >
	result_type operator()(InputIterator start, InputIterator end);
}; // class Mean_neighborhood_ball

#endif MEAN_NEIGHBORHOOD_BALL
