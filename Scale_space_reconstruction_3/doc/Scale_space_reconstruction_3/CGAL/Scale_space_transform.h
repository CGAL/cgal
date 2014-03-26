//A method to transform a point set.
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


#ifndef SCALE_SPACE_PERTURB
#define SCALE_SPACE_PERTURB

#include <omp.h>

#include <map>
#include <vector>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>

#include <Eigen/Dense>

#include "check3264.h"


/// Transform a point set to a rougher scale.
/** This transfor is achieved by projecting each point onto
 *  the best planar fit of the weighted local neighborhood.
 *  \tparam Kernel the geometric traits class. This class specifies, amongst
 *  other things, the number types and predicates used.
 */
template < typename Kernel >
class Scale_space_transform {
	typedef typename Kernel::FT													Scalar;     ///< The number type.

	typedef typename Kernel::Point_3											Point;	    ///< The point type.
	typedef typename std::vector<Point>											Collection; ///< A collection of points.

public:
    /// The constructor.
    /** \param r the radius of the local smoothing neighborhood.
     */
	Scale_space_transform(const Scalar& r);

    /// Accessor for the radius of the local smoothing neighborhood.
    /** \return the radius of the local smoothing neighborhood.
     */
	Scalar get_radius2() const;
    
    /// Mutator for the radius of the local smoothing neighborhood.
    /** \param r the radius of the local smoothing neighborhood.
     */
	void set_radius2(const Scalar& r);

    /// Assign a collection of points to be smoothed.
    /** \tparam InputIterator an iterator over a collection of points.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa operator()(InputIterator start, InputIterator end).
     *  \sa operator()(InputIterator start, InputIterator end, unsigned int iterations).
     */
	template < class InputIterator >
	void assign(InputIterator start, InputIterator end);

    /// Accessor for the collection of points to be smoothed.
    /** \return the collection of points to be smoothed.
     */
	Collection& get_collection();
    
    /// Accessor for the collection of points to be smoothed.
    /** \return the collection of points to be smoothed.
     */
	const Collection& get_collection() const;

    /// Compute one iteration of scale-space transforming.
    /** The result of this iteration is stored in the point collection.
     *  If earlier iterations have been computed, calling iterate()
     *  will add the next iteration.
     *  \sa iterate(unsigned int iterations).
     *  \sa get_collection().
     *  \sa get_collection() const.
     */
	void iterate();
    
    /// Compute a number of iterations of scale-space transforming.
    /** The result of this iteration is stored in the point collection.
     *  \param iterations the number of iterations to perform.
     *  \sa iterate().
     *  \sa get_collection().
     *  \sa get_collection() const.
     */
	void iterate(unsigned int iterations);
	
    /// Compute one transform iteration on a collection of points.
    /** This method is equivalent to running [assign(start, end)](\ref assign)
     *  followed by iterate().
     *  \tparam InputIterator an iterator over a collection of points.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa operator()(InputIterator start, InputIterator end, unsigned int iterations).
     *  \sa assign(InputIterator start, InputIterator end).
     *  \sa iterate().
     */
	template < class InputIterator >
	void operator()(InputIterator start, InputIterator end);
    
    /// Compute a number of transform iterations on a collection of points.
    /** This method is equivalent to running [assign(start, end)](\ref assign)
     *  followed by [iterate(iterations)](\ref iterate).
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
	void operator()(InputIterator start, InputIterator end, unsigned int iterations) {
		// Perturb the points for a number of iterations.
		_points.assign(start, end);
		iterate(iterations);
	}
}; // class Scale_space_transform

#endif // SCALE_SPACE_PERTURB