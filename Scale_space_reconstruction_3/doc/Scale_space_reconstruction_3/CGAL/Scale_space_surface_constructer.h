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


#ifndef SCALE_SPACE_SURFACE_CONSTRUCTER
#define SCALE_SPACE_SURFACE_CONSTRUCTER

#include <iostream>
#include <list>
#include <map>

#include "Mean_neighborhood_ball.h"
#include "Scale_space_transform.h"
#include "Surface_mesher.h"

/// Compute a smoothed surface mesh from a collection of points.
/** This is a convenience class for using Mean_neighborhood_ball,
 *  Scale_space_transform, and Surface_mesher.
 *
 *  An appropriate neighborhood size is estimated, followed by a
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
 *  \tparam Vb the vertex base used in the internal data
 *  structure that spatially orders the point set.
 *  \tparam Cb the cell base used in the internal data
 *  structure that spatially orders the point set.
 *  \sa Mean_neighborhood_ball.
 *  \sa Scale_space_transform.
 *  \sa Surface_mesher.
 */
template < typename Kernel,
		   typename Vb = CGAL::Triangulation_vertex_base_3<Kernel>,
		   typename Cb = CGAL::Triangulation_cell_base_3<Kernel> >
class Scale_space_surface_constructer {
public:
	typedef typename Kernel::FT								Scalar;             ///< The number type.

	typedef typename Kernel::Point_3						Point;              ///< The point type.
	typedef typename Kernel::Triangle_3						Triangle;           ///< The triangle type.
	
	typedef typename Surface_mesher::Vertex_handle			Vertex_handle;      ///< A handle to access the vertices.
	typedef typename Surface_mesher::Facet					Facet;              ///< A facet of the triangulation.

	typedef typename std::vector<Point>						PointCollection;    ///< A collection of points.
    
public:
    /// The constructor.
    /** \param iterations the number of smoothing iterations to perform.
     */
	Scale_space_surface_constructer(unsigned int iterations = 1);
    
    /// Mutator for the number of smoothing iterations.
    /** \param iterations the number of smoothing iterations to perform.
     */
	void set_iterations(unsigned int iterations);

    /// Accessor for the smoothed points.
    /** \return the current collection of smoothed points.
     *  Note that this collection may change as smoothing
     *  iterations are performed.
     */
	const PointCollection& moved() const;

	/// Compute a smoothed surface mesh from a collection of points.
    /** \tparam InputIterator an iterator over the point collection.
     *  The iterator must point to a Point type.
     *  \tparam OutputIterator an output iterator for the surface
     *  triangles. The iterator must point to a Triangle type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \param out an iterator to place to output the triangles.
     */
	template < class InputIterator, class OutputIterator >
	OutputIterator operator()(InputIterator start, InputIterator end, OutputIterator out);
    
    /// Get the vertices of a facet ordered to point towards the outside of the surface.
    /** This orientation is expressed using the 'right-hand rule'
     *  on the ordered vertices of the facet.
     *  \param f a facet of the data structure as seen from the outside.
     *  \param v0 the first vertex of the facet.
     *  \param v1 the second vertex of the facet.
     *  \param v2 the third vertex of the facet.
     */
	void ordered_vertices(const Facet& f, Vertex_handle& v0, Vertex_handle& v1, Vertex_handle& v2);
}; // class Scale_space_surface_constructer

#endif // SCALE_SPACE_SURFACE_CONSTRUCTER