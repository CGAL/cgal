//The surface meshing algorithm implemented using the alpha-shape.
//Copyright (C) 2013 INRIA - Sophia Antipolis
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


#ifndef BALL_PIVOTING
#define BALL_PIVOTING

#include <CGAL/Delaunay_triangulation_3.h>

#define BALL_PIVOTING_FIXED
#define BALL_PIVOTING_CONNECTED

#ifdef BALL_PIVOTING_FIXED
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#else
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#endif

#include "Flagged_triangulation_cell_base_3.h"

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
template < typename Kernel,
		   typename Vb = CGAL::Triangulation_vertex_base_3<Kernel>,
		   typename Cb = CGAL::Triangulation_cell_base_3<Kernel> >
class Surface_mesher {
public:
	typedef typename Kernel::FT											Scalar; ///< The number type.

	// Convert the DT types to be in an alpha-shape.
	/// The vertex type.
#ifdef BALL_PIVOTING_FIXED
	typedef CGAL::Fixed_alpha_shape_vertex_base_3<Kernel, Vb>			AVb;
#else
	typedef CGAL::Alpha_shape_vertex_base_3<Kernel, Vb>					AVb;
#endif

	/// The cell type.
	typedef CGAL::Flagged_triangulation_cell_base_3<Kernel, Cb>			FCb;
#ifdef BALL_PIVOTING_FIXED
	typedef CGAL::Fixed_alpha_shape_cell_base_3<Kernel, FCb>			ACb;
#else
	typedef CGAL::Alpha_shape_cell_base_3<Kernel, FCb>					ACb;
#endif

	// Triangulation.
	typedef CGAL::Triangulation_data_structure_3<AVb, ACb>				Tds;    ///< The data structure that stores the point set.
	typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>					ADT;    ///< The base triangulation that spatially orders the point set.
    /// The structure that identifies the triangles in the surface.
#ifdef BALL_PIVOTING_FIXED
	typedef CGAL::Fixed_alpha_shape_3<ADT>								AS;
#else
	typedef CGAL::Alpha_shape_3<ADT>									AS;
#endif

public:
	typedef ADT															Triangulation;  ///< The triangulation that spatially orders the point set.
	typedef AS															Alpha_shape;    ///< The structure that identifies the triangles in the surface.
	
	typedef typename AS::Vertex_handle									Vertex_handle;  ///< A handle to access the vertices.
	typedef typename AS::Cell_handle									Cell_handle;    ///< A handle to access the cells.
	typedef typename AS::Facet											Facet;          ///< A facet of the triangulation.

	typedef typename Kernel::Point_3									Point;          ///< The point type.
	typedef typename Kernel::Triangle_3									Triangle;       ///< The triangle type.

    /// The constructor that starts the structure from scratch.
    /** \param r2 the squared radius of the ball used to indicate
     *  regions outside the shape.
     *  \sa Surface_mesher(const Triangulation& tr, const Scalar& r2).
     */
	Surface_mesher(const Scalar& r2);

    /// The constructor that uses an existing spatial structure.
    /** \param tr the spatial structure on the point sample.
     *  \param r2 the squared radius of the ball used to indicate
     *  regions outside the shape.
     *  \sa Surface_mesher(const Scalar& r2).
     */
	Surface_mesher(const Triangulation& tr, const Scalar& r2);

    /// Construct a new spatial structure on a set of sample points.
    /** \tparam InputIterator an iterator over the point sample.
     *  The iterator must point to a Point type.
     *  \param start an iterator to the first point of the collection.
     *  \param end a past-the-end iterator for the point collection.
     *  \sa is_constructed().
     */
	template < class InputIterator >
	void construct_shape(InputIterator start, InputIterator end);

    /// Check whether the spatial structure has been constructed.
    /** \return true if the structure exists and false otherwise.
     *  \sa construct_shape(InputIterator start, InputIterator end).
     */
	bool is_constructed() const;

    /// Clear the spatial data structure.
	void clear();

    /// Accessor for the structure indicating the triangles of the mesh.
    /** \return the structure indicating the triangles of the mesh.
     */
	Alpha_shape* shape();
    
    /// Accessor for the ball indicating regions outside the shape.
    /** \return the squared radius of the ball used to indicate
     *  regions outside the shape.
     */
	Scalar get_radius2() const;

    /// Mutator for the ball indicating regions outside the shape.
    /** \param r2 the squared radius of the ball used to indicate
     *  regions outside the shape.
     */
	void set_radius2(const Scalar& r2);

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
	template < class OutputIterator >
	OutputIterator collect_shell(Cell_handle c, unsigned int li, OutputIterator out);
    
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
	template < class OutputIterator >
	OutputIterator collect_shell(const Facet& f, OutputIterator out);

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
	OutputIterator collect_surface(OutputIterator out);
    
    /// Construct the surface trinagles from a set of sample points.
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
	OutputIterator operator()(InputIterator start, InputIterator end, OutputIterator out);

    /// Locate a vertex in the spatial data structure.
    /** \param p the point of the vertex.
     *  \param v the vertex in the data structure.
     *  \param hint where to start looking for the vertex.
     *  Provinding a hint near the vertex will greatly speed
     *  up location.
     *  \return whether the point actually has a vertex in the structure.
     */
	bool locate_vertex(const Point& p, Vertex_handle& v, Cell_handle hint = Cell_handle()) const;

    /// Get the vertices of a facet ordered to point towards the outside of the surface.
    /** This orientation is expressed using the 'right-hand rule'
     *  on the ordered vertices of the facet.
     *  \param f a facet of the data structure as seen from the outside.
     *  \param v0 the first vertex of the facet.
     *  \param v1 the second vertex of the facet.
     *  \param v2 the third vertex of the facet.
     */
	void ordered_vertices(const Facet& f, Vertex_handle& v0, Vertex_handle& v1, Vertex_handle& v2);
    
    /// Get a triangle oriented towards the outside of the surface.
    /** This orientation is expressed using the 'right-hand rule'
     *  on the ordered vertices of the triangle.
     *  \param c a cell touching the triangle from the outside of
     *  the object.
     *  \param li the index of the vertex of the cell opposite to
     *  the triangle touched.
     *  \return the triangle oriented towards the outside of the
     *  surface.
     *  \sa oriented_triangle(const Facet& f) const.
     */
	Triangle oriented_triangle(Cell_handle c, unsigned int li) const;
    
    /// Get a triangle oriented towards the outside of the surface.
    /** This orientation is expressed using the 'right-hand rule'
     *  on the ordered vertices of the triangle.
     *  \param f a facet of the data structure indicating the triangle
     *  as seen from the outside.
     *  \return the triangle oriented towards the outside of the
     *  surface.
     *  \sa oriented_triangle(Cell_handle c, unsigned int li) const.
     */
	Triangle oriented_triangle(const Facet& f) const;
}; // class Surface_mesher

#endif // BALL_PIVOTING
