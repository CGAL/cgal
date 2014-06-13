//A class that follows a straight line through a Delaunay triangulation structure.
//Copyright (C) 2012 - Utrecht University
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
// Author(s): Thijs van Lankveld

#ifndef CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
#define CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H

#include <CGAL/internal/Dummy_gt.h>
#include <CGAL/Random.h>

namespace CGAL {

/** 
 *	\ingroup PkgTriangulation3UtilityClasses
 *
 *	The `Triangulation_segment_traverser_3` iterates over the cells of a
 *  `Triangulation_3` by following a straight line segment
 *  \f$ \overline{st} \f$.
 *
 *	The functionality of this class is closely related to
 *  `Triangulation_3::locate(...)`. 
 *	However, unlike this `locate(...)` method, all the cells traversed by the
 *  `Triangulation_segment_traverser_3` intersect the interior of the line
 *  segment \f$ \overline{st} \f$.
 *
 *	Traversal starts from a cell containing \f$ s \f$ and it ends in a cell
 *  containing \f$ t \f$.
 *	If \f$ \overline{st} \f$ is coplanar with a facet or collinear with an
 *  edge, at most one of the incident cells is traversed.
 *	If \f$ \overline{st} \f$ intersects an edge or vertex (and
 *  \f$ \overline{st} \f$ isn't coplanar with any of the incident facets), at
 *  most two incident cells are traversed: the cells intersecting
 *  \f$ \overline{st} \f$ strictly in their interior.
 *
 *	If \f$ s \f$ lies on the convex hull, traversal starts in an incident cell
 *  inside the convex hull. Similarly, if \f$ t \f$ lies on the convex hull,
 *  traversal ends in an incident cell inside the convex hull.
 *
 *	Both \f$ s \f$ and \f$ t \f$ may lie outside the convex hull of the
 *  triangulation, but they must lie within the affine hull of the
 *  triangulation. In either case, the finite facet of any infinite cells
 *  traversed must intersect \f$ \overline{st} \f$.
 *
 *	The traverser may be applied to triangulations of dimension smaller than 3,
 *  as long as the dimension is larger than 0.
 *	However, for triangulations of dimension 1, the functionality is somewhat
 *  trivial.
 *
 *	The traverser becomes invalid whenever the triangulation is changed.
 *
 *	\tparam Tr is the type of the triangulation to traverse. It must be a
 *  `Triangulation_3` or one of its subclasses.
 *
 *	\sa `Triangulation_3`
 *	\sa `Forward_circulator_base`
 */
template < class Tr >
class Triangulation_segment_traverser_3
: public Forward_circulator_base< typename Tr::Triangulation_data_structure::Cell,
                                  std::ptrdiff_t, std::size_t > {
public:
/// \name Types
/// \{
    typedef Tr                                          Triangulation;  ///< defines the triangulation type.
    typedef Triangulation_segment_traverser_3<Tr>       TST;            ///< defines the segment traverser type.
    typedef typename Tr::Point                          Point;          ///< defines the point type.
    typedef typename Tr::Segment                        Segment;        ///< defines the segment type.
    typedef typename Tr::Cell                           Cell;           ///< defines the cell type.
    typedef typename Tr::Edge                           Edge;           ///< defines the edge type.
    typedef typename Tr::Vertex_handle                  Vertex_handle;  ///< defines the handle type for a vertex.
    typedef typename Tr::Cell_handle                    Cell_handle;    ///< defines the handle type for a cell.
    typedef typename Tr::Locate_type                    Locate_type;    ///< defines the simplex types returned when locating a point.

/// \}
    
protected:
/// \internal \name Protected Attributes
/// \{
    const Tr&       _tr; ///< \internal stores the triangulation to traverse.

/// \}

public:
/// \name Constructors
/// \{
    /// constructs a traverser.
    /** \param tr the triangulation to traverse. 
     *	\param s the source vertex.
     *	\param t the target vertex.
     *  \pre both `s` and `t` must be initialized and cannot be the infinite
     *  vertex.
     *	\pre `s != t`.
     *  \pre `tr` must have dimension > 0.
     */
	Triangulation_segment_traverser_3( const Tr& tr, Vertex_handle s, Vertex_handle t );

    /// constructs a traverser.
    /** \param tr the triangulation to traverse. 
     *	\param s the source vertex.
     *	\param t the target point.
     *  \pre both `s` and `t` must be initialized and `s` cannot be the
     *  infinite vertex.
     *	\pre `s != t`.
     *  \pre `tr` must have dimension > 0. If `tr` has dimension < 3, `t` must
     * lie in the affine hull of `tr`.
     */
	Triangulation_segment_traverser_3( const Tr& tr, Vertex_handle s, const Point& t );

    /// constructs a traverser.
    /** \param tr the triangulation to traverse. 
     *	\param s the source point.
     *	\param t the target point.
     *	\param hint the starting point to locate `s`.
     *  \pre both `s` and `t` must be initialized.
     *	\pre `s != t`.
     *  \pre `tr` must have dimension > 0. If `tr` has dimension < 3, both `s`
     *  and `t` must lie in the affine hull of `tr`.
     */
    Triangulation_segment_traverser_3( const Tr& tr, const Point& s, const Point& t, Cell_handle hint = Cell_handle() );
    
    /// constructs a traverser.
    /** \param tr the triangulation to traverse. 
     *	\param seg the line segment to traverse.
     *	\param hint the starting point to locate the source of `seg`.
     *  \pre `seg` must not be degenerate, i.e. its source and target must not
     *  be equal.
     *  \pre `tr` must have dimension > 0. If `tr` has dimension < 3, `seg`
     *  must lie in the affine hull of `tr`.
     */
    Triangulation_segment_traverser_3( const Tr& tr, const Segment& seg, Cell_handle hint = Cell_handle() );

/// \}

public:
/// \name Accessors
/// \{
    /// gives the source point of the segment being traversed.
    const Point&    source() const      { return _source; }

    /// gives the target point of the segment being traversed.
    const Point&    target() const      { return _target; }

    /// gives the cell currently traversed.
    /** By invariance, this cell is intersected by the segment between
     *  `source()` and `target()`. However, it need not be intersected in its
     *  iterior.
	 *	\return the cell currently traversed.
	 *	\sa `handle()`.
	 */
    const Cell      cell() const        { return *_pos; }

    /// gives a handle to the cell currently traversed.
    /** By invariance, this cell is intersected by the segment between
     *  `source()` and `target()`. However, it need not be intersected in its
     *  iterior.
	 *	\return a handle to the cell currently traversed.
	 *	\sa `cell()`.
	 */
    Cell_handle     handle()            { return _pos; }

    /// gives the cell traversed previously.
    /** This cell is uninitialized until the traverser leaves the initial
	 *	cell containing the `source()`.
	 *	By invariance, once initialized, this cell must be intersected by the
     *  segment between `source()` and `target()`.
	 *	\return the cell traversed previously.
	 *	\sa `handle()`.
	 */
    Cell_handle     previous() const    { return _prev; }

    /// is the dereference operator.
    /**	\return a pointer to the cell currently traversed.
	 */
    Cell*           operator->()        { return &*_pos; }

    /// is the indirection operator.
    /** \return the cell currently traversed.
	 */
    Cell&           operator*()         { return *_pos; }

    /// is the conversion operator.
    operator const  Cell_handle() const { return _pos; }

    /// checks if the traverser has reached the final cell, which contains the `target()`.
    /** If the `target()` lies on a facet, edge, or vertex, the final cell is
     *  the cell containing the interior of the segment between `source()` and
     *  `target()`.
	 *	\return true iff the current cell contains the `target()`.
	 */
    bool            has_next() const    { return !_done; }

    /// indicates how the current cell was entered.
    /** For the first cell, containing the `source()` \f$ s \f$, this indicates
     *  the location of \f$ s \f$ in this cell.
	 */
    void            entry( Locate_type& lt, int& li, int& lj ) const    { lt = _lt; li = _li; lj = _lj; }

/// \}

public:
/// \name Mutators
/// \{
    /// is the increment postfix operator.
    /**	After incrementing the traverser, the current cell intersects the
     *  segment between `source()` and `target()`, incident to the previous
     *  cell, and closer to the `target()` than the previous cell.
     *	\sa `operator++(int)`.
     *  \pre The current cell does not contain the `target()`.
     */
    TST&            operator++();
		
    /// is the increment prefix operator.
    /**	After incrementing the traverser, the current cell intersects the
     *  segment between `source()` and `target()`, incident to the previous
     *  cell, and closer to the `target()` than the previous cell.
     *	\sa `operator++()`.
     *  \pre The current cell does not contain the `target()`.
     */
    TST             operator++( int );

    /// traverses to the final cell, which contains the `target()`.
    /** This is equivalent to calling `operator++()` until `has_next()` is
     *  `false`.
     *  \return the final cell.
     */
    Cell_handle     traverse();

/// \}
	
public:
/// \name Comparison
/// \{
    /// compares this traverser with `ct`.
    /** \param ct the other traverser.
     *	\return `true` iff the other traverser traverses the same triangulation
     *  along the same line segment and has the same current cell.
     *	\sa `operator!=( const TST& t )`.
     */
    bool            operator==( const TST& ct ) const;
		
    /// compares this traverser with `ct`.
    /** \param ct the other traverser.
     *	\return `false` iff the other traverser traverses the same
     *  triangulation along the same line segment and has the same current cell.
     *	\sa `operator==( const TST& t ) const`.
     */
    bool            operator!=( const TST& ct ) const;

    /// compares the cell currently traversed with `ch`.
    /** \param ch a handle to the cell to compare to.
     *	\return true iff the cell currently traversed is the same as the one
     *  pointed to by `ch`. 
     *	\sa `operator!=( const Cell_handle& ch ) const`.
     *	\sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
     */
    bool            operator==( const Cell_handle& ch ) const   { return ch == _pos; }
		
    /// compares the cell currently traversed with `ch`.
    /** \param ch a handle to the cell to compare to.
     *	\return `false` iff the cell currently traversed is the same as the one
     *  pointed to by `ch`.
     *	\sa `operator==( const Cell_handle& ch )`.
     *	\sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
     */
    bool            operator!=( const Cell_handle& ch ) const   { return ch != _pos; }

/// \}

protected:
/// \internal \name Protected Member Functions
/// \{
    /// \internal traverses the triangulation to the next cell.
    /** \internal \sa `traverse()`.
	 */
    virtual void    increment();

/// \}
}; // class Triangulation_segment_traverser_3
	
/// compares a handle to a cell to a traverser.
/** \param ch the handle to a cell.
 *	\param t the traverser.
 *	\return true iff the cell currently traversed by `t` is the same as the one
 *  pointed to by `ch`.
 *	\sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
 *	\sa `Triangulation_segment_traverser_3::operator==( const Cell_handle& ch )`.
 */
template <class Tr>
inline bool operator==( typename Tr::Cell_handle ch, Triangulation_segment_traverser_3<Tr> tci ) { return tci == ch; }

/// compares a handle to a cell to a traverser.
/** \param ch the handle to a cell.
 *	\param t the traverser.
 *	\return `false` iff the cell currently traversed by `t` is the same as the
 *  one pointed to by `ch`.
 *	\sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
 *	\sa `Triangulation_segment_traverser_3::operator!=( const Cell_handle& ch )`.
 */
template <class Tr>
inline bool operator!=( typename Tr::Cell_handle ch, Triangulation_segment_traverser_3<Tr> tci ) { return tci != ch; }

} // namespace CGAL

#include <CGAL/Triangulation_segment_traverser_3_impl.h>

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
