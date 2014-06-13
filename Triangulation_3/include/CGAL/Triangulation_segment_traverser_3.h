//A class that follows a straight line through a Delaunay triangulation structure.
//Copyright (C) 2012  Utrecht University
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

#ifndef DOXYGEN_RUNNING
/** 
 *	\ingroup PkgConformingTriangulation3Classes
 *
 *	The `Triangulation_segment_traverser_3` iterates over the cells
 *	of a `Triangulation_3` by following a straight line segment \f$ st \f$.
 *
 *	This class is closely related to `Triangulation_3::locate(...)`.
 *	However, unlike this `locate(...)` method, all the cells traversed
 *	by the `Triangulation_segment_traverser_3` intersect the interior of the line
 *	segment \f$ st \f$.
 *
 *	Traversal starts from a cell containing \f$ s \f$ and it ends in a cell containing
 *	\f$ t \f$.
 *	If \f$ st \f$ is coplanar with a facet or collinear with an edge, at most one of the
 *	incident cells is traversed.
 *	If \f$ st \f$ intersects an edge or vertex, at most two incident cells are traversed:
 *	the cells intersecting \f$ st \f$ strictly in their interior.
 *
 *	If \f$ s \f$ lies on the convex hull, traversal starts in an incident cell inside
 *	the convex hull. Similarly, if \f$ t \f$ lies on the convex hull, traversal ends in
 *	an adjacent cell inside the convex hull.
 *
 *	Both \f$ s \f$ and \f$ t \f$ may lie outside the convex hull of the triangulation,
 *	but they must lie within the affine hull of the triangulation. In either case, the
 *	finite facet of any infinite cells traversed must intersect \f$ st \f$.
 *
 *	The traverser may be applied to any triangulation of dimension > 0.
 *	However, for triangulations of dimension 1, the functionality is somewhat trivial.
 *
 *	The traverser becomes invalid whenever the triangulation is changed.
 *
 *	\tparam Tr_ is the triangulation type to traverse.
 *
 *	\sa `Triangulation_3`
 *	\sa `Forward_circulator_base`
 */
#endif // DOXYGEN_RUNNING
template < class Tr_ = void >
class Triangulation_segment_traverser_3
: public Forward_circulator_base< typename Tr_::Triangulation_data_structure::Cell,
                                  std::ptrdiff_t, std::size_t > {
    typedef	Tr_                                         Tr;
    typedef typename Tr::Triangulation_data_structure   Tds;
    typedef typename Tr::Geom_traits                    Gt;

public:
/// \name Types
/// \{
    typedef Tr                                          Triangulation;  ///< The triangulation type.
    typedef Triangulation_segment_traverser_3<Tr>       TST;            ///< A triangulation segment traverser.
    typedef typename Tr::Point                          Point;          ///< A point embedded in the vertices of the triangulation.
    typedef typename Tr::Segment                        Segment;        ///< A segment connecting two points.
    typedef typename Tr::Cell                           Cell;           ///< A cell in the triangulation.
    typedef typename Tr::Edge                           Edge;           ///< An edge in the triangulation.
    typedef typename Tr::Vertex_handle                  Vertex_handle;  ///< A handle for a vertex in the triangulation.
    typedef typename Tr::Cell_handle                    Cell_handle;    ///< A handle for a cell in the triangulation.
    typedef typename Tr::Locate_type                    Locate_type;    ///< The simplex containing a geometric object.
/// \}
    
    template < class Tr2 >
    struct Rebind_Tr {
        typedef Triangulation_segment_traverser_3<Tr2>  Other;
    };

protected:
/// \internal \name Protected Attributes
/// \{
    /// \internal The triangulation to traverse.
    const Tr&       _tr;
/// \}

    // The source and target points of the traversal.
    Point           _source;
    Point           _target;
    Vertex_handle   _s_vertex;
    Vertex_handle   _t_vertex;

    // The cell currently traversed and the previous one.
    Cell_handle	_pos, _prev;

    // Otherwise, they indicate the simplex through which this cell was entered,
    // or the location of the source if it is in this cell.
    int             _li, _lj;
    Locate_type     _lt;

    // This bit signifies when a cell containing the target is found.
    bool            _done;

    // Where possible, facets are checked in a random order.
    mutable Random rng;

private:
    Tds             _tds2;
    Vertex_handle   _s_vert, _t_vert;

public:
/// \name Constructors
/// \{
    /// Constructs a traverser.
    /** \param tr the triangulation to traverse. This triangulation must have dimension > 0.
     *	\param s the source vertex. This vertex must be initialized and cannot be the infinite vertex.
     *	\param t the target vertex. This vertex must be initialized and cannot be the infinite vertex.
     *	It cannot equal `s`.
     */
	Triangulation_segment_traverser_3( const Tr& tr, Vertex_handle s, Vertex_handle t );

    /// Constructs a traverser.
    /** \param tr the triangulation to traverse. This triangulation must have dimension > 0.
     *	\param s the source vertex. This vertex must be initialized and cannot be the infinite vertex.
     *	\param t the target point. This point must be initialized and it cannot be be at the same location as `s`.
     *	If `tr` has dimension < 3, `t` must lie inside the affine hull of `tr`.
     */
	Triangulation_segment_traverser_3( const Tr& tr, Vertex_handle s, const Point& t );

    /// Constructs a traverser.
    /** \param tr the triangulation to traverse. This triangulation must have dimension > 0.
     *	\param s the source point. This point must be initialized. If `tr` has dimension < 3, `s` must lie inside
     *	the affine hull of `tr`.
     *	\param t the target point. This point must be initialized and it cannot be be at the same location as `s`.
     *	If `tr` has dimension < 3, `t` must lie inside the affine hull of `tr`.
     *	\param hint the starting point to search for `s`.
     */
    Triangulation_segment_traverser_3( const Tr& tr, const Point& s, const Point& t, Cell_handle hint = Cell_handle() );
    
    /// Constructs a traverser.
    /** \param tr the triangulation to traverse. This triangulation must have dimension > 0.
     *	\param S the segment to be traversed. If `tr` has dimension < 3, `S` must lie inside
     *	the affine hull of `tr`. `S` must not be degenerate, i.e. its source and target must not be equal.
     *	\param hint the starting point to search for `S`.
     */
    Triangulation_segment_traverser_3( const Tr& tr, const Segment& S, Cell_handle hint = Cell_handle() );
/// \}

public:
/// \name Accessors
/// \{
    /// Gives the source point of the traversed segment.
    /** \return the source point.
     */
    const Point&    source() const      { return _source; }

    /// Gives the target point of the traversed segment.
    /** \return the target point.
	 */
    const Point&    target() const      { return _target; }

    /// Gives the cell currently traversed.
    /** By invariance, this cell is intersected by the segment
	 *	between `source()` and `target()`.
	 *	\return the cell currently traversed.
	 *	\sa `handle()`.
	 */
    Cell_handle     cell() const        { return _pos; }

    /// Gives a handle to the cell currently traversed.
    /** By invariance, this cell is intersected by the segment
	 *	between `source()` and `target()`.
	 *	\return a handle to the cell currently traversed.
	 *	\sa `cell()`.
	 */
    Cell_handle     handle()            { return _pos; }

    /// Gives the cell traversed previously.
    /** This cell is uninitialized until the traverser leaves the initial
	 *	cell containing the `source()`.
	 *	By invariance, once initialized, this cell must be intersected by the segment
	 *	between `source()` and `target()`.
	 *	\return the cell traversed previously.
	 *	\sa `handle()`.
	 */
    Cell_handle     previous() const    { return _prev; }

    /// Dereference operator.
    /**	\return a pointer to the cell currently traversed.
	 */
    Cell*           operator->()        { return &*_pos; }

    /// Indirection operator.
    /** \return the cell currently traversed.
	 */
    Cell&           operator*()         { return *_pos; }

    /// Conversion operator.
    /**	\return a handle to the cell currently traversed.
	 */
    operator const  Cell_handle() const { return _pos; }

    /// Checks if the traverser has reached the final cell, which contains the `target()`.
    /** If the `target()` lies on a facet, edge, or vertex, the final cell is the cell containing
	 *	the interior of the segment between `source()` and `target()`.
	 *	\return true iff the current cell contains the `target()`.
	 */
    bool            has_next() const    { return !_done; }

    /// Gives the type of simplex traversed last.
    /**	This simplex indicates where traversal has entered the cell.
	 *	For the first cell, containing the `source()` \f$ s \f$, this indicates the location of \f$ s \f$ in this cell.
	 */
    void            entry( Locate_type& lt, int& li, int& lj ) const    { lt = _lt; li = _li; lj = _lj; }
/// \}

public:
/// \name Mutators
/// \{
    /// Increment postfix operator.
    /**	The current cell must not contain the `target()`.
     *	After incrementing the traverser, the current cell intersects the segment
     *	between `source()` and `target()` closer to the `target()` than the previous cell.
     *	\sa `operator++(int)`.
     */
    TST&            operator++();
		
    /// Increment prefix operator.
    /**	The current cell must not contain the `target()`.
     *	After incrementing the traverser, the current cell intersects the segment
     *	between `source()` and `target()` closer to the `target()` than the previous cell.
     *	than the previous cell.
     *	\sa `operator++()`.
     */
    TST             operator++( int );

    /// Traverse to the final cell, which contains the `target()`.
    /**	\return the final cell.
     */
    Cell_handle     traverse();
/// \}
	
public:
/// \name Comparison
/// \{
    /// Compares this traverser with `ct`.
    /** \param ct the other traverser.
     *	\return true iff the other traverser traverses the same triangulation along the same line segment
     *	and has the same current cell.
     *	\sa `operator!=( const TST& t )`.
     */
    bool            operator==( const TST& ct ) const;
		
    /// Compares this traverser with `ct`.
    /** \param ct the other traverser.
     *	\return `false` iff the other traverser traverses the same triangulation along the same line segment
     *	and has the same current cell.
     *	\sa `operator==( const TST& t ) const`.
     */
    bool            operator!=( const TST& ct ) const;

    /// Compares the cell currently traversed with `ch`.
    /** \param ch a handle to the other cell.
     *	\return true iff the cell currently traversed is the same as the one pointed to by `ch`.
     *	\sa `operator!=( const Cell_handle& ch ) const`.
     *	\sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
     */
    bool            operator==( const Cell_handle& ch ) const   { return ch == _pos; }
		
    /// Compares the cell currently traversed with `ch`.
    /** \param ch a handle to the other cell.
     *	\return `false` iff the cell currently traversed is the same as the one pointed to by `ch`.
     *	\sa `operator==( const Cell_handle& ch )`.
     *	\sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
     */
    bool            operator!=( const Cell_handle& ch ) const   { return ch != _pos; }
/// \}

	bool            operator==( Nullptr_t CGAL_triangulation_assertion_code(n) ) const;
	bool            operator!=( Nullptr_t n ) const;

protected:
/// \internal \name Protected Member Functions
/// \{
    /// \internal Traverse to the next cell.
    /** \internal \sa `traverse()`.
	 */
    virtual void    increment();
/// \}

private:
    // increment(), if the triangulation is 3D.
    void            increment_3();
    void            increment_3_inf( int inf );
	
    // increment(), if the triangulation is 2D.
    void            increment_2();
    void            increment_2_inf( int inf );

private:
    inline int      edgeIndex( int i, int j ) const {
        CGAL_triangulation_precondition( i>=0 && i<=3 );
        CGAL_triangulation_precondition( j>=0 && j<=3 );
        CGAL_triangulation_precondition( i != j );
        return ( i==0 || j==0 ) ? i+j-1 : i+j;
    }
}; // class Triangulation_segment_traverser_3
	
/// Compares a handle to a cell to a traverser.
/** \param ch the handle to a cell.
 *	\param t the traverser.
 *	\return true iff the cell currently traversed by `t` is the same as the one pointed to by `ch`.
 *	\sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
 *	\sa `Triangulation_segment_traverser_3::operator==( const Cell_handle& ch )`.
 */
template <class Tr>
inline bool operator==( typename Tr::Cell_handle ch, Triangulation_segment_traverser_3<Tr> tci ) { return tci == ch; }

/// Compares a handle to a cell to a traverser.
/** \param ch the handle to a cell.
 *	\param t the traverser.
 *	\return `false` iff the cell currently traversed by `t` is the same as the one pointed to by `ch`.
 *	\sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_traverser_3<TriangulationTraits_3> t )`.
 *	\sa `Triangulation_segment_traverser_3::operator!=( const Cell_handle& ch )`.
 */
template <class Tr>
inline bool operator!=( typename Tr::Cell_handle ch, Triangulation_segment_traverser_3<Tr> tci ) { return tci != ch; }

// Specialization for void
template <>
class Triangulation_segment_traverser_3<void> {
    typedef internal::Dummy_gt        Gt;
    typedef internal::Dummy_tds_3     Tds;
public:
    typedef Triangulation_3<Gt,Tds> Triangulation;
    template <typename Tr2>
    struct Rebind_Tr { typedef Triangulation_segment_traverser_3<Tr2> Other; };
}; // class Triangulation_segment_traverser_3<void>

} // namespace CGAL

#include <CGAL/Triangulation_segment_traverser_3_impl.h>

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
