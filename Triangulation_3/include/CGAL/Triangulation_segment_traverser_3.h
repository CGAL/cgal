// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Thijs van Lankveld, Jane Tournois


#ifndef CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
#define CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H

#include <CGAL/license/Triangulation_3.h>

#include <iostream>
#include <utility>

#include <CGAL/assertions.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_simplex_3.h>

#include <boost/optional.hpp>

// If defined, type casting is done statically,
// reducing type-safety overhead.
#define CGAL_TST_ASSUME_CORRECT_TYPES

namespace CGAL {

template < class Tr, class Inc >
class Triangulation_segment_cell_iterator_3;

namespace internal {
template < class Tr >
struct Incrementer {
    typedef Incrementer<Tr>                                 Self;
    typedef Triangulation_segment_cell_iterator_3<Tr,Self>  SCI;    // describes the type of iterator expected by the incrementer.
    Incrementer() {}
    void increment( SCI& sci ) { sci.walk_to_next(); }
}; // struct Incrementer

} // namespace internal


//  provides an iterator over the cells intersected by a line segment.
/*
 *        The `Triangulation_segment_traverser_3` iterates over the cells
 *        of a `Triangulation_3` by following a straight line segment \f$ st \f$.
 *
 *        This class is closely related to `Triangulation_3::locate(...)`.
 *        However, unlike this `locate(...)` method, all the cells traversed
 *        by the `Triangulation_segment_traverser_3` intersect the interior of the line
 *        segment \f$ st \f$.
 *
 *        Traversal starts from a cell containing \f$ s \f$ and it ends in a cell containing
 *        \f$ t \f$.
 *        If \f$ st \f$ is coplanar with a facet or collinear with an edge, at most one of the
 *        incident cells is traversed.
 *        If \f$ st \f$ intersects an edge or vertex, at most two incident cells are traversed:
 *        the cells intersected by \f$ st \f$ strictly in their interior.
 *
 *        If \f$ s \f$ lies on the convex hull, traversal starts in an incident cell inside
 *        the convex hull. Similarly, if \f$ t \f$ lies on the convex hull, traversal ends in
 *        an adjacent cell inside the convex hull.
 *
 *        Both \f$ s \f$ and \f$ t \f$ may lie outside the convex hull of the triangulation,
 *        but they must lie within the affine hull of the triangulation. In either case, the
 *        finite facet of any infinite cells traversed must intersect \f$ st \f$.
 *
 *        The traverser may be applied to any triangulation of dimension > 0.
 *        However, for triangulations of dimension 1, the functionality is somewhat trivial.
 *
 *        The traverser becomes invalid whenever the triangulation is changed.
 *
 *        \tparam Tr_ is the triangulation type to traverse.
 *
 *        \cgalModels{ForwardIterator}
 *
 *        \sa `Triangulation_3`
 *        \sa `Forward_circulator_base`
 */
template < class Tr_, class Inc = internal::Incrementer<Tr_> >
class Triangulation_segment_cell_iterator_3
{
    typedef Tr_                                         Tr;
    typedef typename Tr::Triangulation_data_structure   Tds;
    typedef typename Tr::Geom_traits                    Gt;
    typedef Inc                                         Incrementer;

public:
// \name Types
// \{
    typedef Tr                                          Triangulation;          //< defines the triangulation type.
    typedef Triangulation_segment_cell_iterator_3<Tr,Inc>
                                                        Segment_cell_iterator;  //< defines the segment cell iterator type.

    typedef typename Tr::Point                          Point;                  //< defines the point type.
    typedef typename Tr::Segment                        Segment;                //< defines the line segment type.

    typedef typename Tr::Cell                           Cell;                   //< defines the type of a cell of the triangulation.
    typedef typename Tr::Edge                           Edge;                   //< defines the type of an edge of the triangulation.
    typedef typename Tr::Facet                          Facet;                  //< defines the type of a facet of the triangulation.

    typedef typename Tr::Vertex_handle                  Vertex_handle;          //< defines the type of a handle for a vertex in the triangulation.
    typedef typename Tr::Cell_handle                    Cell_handle;            //< defines the type of a handle for a cell in the triangulation.

    typedef typename Tr::Locate_type                    Locate_type;            //< defines the simplex type returned from location.

    struct Simplex                                                              //< defines the simplex type
    {
      Cell_handle cell = {};
      Locate_type lt = Locate_type::OUTSIDE_AFFINE_HULL;
      int li = -1;
      int lj = -1;
    };

    typedef Cell                                        value_type;             //< defines the value type the iterator refers to.
    typedef Cell&                                       reference;              //< defines the reference type of the iterator.
    typedef Cell*                                       pointer;                //< defines the pointer type of the iterator.
    typedef std::size_t                                 size_type;              //< defines the integral type that can hold the size of a sequence.
    typedef std::ptrdiff_t                              difference_type;        //< defines the signed integral type that can hold the distance between two iterators.
    typedef std::forward_iterator_tag                   iterator_category;      //< defines the iterator category.
// \}

    // describes the iterator type when applied to another type of triangulation or incrementer.
    template < class Tr2, class Inc2 >
    struct Rebind { typedef Triangulation_segment_cell_iterator_3<Tr2,Inc2>  Other; };

#if CGAL_DEBUG_TRIANGULATION_SEGMENT_TRAVERSER_3
    static auto display_vert(Vertex_handle v)
    {
      std::stringstream os;
      os.precision(17);
      if(v->time_stamp() == 0) {
        os << "inf";
      } else {
        os << '#' << v->time_stamp() << "=(" << v->point() << ")";
      }
      return os.str();
    };

    static auto display_lt(Locate_type lt) {
      std::stringstream os;
      switch(lt) {
        case Locate_type::VERTEX: os << " VERTEX"; break;
        case Locate_type::EDGE: os << " EDGE"; break;
        case Locate_type::FACET: os << " FACET"; break;
        case Locate_type::CELL: os << " CELL"; break;
        case Locate_type::OUTSIDE_CONVEX_HULL: os << " OUTSIDE_CONVEX_HULL"; break;
        case Locate_type::OUTSIDE_AFFINE_HULL: os << " OUTSIDE_AFFINE_HULL"; break;
      }
      return os.str();
    }

    static auto debug_simplex(Simplex s) {
      std::stringstream os;
      os.precision(17);
      const auto [c, lt, i, j] = s;
      if(c == Cell_handle{}) {
        os << "end()";
      } else {
        os << display_vert(c->vertex(0)) << " - " << display_vert(c->vertex(1)) << " - "
           << display_vert(c->vertex(2)) << " - " << display_vert(c->vertex(3));
        os << display_lt(lt) << " " << i << " " << j;
      }
      return os.str();
    }

    auto debug_iterator() const
    {
      std::stringstream os;
      os.precision(17);
      os << "  prev: " << debug_simplex(_prev) << "\n  cur: " << debug_simplex(_cur);
      return os.str();
    }
#endif // CGAL_DEBUG_TRIANGULATION_SEGMENT_TRAVERSER_3

private:
    typedef Segment_cell_iterator                       SCI;

    friend internal::Incrementer<Tr>;

protected:
// \internal \name Protected Attributes
// \{
    // \internal The triangulation to traverse.
    const Tr*       _tr;

// \}

    // The source and target points of the traversal.
    // These are also stored as vertices for cheaper equality computation.
    Point           _source;
    Point           _target;
    Vertex_handle   _s_vertex;
    Vertex_handle   _t_vertex;

    // The current cell with its entry point and the previous cell with its
    // exit point.
    // Note that the current cell will be Cell_handle() after incrementing past
    // the first cell containing the target.
    Simplex _cur, _prev;

public:
// \name Constructors
// \{
    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *        \param s the source vertex. This vertex must be initialized and cannot be the infinite vertex.
     *        \param t the target vertex. This vertex must be initialized and cannot be the infinite vertex.
     *        It cannot equal `s`.
     */
        Triangulation_segment_cell_iterator_3( const Tr* tr, Vertex_handle s, Vertex_handle t );

    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *        \param s the source vertex. This vertex must be initialized and cannot be the infinite vertex.
     *        \param t the target point. This point must be initialized and it cannot be be at the same location as `s`.
     *        If `tr` has dimension < 3, `t` must lie inside the affine hull of `tr`.
     */
        Triangulation_segment_cell_iterator_3( const Tr* tr, Vertex_handle s, const Point& t );

    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *        \param s the source point. This point must be initialized and it cannot be be at the same location as `t`.
     *        \param t the target vertex. This vertex must be initialized and cannot be the infinite vertex.
     *        If `tr` has dimension < 3, `s` must lie inside the affine hull of `tr`.
     *        \param hint the starting point to search for `s`.
     */
        Triangulation_segment_cell_iterator_3( const Tr* tr, const Point& s, Vertex_handle t, Cell_handle hint = Cell_handle() );

    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *        \param s the source point. This point must be initialized. If `tr` has dimension < 3, `s` must lie inside
     *        the affine hull of `tr`.
     *        \param t the target point. This point must be initialized and it cannot be be at the same location as `s`.
     *        If `tr` has dimension < 3, `t` must lie inside the affine hull of `tr`.
     *        \param hint the starting point to search for `s`.
     */
    Triangulation_segment_cell_iterator_3( const Tr* tr, const Point& s, const Point& t, Cell_handle hint = Cell_handle() );

    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *        \param S the segment to be traversed. If `tr` has dimension < 3, `S` must lie inside
     *        the affine hull of `tr`. `S` must not be degenerate, i.e. its source and target must not be equal.
     *        \param hint the starting point to search for `S`.
     */
    Triangulation_segment_cell_iterator_3( const Tr* tr, const Segment& S, Cell_handle hint = Cell_handle() );
// \}

    // private constructor that does not initialize the source and target.
    // used for the end()
    Triangulation_segment_cell_iterator_3(const Tr* tr);

#ifndef CGAL_TST_ASSUME_CORRECT_TYPES
    // The virtual destructor is mainly defined to indicate to the casting
    // operators that this is a dynamic type.
    virtual
#endif
    ~Triangulation_segment_cell_iterator_3() {}


public:
// \name Accessors
// \{

    const Tr* triangulation() const     { return _tr; }

    //  gives the source point of the segment followed.
    /*  \return the source point.
     */
    const Point&    source() const      { return _source; }

    //  gives the target point of the segment followed.
    /*  \return the target point.
         */
    const Point&    target() const      { return _target; }

    Vertex_handle target_vertex() const { return _t_vertex; }

    //  gives a handle to the current cell.
    /*  By invariance, this cell is intersected by the segment
         *        between `source()` and `target()`.
         *        \return a handle to the current cell.
         *        \sa `cell()`.
         */
    Cell_handle     handle() const
    {
      return _cur.cell;
    }

    //  gives the previous cell.
    /*  This cell is uninitialized until the iterator leaves the initial
         *        cell.
         *        By invariance, once initialized, this cell must be intersected by the segment
         *        between `source()` and `target()`.
         *        \return the previous cell.
         *        \sa `handle()`.
         */
    Cell_handle     previous() const
    {
      return prev_cell();
    }

    //  provides a dereference operator.
    /*         \return a pointer to the current cell.
         */
    Cell*           operator->()
    {
      return &*(_cur.cell);
    }

    //  provides an indirection operator.
    /*  \return the current cell.
         */
    Cell&           operator*()
    {
      return *(_cur.cell);
    }

    //  provides a conversion operator.
    /*         \return a handle to the current cell.
         */
    operator Cell_handle() const
    {
      return _cur.cell;
    }

    //  provides a conversion operator.
    /*         \return the simplex through which the current cell was entered.
         */
    operator Simplex() const { return _cur; }

    //  checks whether the iterator has reached the final cell, which contains the `target()`.
    /*  If the `target()` lies on a facet, edge, or vertex, the final cell is the cell containing
         *        the interior of the segment between `source()` and `target()`.
         *        \return true iff the current cell contains the `target()`.
         */
    bool            has_next() const
    {
      return this->cell() != Cell_handle();
    }

    //  gives the simplex through which the current cell was entered.
    /*         For the first cell, containing the `source()` \f$ s \f$,
     *  this indicates the location of \f$ s \f$ in this cell.
         */
    void            entry( Locate_type& lt, int& li, int& lj ) const
    {
      lt = this->lt(); li = this->li(); lj = this->lj();
    }
    std::tuple<Locate_type, int, int> entry() const
    {
      return { lt(), li(), lj() };
    }
    //  gives the simplex through which the previous cell was exited.
    /*         \pre the current cell is not the initial cell.
         */
    void            exit( Locate_type& lt, int& li, int& lj ) const
    {
      lt = prev_lt(); li = prev_li(); lj = prev_lj();
    }
    std::tuple<Locate_type, int, int> exit() const
    {
      return { prev_lt(), prev_li(), prev_lj() };
    }

    //  gives the past-the-end iterator associated with this iterator.
    SCI             end() const;
// \}

public:
// \name Mutators
// \{
    //  provides the increment postfix operator.
    /*         After incrementing the iterator, the current cell intersects the segment
     *        between `source()` and `target()` closer to the `target()` than the previous cell.
     *        \sa `operator++(int)`.
     *  \pre The current cell does not contain the `target()`.
     */
    SCI&            operator++();

    //  provides the increment prefix operator.
    /*         After incrementing the iterator, the current cell intersects the segment
     *        between `source()` and `target()` closer to the `target()` than the previous cell.
     *        than the previous cell.
     *        \sa `operator++()`.
     *  \pre The current cell does not contain the `target()`.
     */
    SCI             operator++( int );

    //  iterates to the final cell, which contains the `target()`.
    /*         \return the final cell.
     */
    Cell_handle     complete();
// \}

public:
// \name Comparison
// \{
    //  compares this iterator with `sci`.
    /*  \param sci the other iterator.
     *        \return true iff the other iterator iterates the same triangulation along the same line segment
     *        and has the same current cell.
     *        \sa `operator!=( const SCI& t )`.
     */
    bool            operator==( const SCI& sci ) const;

    //  compares this iterator with `sci`.
    /*  \param sci the other iterator.
     *        \return `false` iff the other iterator iterates the same triangulation along the same line segment
     *        and has the same current cell.
     *        \sa `operator==( const SCI& t ) const`.
     */
    bool            operator!=( const SCI& sci ) const;

    //  compares the current cell with `ch`.
    /*  \param ch a handle to the other cell.
     *        \return true iff the current cell is the same as the one pointed to by `ch`.
     *        \sa `operator!=( const Cell_handle& ch ) const`.
     *        \sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
     */
    bool            operator==( const Cell_handle& ch ) const
    {
      return ch == _cur.cell;
    }

    //  compares the current cell with `ch`.
    /*  \param ch a handle to the other cell.
     *        \return `false` iff the current cell is the same as the one pointed to by `ch`.
     *        \sa `operator==( const Cell_handle& ch )`.
     *        \sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
     */
    bool            operator!=( const Cell_handle& ch ) const
    {
      return ch != _cur.cell;
    }
// \}

        bool            operator==( Nullptr_t CGAL_assertion_code(n) ) const;
        bool            operator!=( Nullptr_t n ) const;

protected:
// \internal \name Protected Member Functions
// \{
    //  walks to the next cell.
    /*  \sa `complete()`.
         */
    void            walk_to_next();

    //  increments the iterator.
    /*  This method may perform more actions based on the superclass.
     *  \sa `complete()`.
         */
    void            increment() {
        typedef typename Incrementer::SCI    Expected;
#ifdef CGAL_TST_ASSUME_CORRECT_TYPES
        Expected& sci = static_cast<Expected&>( *this );
#else // CGAL_TST_ASSUME_CORRECT_TYPES
        Expected& sci = dynamic_cast<Expected&>( *this );
#endif // CGAL_TST_ASSUME_CORRECT_TYPES
        Incrementer().increment( sci );
    }
// \}

private:
    // at the end of the constructors, entry() is a vertex, edge or facet,
    // we need to circulate/iterate over its incident cells to
    // make sure that the current cell intersects the input query
    void jump_to_intersecting_cell();

    //  walk_to_next(), if the triangulation is 3D.
    std::pair<Simplex, Simplex> walk_to_next_3(const Simplex& prev,
                                               const Simplex& cur) const;
    void            walk_to_next_3_inf( int inf );

    //  walk_to_next(), if the triangulation is 2D.
    void            walk_to_next_2();
    void            walk_to_next_2_inf( int inf );

private:
    inline int      edgeIndex( int i, int j ) const {
        CGAL_precondition( i>=0 && i<=3 );
        CGAL_precondition( j>=0 && j<=3 );
        CGAL_precondition( i != j );
        return ( i==0 || j==0 ) ? i+j-1 : i+j;
    }

    bool have_same_entry(const Simplex& s1, const Simplex& s2) const;

    // Compute the orientation of a point compared to the oriented plane supporting a half-facet.
    CGAL::Orientation orientation(const Facet& f, const Point& p) const;

    bool coplanar(const Facet &f, const Point &p) const;

    // Gives the edge incident to the same cell that is not incident to any of the input vertices.
    Edge opposite_edge(Cell_handle c, int li, int lj) const;
    Edge opposite_edge(const Edge& e) const;

protected:
    // ref-accessors to the simplex, for use in internal code
    // access _cur
    Cell_handle& cell()             { return _cur.cell; }
    Cell_handle const& cell() const { return _cur.cell; }

    Locate_type& lt()             { return _cur.lt; }
    Locate_type const& lt() const { return _cur.lt; }

    int& li()             { return _cur.li; }
    int const& li() const { return _cur.li; }

    int& lj()             { return _cur.lj; }
    int const& lj() const { return _cur.lj; }

    // access _prev
    Cell_handle& prev_cell()             { return _prev.cell; }
    Cell_handle const& prev_cell() const { return _prev.cell; }

    Locate_type& prev_lt()             { return _prev.lt; }
    Locate_type const& prev_lt() const { return _prev.lt; }

    int& prev_li()             { return _prev.li; }
    int const& prev_li() const { return _prev.li; }

    int& prev_lj()             { return _prev.lj; }
    int const& prev_lj() const { return _prev.lj; }

}; // class Triangulation_segment_cell_iterator_3

//  compares a handle to a cell to a traverser.
/*  \param ch the handle to a cell.
 *        \param t the traverser.
 *        \return true iff the cell currently traversed by `t` is the same as the one pointed to by `ch`.
 *        \sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
 *        \sa `Triangulation_segment_cell_iterator_3::operator==( const Cell_handle& ch )`.
 */
template < class Tr, class Inc >
inline bool operator==( typename Tr::Cell_handle ch, Triangulation_segment_cell_iterator_3<Tr,Inc> tci ) { return tci == ch; }

//  compares a handle to a cell to a traverser.
/*  \param ch the handle to a cell.
 *        \param t the traverser.
 *        \return `false` iff the cell currently traversed by `t` is the same as the one pointed to by `ch`.
 *        \sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
 *        \sa `Triangulation_segment_cell_iterator_3::operator!=( const Cell_handle& ch )`.
 */
template < class Tr, class Inc >
inline bool operator!=( typename Tr::Cell_handle ch, Triangulation_segment_cell_iterator_3<Tr,Inc> tci ) { return tci != ch; }



/********************************************************************/
/********************************************************************/
/********************************************************************/
template < class Tr_, class Inc = internal::Incrementer<Tr_> >
class Triangulation_segment_simplex_iterator_3
{
  typedef Tr_                                         Tr;
  typedef typename Tr::Triangulation_data_structure   Tds;
  typedef typename Tr::Geom_traits                    Gt;
  typedef Inc                                         Incrementer;

private:
  typedef Triangulation_segment_simplex_iterator_3<Tr_, Inc> Simplex_iterator;
  typedef Triangulation_segment_cell_iterator_3<Tr_, Inc>    SCI;

private:
  typedef typename SCI::Point    Point;
  typedef typename SCI::Segment  Segment;

public:
  // \{
  typedef typename SCI::Vertex_handle Vertex_handle;//< defines the type of a handle for a vertex in the triangulation
  typedef typename SCI::Cell_handle   Cell_handle;  //< defines the type of a handle for a cell in the triangulation.
  typedef typename SCI::Cell          Cell;  //< defines the type of a handle for a cell in the triangulation.
  typedef typename SCI::Triangulation::Edge  Edge;  //< defines the type of an edge in the triangulation.
  typedef typename SCI::Triangulation::Facet Facet; //< defines the type of a facet in the triangulation.
  typedef typename SCI::Locate_type   Locate_type;  //< defines the simplex type returned from location.

  typedef CGAL::Triangulation_simplex_3<Tds>  Simplex_3;

  typedef Simplex_3   value_type;       //< defines the value type  the iterator refers to.
  typedef const Simplex_3&  reference;  //< defines the reference type of the iterator.
  typedef const Simplex_3*  pointer;    //< defines the pointer type of the iterator.
  typedef std::size_t    size_type;        //< defines the integral type that can hold the size of a sequence.
  typedef std::ptrdiff_t difference_type;  //< defines the signed integral type that can hold the distance between two iterators.
  typedef std::forward_iterator_tag iterator_category;      //< defines the iterator category.
  // \}

private:
  SCI _cell_iterator;
  Simplex_3 _curr_simplex;

public:
  Triangulation_segment_simplex_iterator_3(const Tr* tr
    , Vertex_handle s, Vertex_handle t)
    : _cell_iterator(tr, s, t)
  { set_curr_simplex_to_entry(); }
  Triangulation_segment_simplex_iterator_3(const Tr* tr
    , Vertex_handle s, const Point& t)
    : _cell_iterator(tr, s, t)
  { set_curr_simplex_to_entry(); }
  Triangulation_segment_simplex_iterator_3(const Tr* tr
    , const Point& s, Vertex_handle t, Cell_handle hint = Cell_handle())
    : _cell_iterator(tr, s, t, hint)
  { set_curr_simplex_to_entry(); }
  Triangulation_segment_simplex_iterator_3(const Tr* tr
    , const Point& s, const Point& t, Cell_handle hint = Cell_handle())
    : _cell_iterator(tr, s, t, hint)
  { set_curr_simplex_to_entry(); }
  Triangulation_segment_simplex_iterator_3(const Tr* tr
    , const Segment& seg, Cell_handle hint = Cell_handle())
    : _cell_iterator(tr, seg, hint)
  { set_curr_simplex_to_entry(); }
  Triangulation_segment_simplex_iterator_3(const Tr* tr)
    : _cell_iterator(tr)
    , _curr_simplex()
  {}

  bool operator==(const Simplex_iterator& sit) const
  {
    return sit._cell_iterator == _cell_iterator
        && sit._curr_simplex == _curr_simplex;
  }
  bool operator!=(const Simplex_iterator& sit) const
  {
    return sit._cell_iterator != _cell_iterator
        || sit._curr_simplex != _curr_simplex;
  }

  const Point&    source() const      { return _cell_iterator.source(); }
  const Point&    target() const      { return _cell_iterator.target(); }

  const Tr& triangulation() const     { return *_cell_iterator.triangulation(); }

private:
  Triangulation_segment_simplex_iterator_3
    (const SCI& sci)
    : _cell_iterator(sci)
    , _curr_simplex()
  {}

private:
  void set_curr_simplex_to_entry()
  {
#if CGAL_DEBUG_TRIANGULATION_SEGMENT_TRAVERSER_3
    std::cerr << "cell iterator is:\n" << _cell_iterator.debug_iterator() << std::endl;
#endif // #if CGAL_DEBUG_TRIANGULATION_SEGMENT_TRAVERSER_3

    Locate_type lt;
    int li, lj;
    Cell_handle cell = Cell_handle(_cell_iterator);

    //check what is the entry type of _cell_iterator
    if (cell == Cell_handle())
    {
      //where did the segment get out from previous cell
      cell = _cell_iterator.previous();
      _cell_iterator.exit(lt, li, lj);
    }
    else
    {
      _cell_iterator.entry(lt, li, lj);
    }

    switch (lt)
    {
    case Locate_type::VERTEX:
      _curr_simplex = cell->vertex(li);
      break;
    case Locate_type::EDGE:
      _curr_simplex = Edge(cell, li, lj);
      break;
    case Locate_type::FACET:
      _curr_simplex = Facet(cell, li);
      break;
      //the 3 cases below correspond to the case when _cell_iterator
      //is in its initial position: _cur is locate(source)
    case Locate_type::CELL:
    case Locate_type::OUTSIDE_CONVEX_HULL:
    case Locate_type::OUTSIDE_AFFINE_HULL:
      if (Cell_handle(_cell_iterator) == Cell_handle())
        _curr_simplex = Simplex_3();
      else
        _curr_simplex = cell;
      break;
    default:
      CGAL_unreachable();
    };
  }

public:
  Simplex_iterator end() const
  {
    Simplex_iterator sit(_cell_iterator.end());
    return sit;
  }

  //  provides the increment postfix operator.
  Simplex_iterator& operator++()
  {
    auto increment_cell_iterator = [&]() {
      ++_cell_iterator;
#if CGAL_DEBUG_TRIANGULATION_SEGMENT_TRAVERSER_3
      std::cerr << "increment cell iterator to:\n" << _cell_iterator.debug_iterator() << '\n';
#endif
    };
    CGAL_assertion(_curr_simplex.incident_cell() != Cell_handle());

    if(!cell_iterator_is_ahead()) {
      increment_cell_iterator(); // cell_iterator needs to be ahead
    }

    Cell_handle ch_next = Cell_handle(_cell_iterator);
    Cell_handle ch_prev = _cell_iterator.previous();
    Locate_type lt_prev;
    int li_prev, lj_prev;
    _cell_iterator.exit(lt_prev, li_prev, lj_prev);

    if(_curr_simplex.dimension() == 3) {
      set_curr_simplex_to_entry();
      return *this;
    }
    if(lt_prev == Locate_type::CELL ||
       lt_prev == Locate_type::OUTSIDE_CONVEX_HULL ||
       lt_prev == Locate_type::OUTSIDE_AFFINE_HULL)
    {
      CGAL_assertion(ch_next == Cell_handle());
      _curr_simplex = ch_prev;
      return *this;
    }

    switch(_curr_simplex.dimension()) {
    case 2: { /*Facet*/
      CGAL_assertion((ch_next == Cell_handle()) == (_cell_iterator == _cell_iterator.end()));

      switch(lt_prev) {
      case Locate_type::VERTEX: { // facet-cell?-vertex-outside
        Vertex_handle v_prev{ch_prev->vertex(li_prev)};
        if(facet_has_vertex(get_facet(), v_prev))
          _curr_simplex = v_prev;
        else
          _curr_simplex = ch_prev;
      } break;
      case Locate_type::EDGE: { // facet-cell?-edge-outside
        Edge edge_prev{ch_prev, li_prev, lj_prev};
        if(facet_has_edge(get_facet(), edge_prev))
          _curr_simplex = edge_prev;
        else
          _curr_simplex = ch_prev;
      } break;
      case Locate_type::FACET: { // facet-cell-facet-outside
        Facet f_prev{ch_prev, li_prev};
        if(is_same_facet(f_prev, get_facet())) {
          if(ch_next == Cell_handle())
            _curr_simplex = Simplex_3();
          else
            _curr_simplex = ch_next;
        } else
          _curr_simplex = ch_prev;
      } break;
      default:
        CGAL_unreachable();
      }
    } break;
    case 1: {/*Edge*/
      switch(lt_prev) {
      case Locate_type::VERTEX: { //edge-vertex-outside
        Vertex_handle v_prev{ch_prev->vertex(li_prev)};
        if(edge_has_vertex(get_edge(), v_prev))
          _curr_simplex = v_prev;
        else
          _curr_simplex = shared_facet(get_edge(), v_prev);
      } break;
      case Locate_type::EDGE: { //edge-outside or edge-cell-edge-outside
        const Edge e_prev(ch_prev, li_prev, lj_prev);
        if(is_same_edge(get_edge(), e_prev)) {
          if(ch_next == Cell_handle()) {
            _curr_simplex = Simplex_3();
          } else {
            _curr_simplex = ch_next;
          }
        } else {
          auto facet_opt = shared_facet(get_edge(), e_prev);
          if(static_cast<bool>(facet_opt)) {
            _curr_simplex = *facet_opt;
          }
          else {
            _curr_simplex = shared_cell(get_edge(), e_prev);
          }
        }
      } break;
      case Locate_type::FACET: {
        Facet f_prev{ch_prev, li_prev};
        if(facet_has_edge(f_prev, get_edge()))
          _curr_simplex = f_prev; //edge-facet-outside
        else
          _curr_simplex = ch_prev; //query goes through the cell
      } break;
      default:
        CGAL_unreachable();
      }
    } break;
    case 0 :/*Vertex_handle*/
    {
      switch(lt_prev) {
      case Locate_type::VERTEX: {
        if(ch_prev->vertex(li_prev) != get_vertex()) // avoid infinite loop edge-vertex-same edge-...
          _curr_simplex = Edge(ch_prev, li_prev, ch_prev->index(get_vertex()));
        else {
          if(ch_next == Cell_handle()) {
            _curr_simplex = Simplex_3();
          } else {
            _curr_simplex = ch_next;
          }
        }
      } break;
      case Locate_type::EDGE: {
        const Edge e_prev(ch_prev, li_prev, lj_prev);
        if(edge_has_vertex(e_prev, get_vertex()))
          _curr_simplex = e_prev;
        else
          _curr_simplex = shared_facet(Edge(ch_prev, li_prev, lj_prev), get_vertex());
      } break;
      case Locate_type::FACET: {
        if(ch_prev->vertex(li_prev) != get_vertex()) // vertex-facet-outside
          _curr_simplex = Facet(ch_prev, li_prev);
        else // vertex-cell-facet-outside
          _curr_simplex = ch_prev;
      } break;
      default:
        CGAL_unreachable();
      }
    } break;
    default:
      CGAL_unreachable();
    };
    return *this;
  }

  //  provides the increment prefix operator.
  Simplex_iterator operator++(int)
  {
    Simplex_iterator tmp(*this);
    ++(*this);
    return tmp;
  }

  //  provides a dereference operator.
  /*         \return a pointer to the current cell.
  */
  const Simplex_3*   operator->()        { return &_curr_simplex; }

  //  provides an indirection operator.
  /*  \return the current cell.
  */
  const Simplex_3&   operator*()         { return _curr_simplex; }

  //  provides a conversion operator.
  /*  \return the current simplex
  */
  operator const Simplex_3() const { return _curr_simplex; }

  bool is_vertex() const { return _curr_simplex.dimension() == 0; }
  bool is_edge()   const { return _curr_simplex.dimension() == 1; }
  bool is_facet()  const { return _curr_simplex.dimension() == 2; }
  bool is_cell()   const { return _curr_simplex.dimension() == 3; }

  const Cell cell() const
  {
    return _cell_iterator.cell();
  }

  const Simplex_3& get_simplex() const { return _curr_simplex; }
  Vertex_handle get_vertex() const
  {
    CGAL_assertion(is_vertex());
    return Vertex_handle(_curr_simplex);
  }
  Edge get_edge() const
  {
    CGAL_assertion(is_edge());
    return Edge(_curr_simplex);
  }
  Facet get_facet() const
  {
    CGAL_assertion(is_facet());
    return Facet(_curr_simplex);
  }
  Cell_handle get_cell() const
  {
    CGAL_assertion(is_cell());
    return Cell_handle(_curr_simplex);
  }

public:
  //returns true in any of the degenerate cases,
  //i.e. when _curr_simplex has the following values successively
  // edge / facet / edge
  // edge / facet / vertex
  // vertex / facet / edge
  // vertex / edge / vertex
  // TODO : rename this function
  bool is_collinear() const
  {
    int curr_dim = _curr_simplex.dimension();
    //this concerns only edges and facets
    if (curr_dim == 1 || curr_dim == 2)
      return cell_iterator_is_ahead();
      //the degeneracy has been detected by moving cell_iterator forward
    else
      return false;
  }

  int simplex_dimension() const
  {
    return _curr_simplex.dimension();
  }

private:
  bool cell_iterator_is_ahead() const
  {
    Cell_handle ch = Cell_handle(_cell_iterator);
    if(ch == Cell_handle())
      return true;

    switch (_curr_simplex.dimension())
    {
    case 0 ://vertex
      return !ch->has_vertex(get_vertex());
    case 1 ://edge
      return !cell_has_edge(ch, get_edge());
    case 2 ://facet
      return !cell_has_facet(ch, get_facet());
    case 3 ://cell
      return ch != get_cell();
    default:
      CGAL_unreachable();
    }
    //should not be reached
    CGAL_unreachable();
    return false;
  }

  bool cell_has_edge(const Cell_handle ch, const Edge& e) const
  {
    Vertex_handle v1 = e.first->vertex(e.second);
    Vertex_handle v2 = e.first->vertex(e.third);
    return ch->has_vertex(v1) && ch->has_vertex(v2);
  }
  bool cell_has_facet(const Cell_handle c, const Facet& f) const
  {
    return f.first == c
        || f.first->neighbor(f.second) == c;
  }

  bool facet_has_edge(const Facet& f, const Edge& e) const
  {
    Vertex_handle v1 = e.first->vertex(e.second);
    Vertex_handle v2 = e.first->vertex(e.third);
    Cell_handle c = f.first;
    const int fi = f.second;

    unsigned int count = 0;
    for (int i = 1; i < 4; ++i)
    {
      Vertex_handle vi = c->vertex((fi + i) % 4);
      if (vi == v1 || vi == v2)
        ++count;
      if (count == 2)
        return true;
    }
    return false;
  }

  bool facet_has_vertex(const Facet& f, const Vertex_handle v) const
  {
    return triangulation().tds().has_vertex(f, v);
  }

  bool edge_has_vertex(const Edge& e, const Vertex_handle v) const
  {
    return e.first->vertex(e.second) == v
        || e.first->vertex(e.third) == v;
  }

  bool is_same_edge(const Edge& e1, const Edge& e2) const
  {
    return edge_has_vertex(e1, e2.first->vertex(e2.second))
        && edge_has_vertex(e1, e2.first->vertex(e2.third));
  }

  bool is_same_facet(const Facet& f1, const Facet& f2) const
  {
    return f1 == f2 || triangulation().mirror_facet(f1) == f2;
  }

  boost::optional<Vertex_handle> shared_vertex(const Edge& e1, const Edge& e2) const
  {
    Vertex_handle v1a = e1.first->vertex(e1.second);
    Vertex_handle v1b = e1.first->vertex(e1.third);
    Vertex_handle v2a = e2.first->vertex(e2.second);
    Vertex_handle v2b = e2.first->vertex(e2.third);

    if (v1a == v2a || v1a == v2b)
      return v1a;
    else if (v1b == v2a || v1b == v2b)
      return v1b;
    else
      return {};
  }

  boost::optional<Facet> shared_facet(const Edge& e1, const Edge& e2) const
  {
    Vertex_handle v2a = e2.first->vertex(e2.second);
    Vertex_handle v2b = e2.first->vertex(e2.third);

    auto sv_opt = shared_vertex(e1, e2);
    if(!sv_opt) return {};
    Vertex_handle sv = *sv_opt;
    Vertex_handle nsv2 = (sv == v2a) ? v2b : v2a;

    typename Tr::Facet_circulator circ
      = triangulation().incident_facets(e1);
    typename Tr::Facet_circulator end = circ;
    do
    {
      Facet f = *circ;
      for (int i = 1; i < 4; ++i)
      {
        if (nsv2 == f.first->vertex((f.second + i) % 4))
          return f;
      }
    } while (++circ != end);

    return {};
  }

  Facet shared_facet(const Edge& e, const Vertex_handle v) const
  {
    typename Tr::Facet_circulator circ
      = triangulation().incident_facets(e);
    typename Tr::Facet_circulator end = circ;
    do
    {
      Facet f = *circ;
      if (facet_has_vertex(f, v))
        return f;
    } while (++circ != end);

    std::cerr << "There is no facet shared by e and v" << std::endl;
    CGAL_unreachable();
    return Facet(Cell_handle(), 0);
  }

  Cell_handle shared_cell(const Edge& e, const Vertex_handle v) const
  {
    typename Tr::Cell_circulator circ
      = triangulation().incident_cells(e);
    typename Tr::Cell_circulator end = circ;
    do
    {
      Cell_handle c = circ;
      if (c->has_vertex(v))
        return c;
    } while (++circ != end);

    std::cerr << "There is no cell shared by e and v" << std::endl;
    CGAL_unreachable();
    return Cell_handle();
  }

  Cell_handle shared_cell(const Facet& f, const Vertex_handle v) const
  {
    Cell_handle c = f.first;
    if (c->has_vertex(v))
      return c;
    else
    {
      c = f.first->neighbor(f.second);
      CGAL_assertion(c->has_vertex(v));
      return c;
    }
  }

  Cell_handle shared_cell(const Edge e1, const Edge e2) const {
    auto facet = shared_facet(e1, e2.first->vertex(e2.second));
    return shared_cell(facet, e2.first->vertex(e2.third));
  }

};//class Triangulation_segment_simplex_iterator_3

} // namespace CGAL

#include <CGAL/Triangulation_3/internal/Triangulation_segment_traverser_3_impl.h>

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
