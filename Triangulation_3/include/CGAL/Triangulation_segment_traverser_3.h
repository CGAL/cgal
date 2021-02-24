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
#include <tuple>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_simplex_3.h>


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
 *  \cgalModels{ForwardIterator}
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

    typedef std::tuple<Cell_handle,Locate_type,int,int> Simplex;                //< defines the simplex type.

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

    //  gives the target point of the segment follwoed.
    /*  \return the target point.
         */
    const Point&    target() const      { return _target; }

    //  gives a handle to the current cell.
    /*  By invariance, this cell is intersected by the segment
         *        between `source()` and `target()`.
         *        \return a handle to the current cell.
         *        \sa `cell()`.
         */
    Cell_handle     handle()
    {
      return std::get<0>(_cur);
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
      return &*std::get<0>(_cur);
    }

    //  provides an indirection operator.
    /*  \return the current cell.
         */
    Cell&           operator*()
    {
      return *std::get<0>(_cur);
    }

    //  provides a conversion operator.
    /*         \return a handle to the current cell.
         */
    operator Cell_handle() const
    {
      return std::get<0>(_cur);
    }

    //  provides a conversion operator.
    /*         \return the simplex through wich the current cell was entered.
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
    //  gives the simplex through which the previous cell was exited.
    /*         \pre the current cell is not the initial cell.
         */
    void            exit( Locate_type& lt, int& li, int& lj ) const
    {
      lt = prev_lt(); li = prev_li(); lj = prev_lj();
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
      return ch == std::get<0>(_cur);
    }

    //  compares the current cell with `ch`.
    /*  \param ch a handle to the other cell.
     *        \return `false` iff the current cell is the same as the one pointed to by `ch`.
     *        \sa `operator==( const Cell_handle& ch )`.
     *        \sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
     */
    bool            operator!=( const Cell_handle& ch ) const
    {
      return ch != std::get<0>(_cur);
    }
// \}

        bool            operator==( Nullptr_t CGAL_triangulation_assertion_code(n) ) const;
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
        CGAL_triangulation_precondition( i>=0 && i<=3 );
        CGAL_triangulation_precondition( j>=0 && j<=3 );
        CGAL_triangulation_precondition( i != j );
        return ( i==0 || j==0 ) ? i+j-1 : i+j;
    }

    bool have_same_entry(const Simplex& s1, const Simplex& s2) const;

    // Compute the orientation of a point compared to the oriented plane supporting a half-facet.
    CGAL::Orientation orientation(const Facet& f, const Point& p) const;

    bool coplanar(const Facet &f, const Point &p) const;

    // Gives the edge incident to the same cell that is not incident to any of the input vertices.
    Edge opposite_edge(Cell_handle c, int li, int lj) const;
    Edge opposite_edge(const Edge& e) const;

    // ref-accessors to the simplex, for use in internal code
    // access _cur
    Cell_handle& cell()             { return std::get<0>(_cur); }
    Cell_handle const& cell() const { return std::get<0>(_cur); }

    Locate_type& lt()             { return std::get<1>(_cur); }
    Locate_type const& lt() const { return std::get<1>(_cur); }

    int& li()             { return std::get<2>(_cur); }
    int const& li() const { return std::get<2>(_cur); }

    int& lj()             { return std::get<3>(_cur); }
    int const& lj() const { return std::get<3>(_cur); }

    // access _prev
    Cell_handle& prev_cell()             { return std::get<0>(_prev); }
    Cell_handle const& prev_cell() const { return std::get<0>(_prev); }

    Locate_type& prev_lt()             { return std::get<1>(_prev); }
    Locate_type const& prev_lt() const { return std::get<1>(_prev); }

    int& prev_li()             { return std::get<2>(_prev); }
    int const& prev_li() const { return std::get<2>(_prev); }

    int& prev_lj()             { return std::get<3>(_prev); }
    int const& prev_lj() const { return std::get<3>(_prev); }

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

  const Tr* triangulation() const     { return _cell_iterator.triangulation(); }

private:
  Triangulation_segment_simplex_iterator_3
    (const SCI& sci)
    : _cell_iterator(sci)
    , _curr_simplex()
  {}

private:
  void set_curr_simplex_to_entry()
  {
    Locate_type lt;
    int li, lj;
    Cell_handle cell;

    //check what is the entry type of _cell_iterator
    if (Cell_handle(_cell_iterator) == Cell_handle())
    {
      //where did the segment std::get out from previous cell
      cell = _cell_iterator.previous();
      _cell_iterator.exit(lt, li, lj);
    }
    else
    {
      cell = Cell_handle(_cell_iterator);
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
      CGAL_assertion(false);
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
    CGAL_assertion(_curr_simplex.incident_cell() != Cell_handle());

    switch(_curr_simplex.dimension())
    {
    case 3 :/*Cell_handle*/
    {
      Cell_handle ch = Cell_handle(_cell_iterator);
      if (ch == Cell_handle())
      {
        if(!triangulation()->is_infinite(Cell_handle(_curr_simplex)))
          set_curr_simplex_to_entry();
        else
          _curr_simplex = Simplex_3();
        break;
      }
      else
      {
        if (!cell_iterator_is_ahead())
          ++_cell_iterator;
        set_curr_simplex_to_entry();
      }
      break;
    }
    case 2 :/*Facet*/
    {
      Cell_handle ch = Cell_handle(_cell_iterator);
      if (!cell_iterator_is_ahead())
      {
        //cell_iterator is not ahead. get_facet() is part of cell_iterator
        //we cannot be in any of the degenerate cases, only detected by
        //taking cell_iterator one step forward
        CGAL_assertion(cell_has_facet(Cell_handle(_cell_iterator), get_facet()));
        ++_cell_iterator;
        if (Cell_handle(_cell_iterator) == Cell_handle())
        {
          _curr_simplex = _cell_iterator.previous();
          break;
        }
      }
      else
        ch = _cell_iterator.previous();

      Cell_handle chnext = Cell_handle(_cell_iterator);
      Locate_type ltnext;
      int linext, ljnext;
      _cell_iterator.entry(ltnext, linext, ljnext);
      switch (ltnext)//entry simplex in next cell
      {
      case Locate_type::VERTEX:
      {
        if (_cell_iterator == _cell_iterator.end())
        {
          _curr_simplex = Simplex_3();
          break;
        }
        //if the entry vertex is a vertex of current facet
        int i;
        if (triangulation()->has_vertex(get_facet(), chnext->vertex(linext), i))
          set_curr_simplex_to_entry();
        else
          _curr_simplex = ch;
        break;
      }

      case Locate_type::EDGE:
        if (facet_has_edge(get_facet(), Edge(chnext, linext, ljnext)))
          set_curr_simplex_to_entry();
        else
          _curr_simplex = ch;
        break;

      case Locate_type::FACET:
        _curr_simplex = ch;
        break;

      case Locate_type::OUTSIDE_AFFINE_HULL:
        _curr_simplex = Simplex_3();
        break;

      default:
        CGAL_assertion(false);
      };
      break;
    }
    case 1:/*Edge*/
    {
      Cell_handle ch = Cell_handle(_cell_iterator);
      if (ch == _cell_iterator.previous())
      {
        _curr_simplex = Simplex_3();
        break;
      }
      Locate_type lt;
      int li, lj;
      _cell_iterator.entry(lt, li, lj);

      if (!cell_iterator_is_ahead())
      {
        ++_cell_iterator;//cell_iterator needs to be ahead to detect degeneracies
        if (Cell_handle(_cell_iterator) == Cell_handle())
        {
          _curr_simplex = _cell_iterator.previous();
          break;
        }
      }

      Cell_handle chnext = Cell_handle(_cell_iterator);
      Locate_type ltnext;
      int linext, ljnext;
      _cell_iterator.entry(ltnext, linext, ljnext);
      switch (ltnext)//entry simplex in next cell
      {
      case Locate_type::VERTEX:
        if (edge_has_vertex(get_edge(), chnext->vertex(linext)))
          _curr_simplex = chnext->vertex(linext);
        else
          _curr_simplex = shared_facet(get_edge(), chnext->vertex(linext));
        break;

      case Locate_type::EDGE:
      {
        CGAL_assertion(_cell_iterator == _cell_iterator.end()
          || triangulation()->is_infinite(chnext)
          || _curr_simplex != Simplex_3(Edge(chnext, linext, ljnext)));

        if (_cell_iterator == _cell_iterator.end())
          _curr_simplex = Simplex_3();
        else if (triangulation()->is_infinite(chnext)
          && _curr_simplex == Simplex_3(Edge(chnext, linext, ljnext)))
          _curr_simplex = chnext;
        else
          _curr_simplex = shared_facet(get_edge(), Edge(chnext, linext, ljnext));
        break;
      }

      case Locate_type::FACET:
        _curr_simplex = Cell_handle(_cell_iterator);//query goes through the cell
        break;

      case Locate_type::OUTSIDE_AFFINE_HULL:
      {
        Cell_handle chprev = _cell_iterator.previous();
        Locate_type ltprev;
        int liprev, ljprev;
        _cell_iterator.exit(ltprev, liprev, ljprev);

        if (ltprev == Locate_type::VERTEX) //edge-vertex-outside
          _curr_simplex = chprev->vertex(liprev);
        else
          _curr_simplex = Simplex_3(); //edge-outside
        break;
      }
      default:
        CGAL_assertion(false);//should not happen
      };
      break;
    }
    case 0 :/*Vertex_handle*/
    {
      Cell_handle ch = Cell_handle(_cell_iterator);
      if (ch == _cell_iterator.previous())
      {
        _curr_simplex = Simplex_3();
        break;
      }
      if (!cell_iterator_is_ahead()) //_curr_simplex does contain v
      {
        ++_cell_iterator;//cell_iterator needs to be ahead to detect degeneracies
      }
      else
        ch = _cell_iterator.previous();

      Cell_handle chnext = Cell_handle(_cell_iterator);
      //_cell_iterator is one step forward _curr_simplex
      CGAL_assertion(ch != chnext);

      Locate_type ltnext;
      int linext, ljnext;
      _cell_iterator.entry(ltnext, linext, ljnext);

      Cell_handle prev;
      Locate_type ltprev;
      int liprev, ljprev;
      prev = _cell_iterator.previous();
      _cell_iterator.exit(ltprev, liprev, ljprev);

      switch (ltnext)
      {
      case Locate_type::VERTEX:
      {
        CGAL_assertion(_cell_iterator == _cell_iterator.end()
                     || get_vertex() != chnext->vertex(linext)
                     || triangulation()->is_infinite(chnext));
        if (_cell_iterator == _cell_iterator.end())
        {
          if (prev == ch && ltprev == Locate_type::VERTEX)
          {
            CGAL_assertion(prev->vertex(liprev) == get_vertex());
            _curr_simplex = ch;
          }
          else
          {
            if(ltprev == Locate_type::FACET)
              _curr_simplex = Facet(prev, liprev);
            else if(ltprev == Locate_type::EDGE)
              _curr_simplex = Edge(prev, liprev, ljprev);
            else
              CGAL_assertion(false);
          }
        }
        else
        {
          if (triangulation()->is_infinite(chnext) && get_vertex() == chnext->vertex(linext))
            _curr_simplex = chnext;
          else
          {
            Cell_handle ec;
            int ei = -1, ej = -1;
            if (!triangulation()->is_edge(get_vertex(), chnext->vertex(linext), ec, ei, ej))
              CGAL_assertion(false);
            _curr_simplex = Edge(ec, ei, ej);
          }
        }
        break;
      }

      case Locate_type::EDGE:
      {
        //facet shared by get_vertex() and the edge
        //none of ch and chnext is certainly shared by both endpoints
        _curr_simplex = shared_facet(Edge(chnext, linext, ljnext), get_vertex());
        break;
      }

      case Locate_type::OUTSIDE_AFFINE_HULL:
      {
        CGAL_assertion(_cell_iterator == _cell_iterator.end());
        if (ltprev == Locate_type::VERTEX) //vertex-edge-vertex-outside
        {
          if(prev->vertex(liprev) != get_vertex())//avoid infinite loop edge-vertex-same edge-...
            _curr_simplex = Edge(prev, liprev, prev->index(get_vertex()));
          else
            _curr_simplex = Simplex_3();
        }
        else if (ltprev == Locate_type::EDGE)//vertex-facet-edge-outside
          _curr_simplex = Facet(prev, prev->index(get_vertex()));
        else if (ltprev == Locate_type::FACET) //vertex-facet-outside
        {
          if(prev->vertex(liprev) != get_vertex()) //vertex-facet-outside
            _curr_simplex = Facet(prev, liprev);
          else //vertex-cell-facet-outside
            _curr_simplex = prev;
        }
        else
        {
          CGAL_assertion(ltprev == Locate_type::CELL);//vertex-cell-outside
          _curr_simplex = prev;
        }
        break;
      }

      default://FACET
        if (chnext == Cell_handle())
          _curr_simplex = Simplex_3();
        else
          _curr_simplex = shared_cell(Facet(chnext, linext), get_vertex());
        break;
      };
    }
    break;

    default:
      CGAL_assertion(false);
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
      CGAL_assertion(false);
    }
    //should not be reached
    CGAL_assertion(false);
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

  bool edge_has_vertex(const Edge& e, const Vertex_handle v) const
  {
    return e.first->vertex(e.second) == v
        || e.first->vertex(e.third) == v;
  }

  Vertex_handle shared_vertex(const Edge& e1, const Edge& e2) const
  {
    Vertex_handle v1a = e1.first->vertex(e1.second);
    Vertex_handle v1b = e1.first->vertex(e1.third);
    Vertex_handle v2a = e2.first->vertex(e2.second);
    Vertex_handle v2b = e2.first->vertex(e2.third);

    if (v1a == v2a || v1a == v2b)
      return v1a;
    else if (v1b == v2a || v1b == v2b)
      return v1b;

    std::cerr << "There is no vertex shared by e1 and e2" << std::endl;
    CGAL_assertion(false);
    return Vertex_handle();
  }

  Facet shared_facet(const Edge& e1, const Edge& e2) const
  {
    Vertex_handle v2a = e2.first->vertex(e2.second);
    Vertex_handle v2b = e2.first->vertex(e2.third);

    Vertex_handle sv = shared_vertex(e1, e2);
    Vertex_handle nsv2 = (sv == v2a) ? v2b : v2a;

    typename Tr::Facet_circulator circ
      = triangulation()->incident_facets(e1);
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

    std::cerr << "There is no facet shared by e1 and e2" << std::endl;
    CGAL_assertion(false);
    return Facet(Cell_handle(), 0);
  }

  Facet shared_facet(const Edge& e, const Vertex_handle v) const
  {
    typename Tr::Facet_circulator circ
      = triangulation()->incident_facets(e);
    typename Tr::Facet_circulator end = circ;
    do
    {
      Facet f = *circ;
      int i;
      if (triangulation()->has_vertex(f, v, i))
        return f;
    } while (++circ != end);

    std::cerr << "There is no facet shared by e and v" << std::endl;
    CGAL_assertion(false);
    return Facet(Cell_handle(), 0);
  }

  Cell_handle shared_cell(const Edge& e, const Vertex_handle v) const
  {
    typename Tr::Cell_circulator circ
      = triangulation()->incident_cells(e);
    typename Tr::Cell_circulator end = circ;
    do
    {
      Cell_handle c = circ;
      if (c->has_vertex(v))
        return c;
    } while (++circ != end);

    std::cerr << "There is no cell shared by e and v" << std::endl;
    CGAL_assertion(false);
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

};//class Triangulation_segment_simplex_iterator_3

} // namespace CGAL

#include <CGAL/internal/Triangulation_segment_traverser_3_impl.h>

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
