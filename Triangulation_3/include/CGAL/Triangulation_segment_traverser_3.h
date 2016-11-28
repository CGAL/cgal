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


#include <CGAL/basic.h>

#include <iostream>
#include <list>
#include <set>
#include <map>
#include <utility>
#include <stack>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

#include <CGAL/tuple.h>

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
 *  \cgalModels{ForwardIterator}
 *
 *	\sa `Triangulation_3`
 *	\sa `Forward_circulator_base`
 */
template < class Tr_, class Inc = internal::Incrementer<Tr_> >
class Triangulation_segment_cell_iterator_3 {
    typedef	Tr_                                         Tr;
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

    typedef typename Tr::Vertex_handle                  Vertex_handle;          //< defines the type of a handle for a vertex in the triangulation.
    typedef typename Tr::Cell_handle                    Cell_handle;            //< defines the type of a handle for a cell in the triangulation.

    typedef typename Tr::Locate_type                    Locate_type;            //< defines the simplex type returned from location.

    typedef CGAL::cpp11::tuple<Cell_handle,Locate_type,int,int>
                                                        Simplex;                //< defines the simplex type.

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
    const Tr&       _tr;

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
     *	\param s the source vertex. This vertex must be initialized and cannot be the infinite vertex.
     *	\param t the target vertex. This vertex must be initialized and cannot be the infinite vertex.
     *	It cannot equal `s`.
     */
	Triangulation_segment_cell_iterator_3( const Tr& tr, Vertex_handle s, Vertex_handle t );

    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *	\param s the source vertex. This vertex must be initialized and cannot be the infinite vertex.
     *	\param t the target point. This point must be initialized and it cannot be be at the same location as `s`.
     *	If `tr` has dimension < 3, `t` must lie inside the affine hull of `tr`.
     */
	Triangulation_segment_cell_iterator_3( const Tr& tr, Vertex_handle s, const Point& t );

    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *	\param s the source point. This point must be initialized and it cannot be be at the same location as `t`.
     *	\param t the target vertex. This vertex must be initialized and cannot be the infinite vertex.
     *	If `tr` has dimension < 3, `s` must lie inside the affine hull of `tr`.
     *	\param hint the starting point to search for `s`.
     */
	Triangulation_segment_cell_iterator_3( const Tr& tr, const Point& s, Vertex_handle t, Cell_handle hint = Cell_handle() );

    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *	\param s the source point. This point must be initialized. If `tr` has dimension < 3, `s` must lie inside
     *	the affine hull of `tr`.
     *	\param t the target point. This point must be initialized and it cannot be be at the same location as `s`.
     *	If `tr` has dimension < 3, `t` must lie inside the affine hull of `tr`.
     *	\param hint the starting point to search for `s`.
     */
    Triangulation_segment_cell_iterator_3( const Tr& tr, const Point& s, const Point& t, Cell_handle hint = Cell_handle() );
    
    //  constructs an iterator.
    /*  \param tr the triangulation to iterate though. This triangulation must have dimension > 0.
     *	\param S the segment to be traversed. If `tr` has dimension < 3, `S` must lie inside
     *	the affine hull of `tr`. `S` must not be degenerate, i.e. its source and target must not be equal.
     *	\param hint the starting point to search for `S`.
     */
    Triangulation_segment_cell_iterator_3( const Tr& tr, const Segment& S, Cell_handle hint = Cell_handle() );
// \}

    // The virtual destructor is mainly defined to indicate to the casting
    // operators that this is a dynamic type.
    virtual ~Triangulation_segment_cell_iterator_3() {}
    
private:
    //  private constructor that does not initialize the source and target.
    Triangulation_segment_cell_iterator_3( const Tr& tr );

public:
// \name Accessors
// \{
    //  gives the source point of the segment followed.
    /*  \return the source point.
     */
    const Point&    source() const      { return _source; }

    //  gives the target point of the segment follwoed.
    /*  \return the target point.
	 */
    const Point&    target() const      { return _target; }

    //  gives the current cell.
    /*  By invariance, this cell is intersected by the segment
	 *	between `source()` and `target()`.
	 *	\return the current cell.
	 *	\sa `handle()`.
	 */
    const Cell      cell() const        { return *get<0>(_cur); }

    //  gives a handle to the current cell.
    /*  By invariance, this cell is intersected by the segment
	 *	between `source()` and `target()`.
	 *	\return a handle to the current cell.
	 *	\sa `cell()`.
	 */
    Cell_handle     handle()            { return get<0>(_cur); }

    //  gives the previous cell.
    /*  This cell is uninitialized until the iterator leaves the initial
	 *	cell.
	 *	By invariance, once initialized, this cell must be intersected by the segment
	 *	between `source()` and `target()`.
	 *	\return the previous cell.
	 *	\sa `handle()`.
	 */
    Cell_handle     previous() const    { return get<0>(_prev); }

    //  provides a dereference operator.
    /* 	\return a pointer to the current cell.
	 */
    Cell*           operator->()        { return &*get<0>(_cur); }

    //  provides an indirection operator.
    /*  \return the current cell.
	 */
    Cell&           operator*()         { return *get<0>(_cur); }

    //  provides a conversion operator.
    /* 	\return a handle to the current cell.
	 */
    operator const  Cell_handle() const { return get<0>(_cur); }

    //  provides a conversion operator.
    /* 	\return the simplex through wich the current cell was entered.
	 */
    operator const  Simplex() const { return _cur; }

    //  checks whether the iterator has reached the final cell, which contains the `target()`.
    /*  If the `target()` lies on a facet, edge, or vertex, the final cell is the cell containing
	 *	the interior of the segment between `source()` and `target()`.
	 *	\return true iff the current cell contains the `target()`.
	 */
    bool            has_next() const    { return get<0>(_cur) != Cell_handle(); }

    //  gives the simplex through which the current cell was entered.
    /* 	For the first cell, containing the `source()` \f$ s \f$,
     *  this indicates the location of \f$ s \f$ in this cell.
	 */
    void            entry( Locate_type& lt, int& li, int& lj ) const    { lt = get<1>(_cur); li = get<2>(_cur); lj = get<3>(_cur); }

    //  gives the simplex through which the previous cell was exited.
    /* 	\pre the current cell is not the initial cell.
	 */
    void            exit( Locate_type& lt, int& li, int& lj ) const    { lt = get<1>(_prev); li = get<2>(_prev); lj = get<3>(_prev); }

    //  gives the past-the-end iterator associated with this iterator.
    SCI             end() const;
// \}

public:
// \name Mutators
// \{
    //  provides the increment postfix operator.
    /* 	After incrementing the iterator, the current cell intersects the segment
     *	between `source()` and `target()` closer to the `target()` than the previous cell.
     *	\sa `operator++(int)`.
     *  \pre The current cell does not contain the `target()`.
     */
    SCI&            operator++();
		
    //  provides the increment prefix operator.
    /* 	After incrementing the iterator, the current cell intersects the segment
     *	between `source()` and `target()` closer to the `target()` than the previous cell.
     *	than the previous cell.
     *	\sa `operator++()`.
     *  \pre The current cell does not contain the `target()`.
     */
    SCI             operator++( int );

    //  iterates to the final cell, which contains the `target()`.
    /* 	\return the final cell.
     */
    Cell_handle     complete();
// \}
	
public:
// \name Comparison
// \{
    //  compares this iterator with `sci`.
    /*  \param sci the other iterator.
     *	\return true iff the other iterator iterates the same triangulation along the same line segment
     *	and has the same current cell.
     *	\sa `operator!=( const SCI& t )`.
     */
    bool            operator==( const SCI& sci ) const;
		
    //  compares this iterator with `sci`.
    /*  \param sci the other iterator.
     *	\return `false` iff the other iterator iterates the same triangulation along the same line segment
     *	and has the same current cell.
     *	\sa `operator==( const SCI& t ) const`.
     */
    bool            operator!=( const SCI& sci ) const;

    //  compares the current cell with `ch`.
    /*  \param ch a handle to the other cell.
     *	\return true iff the current cell is the same as the one pointed to by `ch`.
     *	\sa `operator!=( const Cell_handle& ch ) const`.
     *	\sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
     */
    bool            operator==( const Cell_handle& ch ) const   { return ch == get<0>(_cur); }
		
    //  compares the current cell with `ch`.
    /*  \param ch a handle to the other cell.
     *	\return `false` iff the current cell is the same as the one pointed to by `ch`.
     *	\sa `operator==( const Cell_handle& ch )`.
     *	\sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
     */
    bool            operator!=( const Cell_handle& ch ) const   { return ch != get<0>(_cur); }
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
        typedef Incrementer::SCI    Expected;
#ifdef CGAL_TST_ASSUME_CORRECT_TYPES
        Expected& sci = static_cast<Expected&>( *this );
#else // CGAL_TST_ASSUME_CORRECT_TYPES
        Expected& sci = dynamic_cast<Expected&>( *this );
#endif // CGAL_TST_ASSUME_CORRECT_TYPES
        Incrementer().increment( sci );
    }
// \}

private:
    //  walk_to_next(), if the triangulation is 3D.
    void            walk_to_next_3();
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
}; // class Triangulation_segment_cell_iterator_3
	
//  compares a handle to a cell to a traverser.
/*  \param ch the handle to a cell.
 *	\param t the traverser.
 *	\return true iff the cell currently traversed by `t` is the same as the one pointed to by `ch`.
 *	\sa `operator!=( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
 *	\sa `Triangulation_segment_cell_iterator_3::operator==( const Cell_handle& ch )`.
 */
template < class Tr, class Inc >
inline bool operator==( typename Tr::Cell_handle ch, Triangulation_segment_cell_iterator_3<Tr,Inc> tci ) { return tci == ch; }

//  compares a handle to a cell to a traverser.
/*  \param ch the handle to a cell.
 *	\param t the traverser.
 *	\return `false` iff the cell currently traversed by `t` is the same as the one pointed to by `ch`.
 *	\sa `operator==( typename TriangulationTraits_3::Cell_handle ch, Triangulation_segment_cell_iterator_3<TriangulationTraits_3> t )`.
 *	\sa `Triangulation_segment_cell_iterator_3::operator!=( const Cell_handle& ch )`.
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
  typedef Triangulation_segment_cell_iterator_3<Tr_, Inc> SCI;

private:
  typedef typename SCI::Point    Point;
  typedef typename SCI::Segment  Segment;

public:
  // \{
  typedef typename SCI::Vertex_handle Vertex_handle;//< defines the type of a handle for a vertex in the triangulation
  typedef typename SCI::Cell_handle   Cell_handle;  //< defines the type of a handle for a cell in the triangulation.
  typedef typename SCI::Triangulation::Edge  Edge;  //< defines the type of an edge in the triangulation.
  typedef typename SCI::Triangulation::Facet Facet; //< defines the type of a facet in the triangulation.
  typedef typename SCI::Locate_type   Locate_type;  //< defines the simplex type returned from location.

  typedef boost::variant<Vertex_handle, Edge, Facet, Cell_handle> Simplex_type; //types sorted by dimension

  typedef Simplex_type   value_type;       //< defines the value type  the iterator refers to.
  typedef Simplex_type&  reference;        //< defines the reference type of the iterator.
  typedef Simplex_type*  pointer;          //< defines the pointer type of the iterator.
  typedef std::size_t    size_type;        //< defines the integral type that can hold the size of a sequence.
  typedef std::ptrdiff_t difference_type;  //< defines the signed integral type that can hold the distance between two iterators.
  typedef std::forward_iterator_tag iterator_category;      //< defines the iterator category.
  // \}

private:
  SCI _cell_iterator;
  Simplex_type _curr_simplex;

public:
  Triangulation_segment_simplex_iterator_3(const Tr& tr
    , Vertex_handle s, Vertex_handle t)
    : _cell_iterator(tr, s, t)
  { set_curr_simplex(); }
  Triangulation_segment_simplex_iterator_3(const Tr& tr
    , Vertex_handle s, const Point& t)
    : _cell_iterator(tr, s, t)
  { set_curr_simplex(); }
  Triangulation_segment_simplex_iterator_3(const Tr& tr
    , const Point& s, Vertex_handle t, Cell_handle hint = Cell_handle())
    : _cell_iterator(tr, s, t, hint)
  { set_curr_simplex(); }
  Triangulation_segment_simplex_iterator_3(const Tr& tr
    , const Point& s, const Point& t, Cell_handle hint = Cell_handle())
    : _cell_iterator(tr, s, t, hint)
  { set_curr_simplex(); }
  Triangulation_segment_simplex_iterator_3(const Tr& tr
    , const Segment& seg, Cell_handle hint = Cell_handle())
    : _cell_iterator(tr, seg, hint)
  { set_curr_simplex(); }

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

private:
  Triangulation_segment_simplex_iterator_3
    (const SCI& sci)
    : _cell_iterator(sci)
    , _curr_simplex(Cell_handle())
  {}

private:
  void set_curr_simplex()
  {
    //check what is the entry type of _cell_iterator
    if (Cell_handle(_cell_iterator) == Cell_handle())
    {
      _curr_simplex = Cell_handle();
      return;//end has been reached
    }

    Cell_handle cell = Cell_handle(_cell_iterator);
    Locate_type lt;
    int li, lj;
    _cell_iterator.entry(lt, li, lj);

    switch (lt)
    {
    case Locate_type::VERTEX:
      _curr_simplex = cell->vertex(li);
      break;

    case Locate_type::EDGE:
      _curr_simplex = Edge(cell, li, lj);
      break;
    case Locate_type::FACET: //basic case where segment enters a cell
                             //by crossing a facet
      //the 3 cases below correspond to the case when _cell_iterator
      //is in its initial position: _cur is locate(source)
    case Locate_type::CELL:
    case Locate_type::OUTSIDE_CONVEX_HULL:
    case Locate_type::OUTSIDE_AFFINE_HULL:
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
    CGAL_assertion(!_curr_simplex.empty());

    switch(_curr_simplex.which())
    {
    case 3 :/*Cell_handle*/
    {
      ++_cell_iterator;
      set_curr_simplex();
      break;
    }
    case 2 :/*Facet*/
    {
      if (!cell_iterator_is_ahead())
      {
        //cell_iterator is not ahead. get_facet() is part of cell_iterator
        //we cannot be in any of the degenerate cases, only detected by
        //taking cell_iterator one step forward
        _curr_simplex = Cell_handle(_cell_iterator);
      }
      else
      {
        Cell_handle chnext = Cell_handle(_cell_iterator);
        Locate_type ltnext;
        int linext, ljnext;
        _cell_iterator.entry(ltnext, linext, ljnext);
        switch (ltnext)//entry simplex in next cell
        {
        case Locate_type::VERTEX:
          //if the entry vertex is a vertex of current facet
          if (facet_has_vertex(get_facet(), chnext->vertex(linext)))
            set_curr_simplex();
          else
            _curr_simplex = chnext;
          break;

        case Locate_type::EDGE:
          if (facet_has_edge(get_facet(), Edge(chnext, linext, ljnext)))
            set_curr_simplex();
          else
            _curr_simplex = chnext;
          break;

        case Locate_type::FACET:
          _curr_simplex = chnext;
          break;

        default:
          CGAL_assertion(false);
        };
      }//end else
      break;
    }
    case 1:/*Edge*/
    {
      if (!cell_iterator_is_ahead())
      {
        //we cannot be in any of the degenerate cases, only detected by
        //taking cell_iterator one step forward
        _curr_simplex = Cell_handle(_cell_iterator);
      }
      else
      {
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
            _curr_simplex = get_edge();
          break;

        case Locate_type::EDGE:
        {
          //find the Facet between current and next edges
          Cell_handle ch = Cell_handle(_cell_iterator);
          Locate_type lt;
          int li, lj;
          _cell_iterator.entry(lt, li, lj);
          int index_v3 = -1;//todo : check this is correct
          int tmp = ch->index(chnext->vertex(li));
          if (tmp != li && tmp != lj)
            index_v3 = tmp;
          else
          {
            tmp = ch->index(chnext->vertex(lj));
            index_v3 = tmp;
          }
          CGAL_assertion(index_v3 != -1);
          _curr_simplex = Facet(ch, (6 - li - lj - index_v3));
          break;
        }
        case Locate_type::FACET:
          _curr_simplex = Cell_handle(_cell_iterator);//query goes through the cell
          break;
        };
      }//end else
      break;
    }
    case 0 :/*Vertex_handle*/
    {
      if (degeneracy())
      {
        set_curr_simplex();
        break;
      }
      Cell_handle ch = Cell_handle(_cell_iterator);
      if (cell_iterator_is_ahead())
      {
        _curr_simplex = ch;
        break;
      }

      ++_cell_iterator;
      Cell_handle chnext = Cell_handle(_cell_iterator);
      //_cell_iterator is one step forward _curr_simplex

      Locate_type lt;
      int li, lj;
      _cell_iterator.entry(lt, li, lj);

      int index_v = ch->index(get_vertex());
      int index_vnext = ch->index(chnext->vertex(li));

      if (lt == Locate_type::VERTEX)
      {
        _curr_simplex = Edge(ch, index_v, index_vnext);
      }
      else if (lt == Locate_type::EDGE)
      {
        int index_f = 6 - (index_v + index_vnext + ch->index(chnext->vertex(lj)));
        _curr_simplex = Facet(ch, index_f);
      }
      else
        _curr_simplex = Cell_handle(_cell_iterator);
      break;
    }
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
  /* 	\return a pointer to the current cell.
  */
  Simplex_type*   operator->()        { return &*_curr_simplex; }

  //  provides an indirection operator.
  /*  \return the current cell.
  */
  Simplex_type&   operator*()         { return _curr_simplex; }

  //  provides a conversion operator.
  /* 	\return a handle to the current cell.
  */
  operator const  Cell_handle() const { return Cell_handle(_cell_iterator); }

  bool is_vertex() const { return _curr_simplex.which() == 0; }
  bool is_edge()   const { return _curr_simplex.which() == 1; }
  bool is_facet()  const { return _curr_simplex.which() == 2; }
  bool is_cell()   const { return _curr_simplex.which() == 3; }

  const Simplex_type& get_simplex() const { return _curr_simplex; }
  Vertex_handle get_vertex() const
  {
    CGAL_assertion(is_vertex());
    return boost::get<Vertex_handle>(_curr_simplex);
  }
  Edge get_edge() const
  {
    CGAL_assertion(is_edge());
    return boost::get<Edge>(_curr_simplex);
  }
  Facet get_facet() const
  {
    CGAL_assertion(is_facet());
    return boost::get<Facet>(_curr_simplex);
  }
  Cell_handle get_cell() const
  {
    CGAL_assertion(is_cell());
    return boost::get<Cell_handle>(_curr_simplex);
  }

private:
  /** returns true when _current_simplex does not have the same
  * type as the "entry" of _cell_iterator
  * it is for example the case after vertex - edge.
  * the following value should be vertex, and the corresponding vertex
  * is the entry of _cell_iterator
  */
  bool degeneracy() const
  {
    Locate_type lt;
    int li, lj;
    _cell_iterator.entry(lt, li, lj);
    return (is_vertex() && lt != Locate_type::VERTEX)
        || (is_edge() && lt != Locate_type::EDGE)
        || (is_facet() && lt != Locate_type::FACET);
  }

  bool cell_iterator_is_ahead() const
  {
    Cell_handle ch = Cell_handle(_cell_iterator);
    switch (_curr_simplex.which())
    {
    case 0 ://vertex
      return !ch->has_vertex(get_vertex());
    case 1 ://edge
      return !cell_has_edge(ch, get_edge());
    case 2 ://facet
      return !cell_has_facet(ch, get_facet());
    case 3 ://cell
      return ch != get_cell();
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

  bool facet_has_vertex(const Facet& f, const Vertex_handle v) const
  {
    Cell_handle c = f.first;
    const int fi = f.second;
    for (int i = 1; i < 4; ++i)
    {
      if (c->vertex((fi + i) % 4) == v)
        return true;
    }
    return false;
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
      if (vi == v1) ++count;
      else if (vi == v2) ++count;
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

};//class Triangulation_segment_simplex_iterator_3

} // namespace CGAL

#include <CGAL/Triangulation_segment_traverser_3_impl.h>

#endif // CGAL_TRIANGULATION_SEGMENT_TRAVERSER_3_H
