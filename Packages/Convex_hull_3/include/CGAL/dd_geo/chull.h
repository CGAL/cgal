// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 2000, August 21
// 
// source        : chull.fw
// file          : include/CGAL/dd_geo/chull.h
// package       : Convex_hull_3 (2.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.6
// revision_date : 21 Aug 2000 
// author(s)     : Kurt Mehlhorn
//                 Michael Seel
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 


#ifndef CGAL_DD_GEO_CHULL_H
#define CGAL_DD_GEO_CHULL_H
#define LEP_DDGEO_CHULL_H

#if !defined(LEP_DDGEO_INCL_ID)
#define LEDA_ROOT_INCL_ID NOT_ANY_KERNEL_NUMBER
#define LEP_DDGEO_INCL_ID 21010
#include <LEDA/REDEFINE_NAMES.h>
#endif


/*{\Msubst
EGCSNTA#
}*/

// remark: we templatize interface types directly, all other types
//         and operations are hidden within the traits class

/*{\Manpage {chull} {CHTRAITS} {Convex Hulls} {C}}*/

/*{\Mdefinition An instance |\Mvar| of type |\Mname| is the convex
hull of a multi-set |S| of points in $d$-dimensional space. We call
|S| the underlying point set and $d$ or |dim| the dimension of the
underlying space. We use |dcur| or |ducrrent| to denote the affine
dimension of |S|.  The data type supports incremental construction of
hulls.

The closure of the hull is maintained as a simplicial complex, i.e.,
as a collection of simplices the intersection of any two is a face of
both\footnote{The empty set if a facet of every simplex.}. In the
sequel we reserve the word simplex for the simplices of dimension
|dcur|. For each simplex there is an item of type |ch_simplex| and for
each vertex there is an item of type |ch_vertex|.  Each simplex has
|1+ dcur| vertices indexed from $0$ to |dcur|; for a simplex $s$ and
an index $i$, |C.vertex(s,i)| returns the $i$-th vertex of $s$. For
any simplex $s$ and any index $i$ of $s$ there may or may not be a
simplex $t$ opposite to the $i$-th vertex of $s$.  The function
|C.opposite_simplex(s,i)| returns $t$ if it exists and returns |nil|
otherwise. If $t$ exists then $s$ and $t$ share |dcur| vertices,
namely all but the vertex with index $i$ of $s$ and the vertex with
index |C.index_of_vertex_in_opposite_simplex(s,i)| of $t$. Assume that
$t$ exists and let |j =C.index_of_vertex_in_opposite_simplex(s,i)|.
Then $s =$ |C.opposite_simplex(t,j)| and $i =$
|C.index_of_vertex_in_opposite_simplex(t,j)|.

The boundary of the hull is also a simplicial complex. All simplices
in this complex have dimension |dcur - 1|.  For each boundary simplex
there is an item of type |ch_facet|.  Each facet has |dcur| vertices
indexed from $0$ to |dcur -1|. If |dcur > 1| then for each facet $f$
and each index $i$, $0 \le i < |dcur|$, there is a facet $g$ opposite
to the $i$-th vertex of $f$.  The function |C.opposite_facet(f,i)|
returns $g$.  Two neighboring facets $f$ and $g$ share |dcur - 1|
vertices, namely all but the vertex with index $i$ of $f$ and the
vertex with index |C.index_of_vertex_in_opposite_facet(f,i)| of $g$.
Let |j = C.index_of_vertex_in_opposite_facet(f,i)|. Then
|f = C.opposite_facet(g,j)| and
|i = C.index_of_vertex_in_opposite_facet(g,j)|.}*/

#include <LEDA/array.h>
#include <LEDA/list.h>
#include <LEDA/stack.h>
#include <LEDA/stream.h>
#include <LEDA/window.h>
#include <LEDA/graph.h>
// #include <LEDA/integer.h>
#include <LEDA/map.h>
#include <CGAL/dd_geo/config.h>
#include <CGAL/dd_geo/regl_complex.h>

//#define SELFCHECK
//#define DDGEO_STL_ITERATORS
//#define _DEBUG
#include <CGAL/dd_geo/debug.h>

#ifndef DDGEO_TEMPLATE_DEFAULTS
template <class CHTRAITS,class POINT,class PLANE> class chull;
#else
template <class CHTRAITS,
          class POINT = CGAL_TYPENAME_MSVC_NULL CHTRAITS::POINT,
          class PLANE = CGAL_TYPENAME_MSVC_NULL CHTRAITS::PLANE>
class chull;
#endif


#ifndef DDGEO_TEMPLATE_DEFAULTS
template <class CHTRAITS,class POINT,class PLANE>
#else
template <class CHTRAITS,
          class POINT = CGAL_TYPENAME_MSVC_NULL CHTRAITS::POINT,
          class PLANE = CGAL_TYPENAME_MSVC_NULL CHTRAITS::PLANE>
#endif
class ch_Simplex : public rc_Simplex<CHTRAITS,POINT>
{


  typedef ch_Simplex<CHTRAITS,POINT,PLANE>* ch_simplex;
  typedef ch_Simplex<CHTRAITS,POINT,PLANE>* ch_facet;
  typedef rc_Vertex<CHTRAITS,POINT>*        ch_vertex;
  friend class chull<CHTRAITS,POINT,PLANE>;
  // chull has unrestricted access

  PLANE base_facet;  // hyperplane supporting base facet
  bool  visited;     // visited mark when traversing the triangulation

  ch_Simplex(int dmax)  : rc_Simplex<CHTRAITS,POINT>(dmax)
  { visited = false; }

  ~ch_Simplex() {}  // cleanup in ~rc_Simplex suffices

  LEDA_MEMORY(ch_Simplex);
};


template <class CHTRAITS,class POINT,class PLANE>
class chull : public regl_complex<CHTRAITS,POINT>
{
public:
/*{\Mtypes 4}*/
typedef ch_Simplex<CHTRAITS,POINT,PLANE>* ch_simplex;
/*{\Mtypemember the item type for simplices of the complex.}*/

typedef ch_Simplex<CHTRAITS,POINT,PLANE>* ch_facet;
/*{\Mtypemember the item type for facets of the complex.}*/

typedef rc_Vertex<CHTRAITS,POINT>*        ch_vertex;
/*{\Mtypemember the item type for vertices of the complex.}*/

  typedef typename CHTRAITS::VECTOR       VECTOR;
  typedef typename CHTRAITS::IPOINT       IPOINT;
  typedef typename CHTRAITS::IRAY         IRAY;
  typedef typename CHTRAITS::RT           RT;
  // make traits types locally available

  #ifdef DDGEO_STL_ITERATORS
  /*{\Mtext \textbf{Iterator} types and operations are only
  defined if compiled including [[-DDDGEO_STL_ITERATORS]].}*/

  class ch_point_iterator;
  class ch_simplex_iterator;
  class ch_facet_iterator;
  friend class ch_point_iterator;
  friend class ch_simplex_iterator;

  class ch_point_iterator
  /*{\Mtypemember the const iterator for points.}*/
  {
    friend class chull<CHTRAITS,POINT,PLANE>;
    const chull<CHTRAITS,POINT,PLANE>* CH;
    list_item it;
    ch_point_iterator(const chull<CHTRAITS,POINT,PLANE>* x,
    list_item y) : CH(x), it(y)  {}
  public:
    ch_point_iterator() : CH(0),it(0) {}
    ch_point_iterator& operator++()
    { it = CH->all_pnts.stl_next_item(it); return *this; }
    ch_point_iterator& operator--()
    { it = CH->all_pnts.stl_pred_item(it); return *this; }
    ch_point_iterator  operator++(int)
    { ch_point_iterator tmp = *this; ++(*this); return tmp; }
    ch_point_iterator  operator--(int)
    { ch_point_iterator tmp = *this; --(*this); return tmp; }
    const POINT& operator*() const { return CH->all_pnts[it]; }
    list_item* operator->() { return &it; }
    friend bool operator==(const ch_point_iterator& x,
                           const ch_point_iterator& y)
    { return x.it == y.it; }

    typedef POINT  value_type;
    typedef POINT* pointer;
    typedef POINT& reference;
    typedef bidirectional_iterator_tag iterator_category;
    typedef ptrdiff_t difference_type;
  };


  typedef rc_vertex_iterator ch_vertex_iterator;
  /*{\Mtypemember the const iterator for vertices.}*/

  class ch_simplex_iterator
  /*{\Mtypemember the const iterator for simplices.}*/
  {
    friend class chull<CHTRAITS,POINT,PLANE>;
    const chull<CHTRAITS,POINT,PLANE>* CH;
    list_item it;
    ch_simplex_iterator(const chull<CHTRAITS,POINT,PLANE>* x,
    list_item y) : CH(x), it(y)  {}
  public:
    ch_simplex_iterator() : CH(0),it(0) {}
    ch_simplex_iterator& operator++()
    { it = CH->all_simps.stl_next_item(it); return *this; }
    ch_simplex_iterator& operator--()
    { it = CH->all_simps.stl_pred_item(it); return *this; }
    ch_simplex_iterator  operator++(int)
    { ch_simplex_iterator tmp = *this; ++(*this); return tmp; }
    ch_simplex_iterator  operator--(int)
    { ch_simplex_iterator tmp = *this; --(*this); return tmp; }
    const ch_simplex& operator*() const
    { return (const ch_simplex&)CH->all_simps[it]; }
    list_item* operator->() { return &it; }
    friend bool operator==(const ch_simplex_iterator& x,
                           const ch_simplex_iterator& y)
    { return x.it == y.it; }

    typedef ch_simplex  value_type;
    typedef ch_simplex* pointer;
    typedef ch_simplex& reference;
    typedef bidirectional_iterator_tag iterator_category;
    typedef ptrdiff_t   difference_type;
  };


  class ch_facet_iterator {
  /*{\Mtypemember the const iterator for facets.}*/
    friend class chull<CHTRAITS,POINT,PLANE>;
    const chull<CHTRAITS,POINT,PLANE>* CH;
    ch_facet current;
    stack<ch_facet> candidates;
    map<ch_facet,bool> visited;
    // invars: visited[current] = true
    //         visited[candidates] = true;
    //         all facet neighbors of current are
    //         in candidates or have left it already
    ch_facet_iterator(const chull<CHTRAITS,POINT,PLANE>* x,
    ch_facet y) : CH(x), current(y), visited(false)
    /* if the facet is not nil we set the current marker to
       the facet and insert all it's neighbors into the
       candidates stack */
    { if (current != nil) {
        visited[current] = true;
        ASSERT(!CH->vertex_of_simplex(current,0), facet_iterator);
        if ( CH->dcurrent() > 1 )
          for(int i = 1; i<= CH->dcurrent(); ++i) {
            ch_facet f = CH->opposite_simplex(current,i);
            if (!visited[f]) {
              candidates.push(f);
              visited[f]=true;
            }
           }
      }
    }
  public:
    ch_facet_iterator() : CH(0),current(0) {}
    ch_facet_iterator& operator++()
    /* here we get a new candidate from the stack
       and insert all its not-visited neighbors */
    { if (!candidates.empty()) {
        current = candidates.pop();
        ASSERT(!CH->vertex_of_simplex(current,0), facet_iterator);
        for(int i = 1; i<= CH->dcurrent(); ++i) {
          ch_facet f = CH->opposite_simplex(current,i);
          if (!visited[f]) {
            candidates.push(f);
            visited[f]=true;
          }
        }
      }
      else
        current = nil;
      return *this;
    }
    ch_facet_iterator  operator++(int)
    { ch_facet_iterator tmp = *this; ++(*this); return tmp; }
    const ch_facet& operator*() const
    { return (const ch_facet&) current; }
    list_item* operator->() { return &current; }
    friend bool operator==(const ch_facet_iterator& x,
                           const ch_facet_iterator& y)
    { return x.current == y.current; }

    typedef ch_facet  value_type;
    typedef ch_facet* pointer;
    typedef ch_facet& reference;
    typedef forward_iterator_tag iterator_category;
    typedef ptrdiff_t   difference_type;
  };

  #endif



protected:
  list<POINT> all_pnts;
  VECTOR      quasi_center;    // sum of the vertices of the origin simplex
  ch_simplex  origin_simplex;  // pointer to the origin simplex
  ch_facet    start_facet;     // pointer to some facet on the surface
  ch_vertex   anti_origin;

public: // until there are template friend functions possible
  IPOINT center() const        // compute center from quasi center
  { return CHTRAITS::to_ipoint(quasi_center/RT(dcur + 1)); }
  const VECTOR& qcenter() const
  { return quasi_center; }
  const ch_simplex& osimplex() const
  { return origin_simplex; }
  const PLANE& base_facet_plane(ch_simplex S) const
  { return S -> base_facet; }

protected:
  bool& visited_mark(ch_simplex S) const
  { return S -> visited; }

  void clean_dynamic_memory()
  {
    ch_simplex s;
    while (!all_simps.empty()) {
      s = (ch_simplex)all_simps.pop();
      delete(s);
    }
    ch_vertex v;
    while (!all_verts.empty()) {
      v = (ch_vertex)all_verts.pop();
      delete(v);
    }
  }

  int number_of_bounded_simplices;
  int number_of_unbounded_simplices;
  int number_of_visibility_tests;
  int number_of_points;
  int number_of_vertices;


  void  compute_equation_of_base_facet(ch_simplex s);
  /*{\Xop computes the equation of the base facet of $s$ and sets the
          |base_facet| member of $s$. The equation is normalized such
          that the origin simplex lies in the negative halfspace.}*/

  bool is_base_facet_visible(ch_simplex s, const POINT& x) const
  { return ( CHTRAITS::side(s->base_facet,x) > 0); }
  /*{\Xop returns true if $x$ sees the base facet of $s$, i.e., lies in
          its positive halfspace.}*/

  bool contains_in_base_facet(ch_simplex s, const POINT& x) const;
  /*{\Xop returns true if $x$ lies in the closure of the base facet of
          $s$.}*/


  ch_simplex  new_simplex()
  {
    ch_simplex s = new ch_Simplex<CHTRAITS,POINT,PLANE>(dmax);
    all_simps.append(s);
    return s;
  }
  /* creates a new simplex, adds it to the list of all simplices,
     and returns it. */



  void visibility_search(ch_simplex S, const POINT& x,
                         list<ch_simplex>& visible_simplices,
                         int& number_of_visited_simplices,
                         int& location,
                         ch_facet& f) const;
  /*{\Xop adds all unmarked unbounded simplices with $x$-visible base
          facet to |visible_simplices| and marks them. In |location| the
          procedure returns the position of |x| with respect to the
          current hull: $-1$ for inside, $0$ for on the the boundary,
          and $+1$ for outside; the initial value of |location| for the
          outermost call must be $-1$. If $x$ is contained in the
          boundary of |\Mvar| then a facet incident to $x$ is returned
          in $f$.}*/

  void clear_visited_marks(ch_Simplex<CHTRAITS,POINT,PLANE>* s) const;
  /*{\Xop removes the mark bits from all marked simplices reachable from $s$. }*/


  void dimension_jump(ch_simplex S, ch_vertex x);
  /*{\Xop Adds a new vertex $x$ to the triangulation. The point associated
  with $x$ lies outside the affine hull of the current point set.  }*/

public:

  /*{\Mcreation 3}*/

  chull(int d = 2);
  /*{\Mcreate creates an instance |\Mvar| of type |\Mtype|. The
  dimension of the underlying space is $d$ and |S| is initialized to the
  empty point set. The traits class |CHTRAITS| specifies the models of
  all types and the implementations of all geometric primitives used by
  the convex hull class.  In the following we use further template
  parameters like the point type |POINT=CHTRAITS::POINT| and the
  hyperplane type |PLANE=CHTRAITS::PLANE|.  At this point, it suffices
  to say that |POINT| represents points in $d$-space and that |PLANE|
  represents hyperplanes in $d$-space. The complete specification of the
  traits class is to be found in section \ref{chulltraitsclass}.}*/

  /*{\Mtext The above template instantiation works only for compilers
  which can handle \textit{template default arguments}. The exact
  template prototype is
  [[<CHTRAITS, POINT = CHTRAITS::POINT, PLANE = CHTRAITS::PLANE>]].
  In case you use a compiler without that feature please write
  [[chull<CHTRAITS,POINT,PLANE>]] where we write [[\Mname]] after
  providing [[POINT]] and [[PLANE]] in the global scope.  }*/

  ~chull() { clean_dynamic_memory(); }

  protected:

   /* the default copy constructor and assignment operator for class
      triangulation work incorrectly; it is therefore good practice to
      either implement them correctly or to make them inaccessible. We
      do the latter. */

    chull(const chull&);
    chull& operator=(const chull&);

  /*{\Mtext The data type |\Mtype| offers neither copy constructor nor
            assignment operator.}*/

  public:


  /*{\Moperations 3}*/
  /*{\Mtext All operations below that take a point |x| as argument have
            the common precondition that |x| is a point of ambient
            space.}*/

  /*{\Moptions nextwarning=no}*/
  /* inherited
  int dim(){ return dmax; }
  */
  /*{\Mop returns the dimension of ambient space}*/

  /*{\Moptions nextwarning=no}*/
  /* inherited
  int dcurrent(){ return dcur; }
  */
  /*{\Mop returns the affine dimension of $S$.}*/

  POINT  associated_point(ch_vertex v) const
  { return regl_complex<CHTRAITS,POINT>::associated_point(v);}
  /*{\Mop returns the point associated with vertex $v$.}*/

  ch_vertex  vertex_of_simplex(ch_simplex s, int i) const
  { return regl_complex<CHTRAITS,POINT>::vertex(s,i); }
  /*{\Mop returns the vertex corresponding to the $i$-th vertex of $s$.\\
          \precond $0 \leq i \leq |dcur|$. }*/

  POINT  point_of_simplex(ch_simplex s,int i) const
  { return associated_point(vertex_of_simplex(s,i)); }
  /*{\Mop same as |C.associated_point(C.vertex_of_simplex(s,i))|. }*/

  ch_simplex  opposite_simplex(ch_simplex s,int i) const
  { return (ch_simplex)regl_complex<CHTRAITS,POINT>::opposite_simplex(s,i); }
  /*{\Mop returns the simplex opposite to the $i$-th vertex of $s$
          (|nil| if there is no such simplex).\\
          \precond $0 \leq i \leq |dcur|$. }*/

  int  index_of_vertex_in_opposite_simplex(ch_simplex s,int i) const
  { return regl_complex<CHTRAITS,POINT>::index_of_opposite_vertex(s,i); }
  /*{\Mop returns the index of the vertex opposite to the $i$-th vertex
          of $s$.\\
          \precond $0 \leq i \leq |dcur|$ and there is a
          simplex opposite to the $i$-th vertex of $s$. }*/

  ch_simplex  simplex(ch_vertex v) const
  { return (ch_simplex)regl_complex<CHTRAITS,POINT>::simplex(v); }
  /*{\Mop returns a simplex of which $v$ is a node. Note that this
          simplex is not unique. }*/

  /*{\Moptions nextwarning=no}*/
  /* inherited
  int index(ch_vertex v){ return v -> i; }
  */
  /*{\Mop returns the index of $v$ in |simplex(v)|.}*/

  ch_vertex  vertex_of_facet(ch_facet f, int i) const
  { return vertex_of_simplex(f,i+1); }
  /*{\Mop returns the vertex corresponding to the $i$-th vertex of $f$.
  \precond $0 \leq i < |dcur|$. }*/

  POINT  point_of_facet(ch_facet f, int i) const
  { return point_of_simplex(f,i+1); }
  /*{\Mop same as |C.associated_point(C.vertex_of_facet(f,i))|.}*/

  ch_facet  opposite_facet(ch_facet f, int i) const
  { return opposite_simplex(f,i+1); }
  /*{\Mop returns the facet opposite to the $i$-th vertex of $f$ (|nil|
          if there is no such facet).\\
          \precond $0 \leq i < |dcur|$ and |dcur > 1|. }*/

  int  index_of_vertex_in_opposite_facet(ch_facet f, int i) const
  { return index_of_vertex_in_opposite_simplex(f,i+1) - 1; }
  /*{\Mop returns the index of the vertex opposite to the $i$-th vertex of $f$.\\
          \precond $0 \leq i < |dcur|$ and |dcur > 1|.}*/

  PLANE hyperplane_supporting(ch_facet f) const
  { return f -> base_facet; }
  /*{\Mop returns a hyperplane supporting facet |f|. The hyperplane is
          oriented such that the interior of |\Mvar| is on the negative
          side of it.\\
          \precond |f| is a facet of |\Mvar| and |dcur > 1|.}*/



  ch_vertex insert(const POINT& x);
  /*{\Mop adds point |x| to the underlying set of points.  If $x$ is
          equal to (the point associated with) a vertex of the current
          hull this vertex is returned and its associated point is
          changed to $x$. If $x$ lies outside the current hull, a new
          vertex |v| with associated point $x$ is added to the hull and
          returned. In all other cases, i.e., if $x$ lies in the
          interior of the hull or on the boundary but not on a vertex,
          the current hull is not changed and |nil| is returned.}*/


  bool is_dimension_jump(const POINT& x) const
  /*{\Mop returns true if $x$ is not contained in the affine hull of |S|. }*/
  {
    if (dcur == dmax) return false;
    array<POINT>  A(0,dcur);
    for (int l = 0; l <= dcur; l++)
      A[l] = point_of_simplex(origin_simplex,l);
    return ( !CHTRAITS::contained_in_affine_hull(A,x) );
  }


  list<ch_facet>  facets_visible_from(const POINT& x);
  /*{\Mop returns the list of all facets that are visible from |x|.\\
          \precond |x| is contained in the affine hull of |S|.}*/

  int is_where(const POINT& x);
  /*{\Mop returns |-1| (|0|,|1|) if |x| is contained in the interior
          (lies on the boundary, is contained in the exterior) of
          |\Mvar|.\\
          \precond |x| is contained in the affine hull of |S|.}*/



  bool is_unbounded_simplex(ch_simplex S) const
  { return vertex_of_simplex(S,0) == anti_origin; }

  bool is_bounded_simplex(ch_simplex S) const
  { return vertex_of_simplex(S,0) != anti_origin; }

  void init(int d = 2)
  {
    clean_dynamic_memory();
    quasi_center = CHTRAITS::zero_vector(d);
    anti_origin = nil;
    origin_simplex = nil;
    all_pnts.clear();
    regl_complex<CHTRAITS,POINT>::init(d);
    number_of_points = number_of_vertices = 0;
    number_of_unbounded_simplices = number_of_bounded_simplices = 0;
    number_of_visibility_tests = 0;
  }
  /*{\Mop reinitializes |\Mvar| to an empty hull in $d$-dimensional space.}*/

  #define STATISTIC(t) cout << #t << " = " << t << endl

  void print_statistics()
  {
    newline;
    cout << "chull(" << dim() << ") - statistic" << endl;
    STATISTIC(number_of_points);
    STATISTIC(number_of_vertices);
    STATISTIC(number_of_unbounded_simplices);
    STATISTIC(number_of_bounded_simplices);
    STATISTIC(number_of_visibility_tests);
    newline;
  }
  /*{\Mop gives information about the size of the current hull and the
          number of visibility tests performed. }*/

  #undef STATISTIC


  void check() const;
  /*{\Mop checks the validity of the data structure as described in
          \cite{Mehlhorn-et-al:Checking}. }*/


  /*{\Mtext \headerline{Lists and Iterators}
  \setopdims{3.5cm}{3cm}}*/

  list<POINT>  all_points() const
  { return all_pnts; }
  /*{\Mop returns a list of all points in |\Mvar|.}*/

  /*{\Moptions nextwarning=no}*/
  /*
  list<ch_vertex>  all_vertices() const;
  */
  /*{\Mop returns a list of all vertices of |\Mvar|
  (also interior ones).}*/

  list<ch_simplex> all_simplices() const
  { list_item it;
    list<ch_simplex> result;
    forall_items(it,all_simps) {
      result.append((ch_simplex) all_simps[it]);
    }
    return result;
  }
  /*{\Mop returns a list of all simplices in |\Mvar|.}*/

  list<ch_facet>  all_facets() const;
  /*{\Mop returns a list of all facets of |\Mvar|.}*/

  #ifdef DDGEO_STL_ITERATORS
  /*{\Mtext \textbf{Iterator} types and operations are only
  defined if compiled including [[-DDDGEO_STL_ITERATORS]].
  \setopdims{4cm}{3cm}}*/

  ch_point_iterator points_begin() const
  { return ch_point_iterator(this,all_pnts.first_item()); }
  /*{\Mop returns the start iterator for points in |\Mvar|.}*/

  ch_point_iterator points_end() const
  { return ch_point_iterator(this,
      all_pnts.stl_next_item(all_pnts.last_item())); }
  /*{\Mop returns the past the end iterator for points in |\Mvar|.}*/

  /*{\Moptions nextwarning=no}*/
  /*
  ch_vertex_iterator vertices_begin() const
  */
  /*{\Mop returns the start iterator for vertices of |\Mvar|.}*/

  /*{\Moptions nextwarning=no}*/
  /*
  ch_vertex_iterator vertices_end() const
  */
  /*{\Mop returns the past the end iterator for vertices of |\Mvar|.}*/


  ch_simplex_iterator simplices_begin() const
  { return ch_simplex_iterator(this,all_simps.first_item()); }
  /*{\Mop returns the start iterator for simplices of |\Mvar|.}*/

  ch_simplex_iterator simplices_end() const
  { return ch_simplex_iterator(this,
      all_simps.stl_next_item(all_simps.last_item())); }
  /*{\Mop returns the past the end iterator for simplices of |\Mvar|.}*/

  ch_facet_iterator facets_begin() const
  { return ch_facet_iterator(this,start_facet); }
  /*{\Mop returns the start iterator for simplices of |\Mvar|.}*/

  ch_facet_iterator facets_end() const
  { return ch_facet_iterator(this,nil); }
  /*{\Mop returns the past the end iterator for simplices of |\Mvar|.}*/
  #endif



};


/*{\Mimplementation The implementation of type |\Mtype| is based on
\cite{CMS} and \cite{BMS:degeneracy}.  The details of the
implementation can be found in the implementation document available
at the downloadsite of this package.

The time and space requirement is input dependent.  Let $C_1$, $C_2$,
$C_3$, \ldots be the sequence of hulls constructed and for a point $x$
let $k_i$ be the number of facets of $C_i$ that are visible from $x$
and that are not already facets of $C_{i-1}$. Then the time for
inserting $x$ is $O(|dim| \sum_i k_i)$ and the number of new simplices
constructed during the insertion of $x$ is the number of facets of the
hull which were not already facets of the hull before the insertion.

The data type |\Mtype| is derived from |regl_complex|. The space
requirement of regular complexes is essentially $12(|dim| +2)$ Bytes
times the number of simplices plus the space for the points. |\Mtype|
needs an additional $8 + (4 + x)|dim|$ Bytes per simplex where $x$ is
the space requirement of the underlying number type and an additional
$12$ Bytes per point. The total is therefore $(16 + x)|dim| + 32$
Bytes times the number of simplices plus $28 + x \cdot |dim|$ Bytes
times the number of points.}*/

/*{\Mexample See the section on Delaunay triangulations and convex
hulls of the LEDA book for an extensive example. Or the demo package
of the LEP |dd_geokernel|.}*/


/*{\Mtext \headerline{Low Dimensional Output Routines}
\setopdims{2cm}{3cm}}*/

template <class CHTRAITS,class POINT,class PLANE>
void d2_show(const chull<CHTRAITS,POINT,PLANE>& C, window& W);
/*{\Mfunc draws the convex hull |C| in window |W|.\\
\precond |dim == 2|. }*/

/*{\Moptions nextwarning=no}*/
/* inherited from regl_complex:
void d2_map(const chull<CHTRAITS,POINT,PLANE>& C, GRAPH<PNT,int>& G);
*/
/*{\Mfunc constructs the representation of |C| as a bidirected graph |G|.\\
\precond |dim == 3|.}*/

template <class CHTRAITS,class POINT,class PLANE>
void d3_surface_map(const chull<CHTRAITS,POINT,PLANE>& C,
                    GRAPH<POINT,int>& G);
/*{\Mfunc constructs the representation of the surface of |\Mvar| as a
bidirected graph |G|.\\ \precond |dim == 3|.}*/


#ifdef THIS_SHOULD_NEVER_BE_PARSED

/*{\Moptions
outfile=chulltraits.man
}*/

/*{\Manpage {CHTRAITS} {} {The Traits Class for Convex Hulls} {}}*/
/*{\Mtext \label{chulltraitsclass} }*/

/*{\Mdefinition The class |\Mname| allows the adaptation of the convex
hull class to an environment. It defines the models for the types |POINT|,
|PLANE|, |RT| (ring type), |VECTOR|, |IPOINT| (internal point type) and
|IRAY|, and the implementations for the geometric primitives used by the convex hull class.

In order to adapt the convex hull class to a new environment, the
class |\Mname| must be defined appropriately. We specify the
requirements for the class |\Mname| below. We also present three concrete
instantiations in the example section below:

\begin{itemize}
\item |d2_chull_traits| is an instantiation for convex hulls in
$2$-space based on LEDA's $2d$ geometry kernel.

\item |d3_chull_traits| is an instantiation for convex hulls in
$3$-space based on LEDA's $3d$ geometry kernel.

\item |dd_chull_traits| is an instantiation for convex hulls in
$d$-space based on LEDA's $dd$ geometry kernel (LEP |dd_geo_kernel|).
\end{itemize} }*/


class CHTRAITS {
/*{\Mcreation}*/
/*{\Mtext There are no constructors necessary as this type only provides
necessary operations by inline static operations.}*/

public:
/*{\Mtypes 4}*/
typedef user_point_type POINT;  // vertex type = chull<POINT,..>
/*{\Mtypedef instantiates the point class. POINT must be able to represent
points in $d$-space.}*/

typedef user_hyperplane_type PLANE;  // hyperplane type = chull<..PLANE..>
/*{\Mtypedef instantiates the hyperplane class. |PLANE| must be able
to represent all hyperplanes defined by points in |POINT|.}*/

/*{\Mtext The following types are used inside the class |chull|. They
are not visible to the outside}*/

typedef user_integer_type RT;
/*{\Mtypedef instantiates the arithmetic ring type used for calculations.
|RT| has to be a arbitrary precision integer.}*/

typedef user_vector_type VECTOR; // internal linear algebra
/*{\Mtypedef instantiates the internal vector type which is used for
the calculations of linear algebra. The type has to provide the
operations |v = v1 + v2|, |v += v1|, |v = v1/(RT)|.}*/

typedef user_internal_point_type IPOINT; // internal calculation
/*{\Mtypedef instantiates the point class of chull which is used for
internal calculations.  IPOINT must have the same dimension as
chull.}*/

typedef user_internal_ray_type IRAY;   // internal calculation
/*{\Mtypedef instantiates the ray class of chull which is used for
internal calculations in the affine space.  IRAY must have the same
dimension as chull. We expect IRAY to have a constructor |IRAY(p1,p2)|
where |p1| and |p2| are two IPOINTs.}*/

/*{\Mtext{ The reader may wonder why we distinguish between classes
|POINT| and |IPOINT|.  The inputs come from class |POINT|. The
algorithm will never construct an object of class |POINT|. Thus
|POINT| must only be able to represent the actual input points.  On
the other hand the points in |IPOINT| are constructed by the algorithm
and hence |IPOINT| must be able to represent an arbitrary point in
$d$-space.  }}*/

/*{\Moperations 4}*/

static int side(const PLANE& h, const POINT& p);
/*{\Mstatic computes the side of |h| on which |p| lies.}*/

static bool affinely_independent(const array<POINT>& A);
/*{\Mstatic decides whether the points in $A$ are affinely independent.}*/

static bool contained_in_simplex(const array<POINT>& A,
                                 const POINT& p);
/*{\Mstatic determines whether $p$ is contained in the simplex spanned by
the points in |A|.}*/

static bool contained_in_affine_hull(const array<POINT>& A,
                                     const POINT& p);
/*{\Mstatic determines whether $p$ is contained in the affine hull of the
points in $A$.}*/

static VECTOR normal(const PLANE& h);
/*{\Mstatic returns the normal vector of |h|. It points from the negative
halfspace into the positive halfspace.}*/

static VECTOR zero_vector(int d);
/*{\Mstatic creates a vector of dimension $d$ with zero entries.}*/

static VECTOR to_vector(const POINT& p);
/*{\Mstatic converts the point to a vector.}*/

static IPOINT to_ipoint(const VECTOR& p);
/*{\Mstatic converts the vector to an internal point.}*/

static PLANE hyperplane_through_points(
             const array<POINT>& A, const IPOINT& p, int side);
/*{\Mstatic constructs some hyperplane that passes through the points
in $A$.  If |side| is positive or negative then $p$ is on that side of
the constructed hyperplane.\\ \precond A hyperplane with the stated
properties must exist.}*/


static bool intersection(const PLANE& h, const IRAY& r,
                         IPOINT& p);
/*{\Mstatic returns true if the intersection set $S = h \cap r$ is
empty, else false. |p| is a point on |h| if $S \neq \emptyset$.}*/

static bool contained_in_simplex(const array<POINT>& A,
                                   const IPOINT& p);
/*{\Mstatic determines whether $p$ is contained in the simplex spanned
by the points in |A|. Note that there are two equally named
operations.  You can ommit the second if |POINT == IPOINT|.}*/

static int side(const PLANE& h, const IPOINT& p);
/*{\Mstatic computes the side of |h| on which |p| lies. Note that
there are two equally named operations. You can ommit the second if
|POINT == IPOINT|.}*/

/*{\Mtext
\headerline{Optional - Low-dimensional conversion to LEDA modules}
The following Operations are needed if the visualization and graph
modules of LEDA are used. Users can omit them if they don't use
the corresponding operations of |chull|. }*/

static int orientation(const POINT& p1, const POINT& p2,
                       const POINT& p3, const IPOINT& p4);
/*{\Mstatic determines the orientation of the points $p1,p2,p3,p4$
in $3$-space. Implementation \emph{only necessary} if |d2_surface_map| is
used.}*/

static int orientation(const POINT& p1, const POINT& p2,
                       const POINT& p3);
/*{\Mstatic determines the orientation of the points $p1,p2,p3$
in $2$-space. Implementation \emph{only necessary} if |d2_map| is used.}*/

  static point to_d2_point(const POINT& p);
/*{\Mstatic converts the point to a $2d$ LEDA point. Implementation
\emph{only necessary} if |d2_show| is used.}*/

};

/*{\Mexample
Together with the abstract data type |chull| come three instantiations.
\begin{enumerate}
\item $2$-dimensional convex hulls

The class |d2_chull_traits| is an instantiation for $2d$-convex hulls
based on LEDA's $2d$ geometry kernel. We have
|POINT| $\equiv$ |LEDA/rat_point| and |PLANE| $\equiv$ |LEDA/rat_line|.
In order to use |chull| together with |d2_chull_traits| and
your compiler does know {\bf not} how to handle template default parameters
write

\begin{Mverb}
#include <LEP/dd_geo/chull.h>
#include <LEP/dd_geo/d2_chull_traits.h>

typedef d2_chull_traits D2CHTRAITS;
typedef D2CHTRAITS::POINT D2POINT;
typedef D2CHTRAITS::PLANE D2PLANE;
typedef chull<D2CHTRAITS,D2POINT,D2PLANE> d2_chull;
typedef d2_chull::ch_vertex ch2_vertex;
typedef d2_chull::ch_simplex ch2_simplex;
typedef d2_chull::ch_facet ch2_facet;

main()
{
  d2_chull T(2);
  ch2_vertex v1 = T.insert(D2POINT(2,11));
  ...
}
\end{Mverb}

\item $3$-dimensional convex hulls

The class |d3_chull_traits| is an instantiation for $3d$-convex hulls
based on LEDA's $3d$ geometry kernel. We have |POINT| $\equiv$
|LEDA/d3_rat_point| and |PLANE| $\equiv$ |LEDA/d3_rat_plane|.  In
order to use |chull| together with |d3_chull_traits| and your compiler
does know how to handle template default parameters write

\begin{Mverb}
#include <LEP/dd_geo/chull.h>
#include <LEP/dd_geo/d3_chull_traits.h>

typedef d3_chull_traits D3CHTRAITS;
typedef D3CHTRAITS::POINT D3POINT;
typedef D3CHTRAITS::PLANE D3PLANE;
typedef chull<D3CHTRAITS> d3_chull;
typedef d3_chull::ch3_vertex ch_vertex;
typedef d3_chull::ch3_simplex ch_simplex;
typedef d3_chull::ch3_facet ch_facet;

main()
{
  d3_chull T(3);
  ch3_vertex v1 = T.insert(D3POINT(2,5,11));
  ...
}
\end{Mverb}

\item $d$-dimensional convex hulls

The class |dd_chull_traits| is an instantiation for $dd$-convex hulls
based on the LEP dd\_geokernel's geometry kernel. We have |POINT|
$\equiv$ |LEP/dd_geo/hpoint<integer>| and |PLANE| $\equiv$
|LEP/dd_geo/hhyperplane<integer>|.  In order to use |chull| together
with |dd_chull_traits| and your compiler does {\bf not} know how to
handle template default parameters write

\begin{Mverb}
#include <LEP/dd_geo/chull.h>
#include <LEP/dd_geo/dd_chull_traits.h>

typedef dd_chull_traits DDCHTRAITS;
typedef DDCHTRAITS::POINT DDPOINT;
typedef DDCHTRAITS::PLANE DDPLANE;
typedef chull<DDCHTRAITS,DDPOINT,DDPLANE> dd_chull;
typedef dd_chull::chd_vertex ch_vertex;
typedef dd_chull::chd_simplex ch_simplex;
typedef dd_chull::chd_facet ch_facet;

main()
{
  dd_chull T(5);
  chd_vertex v1 = T.insert(DDPOINT(...));
  ...
}
\end{Mverb}

\end{enumerate}

For the concrete implementation of the traits classes please refer
to the include files of the LEP. For more extensive examples please
have a look into the demo directory of the LEP.
}*/

#endif // THIS_SHOULD_NEVER_BE_PARSED



#include <CGAL/dd_geo/chull.C>


#if LEP_DDGEO_INCL_ID == 21010
#undef LEDA_ROOT_INCL_ID
#undef LEP_DDGEO_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif

#endif // CGAL_DD_GEO_CHULL_H


