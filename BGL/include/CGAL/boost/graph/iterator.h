// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BGL_ITERATORS_H
#define CGAL_BGL_ITERATORS_H

#include <stdexcept>

#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <CGAL/Iterator_range.h>
#include <CGAL/assertions.h>
#include <CGAL/circulator_bases.h>
#include <CGAL/boost/graph/internal/helpers.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

namespace internal {

template <typename G>
struct Edge {
  const G* g; 

  Edge()
    : g(NULL)
  {}

  Edge(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::edge_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return edge(h, *g);
  }
};

template <typename G>
struct Opposite_edge {
  const G* g; 

  Opposite_edge()
    : g(NULL)
  {}

  Opposite_edge(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::edge_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return edge(opposite(h,*g), *g);
  }
};

template <typename G>
struct Opposite_halfedge {
  const G* g; 

  Opposite_halfedge()
    : g(NULL)
  {}

  Opposite_halfedge(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::halfedge_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return opposite(h,*g);
  }
};

template <typename G>
struct Target {
  const G* g; 

  Target()
    : g(NULL)
  {}

  Target(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::vertex_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return target(h,*g);
  }
};

template <typename G>
struct Source {
  const G* g; 

  Source()
    : g(NULL)
  {}

  Source(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::vertex_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return source(h,*g);
  }
};

template <typename G>
struct Face {
  const G* g; 

  Face()
    : g(NULL)
  {}

  Face(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::face_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return face(h,*g);
  }
};
template <typename G>
struct Opposite_face {
  const G* g; 

  Opposite_face()
    : g(NULL)
  {}

  Opposite_face(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::face_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return face(opposite(h,*g),*g);
  }
};
} // namespace internal
/// \endcond

/**
 * \ingroup PkgBGLIterators
 * A bidirectional iterator with value type `boost::graph_traits<Graph>::%halfedge_descriptor` over all halfedges having the same vertex as source.
 * Let `h` be a halfedge of graph `g`. For a `Halfedge_around_source_iterator` `havi` with `h = *havi;` 
 * the following holds: Either `++havi` is the past the end iterator, or `next(opposite(h,g),g) == *++havi`. 
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalIterator`
 */
template <typename Graph>
class Halfedge_around_source_iterator {
  typedef Halfedge_around_source_iterator Self;
public:
#ifndef DOXYGEN_RUNNING
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor vertex_descriptor;
  typedef std::bidirectional_iterator_tag     iterator_category;
  typedef halfedge_descriptor                 value_type;
  typedef value_type*                         pointer;
  typedef const value_type&                   reference;
  typedef std::ptrdiff_t                      difference_type;

private:  
  halfedge_descriptor anchor, pos;
  const Graph* g;
  int      winding;
#endif 
  
public:
  Halfedge_around_source_iterator()
    : anchor(), pos(), g(0)
  {}

  Halfedge_around_source_iterator(halfedge_descriptor hd, const Graph& g, int n=0)
    : anchor(hd), pos(hd), g(&g), winding((hd==halfedge_descriptor())?1:n)
  {}

#ifndef DOXYGEN_RUNNING
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Halfedge_around_source_iterator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base() == NULL)) ?
      &Halfedge_around_source_iterator::this_type_does_not_support_comparisons : 0;
  }

  bool operator==( const Self& i) const {
    CGAL_assertion( anchor == anchor);
    return  ( g == i.g) && ( pos == i.pos) && ( winding == i.winding);
  }
  
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  bool operator== (void* ) const
  {
    return g == NULL;
  }

  reference operator*() const
  {
    return pos;
  }

  pointer operator->() const
  {
    return &pos;
  }

  Self& operator++() {
    pos = next(opposite(pos,*g),*g);
    if ( pos == anchor)
      ++winding;
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
  Self& operator--() {
    if ( pos == anchor)
      --winding;
    pos = opposite(prev(pos,*g),*g);
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }
#endif
};

/**
 * \ingroup PkgBGLIterators
 * A bidirectional iterator with value type `boost::graph_traits<Graph>::%halfedge_descriptor` over all halfedges having the same vertex as target.
 * Let `h` be a halfedge of graph `g`. For a `Halfedge_around_target_iterator` `havi` with `h = *havi;` 
 * the following holds: Either `++havi` is the past the end iterator, or `opposite(next(h,g),g) == *++havi`. 
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalIterator`
 */

template <typename Graph>
class Halfedge_around_target_iterator {
  typedef Halfedge_around_target_iterator Self;

public:
#ifndef DOXYGEN_RUNNING
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor vertex_descriptor;
  typedef std::bidirectional_iterator_tag     iterator_category;
  typedef halfedge_descriptor                 value_type;
  typedef value_type*                         pointer;
  typedef const value_type&                   reference;
  typedef std::ptrdiff_t                      difference_type;

private:  
  halfedge_descriptor anchor, pos;
  const Graph* g;
  int      winding;
#endif 
  
public:
  Halfedge_around_target_iterator()
    : anchor(), pos(), g(0)
  {}

  Halfedge_around_target_iterator(halfedge_descriptor hd, const Graph& g, int n=0)
    : anchor(hd), pos(hd), g(&g), winding((hd==halfedge_descriptor())?1:n)
  {}

#ifndef DOXYGEN_RUNNING
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Halfedge_around_target_iterator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base() == NULL)) ?
      &Halfedge_around_target_iterator::this_type_does_not_support_comparisons : 0;
  }

  bool operator==( const Self& i) const {
    CGAL_assertion( anchor == anchor);
    return  ( g == i.g) && ( pos == i.pos) && ( winding == i.winding);
  }
  
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  bool operator== (void* ) const
  {
    return g == NULL;
  }

  reference operator*() const
  {
    return pos;
  }

  pointer operator->() const
  {
    return &pos;
  }

  Self& operator++() {
    pos = opposite(next(pos,*g),*g);
    if ( pos == anchor)
      ++winding;
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
  Self& operator--() {
    if ( pos == anchor)
      --winding;
    pos = prev(opposite(pos,*g),*g);
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }
#endif
};


/**
 * \ingroup PkgBGLIterators
 * A bidirectional iterator  with value type `boost::graph_traits<Graph>::%halfedge_descriptor` over all halfedges incident to the same face or border.
 * Let `h` be a halfedge of graph `g`. For a `Halfedge_around_face_iterator` `hafi` with  `h = *hafi` 
 * the following holds: Either `++hafi` is the past the end iterator, or `next(h,g) == *++hafi`.
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalIterator`
 */

template <typename Graph>
class Halfedge_around_face_iterator {
#ifndef DOXYGEN_RUNNING 
  typedef Halfedge_around_face_iterator Self;

public:
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef std::bidirectional_iterator_tag      iterator_category;
  typedef halfedge_descriptor                 value_type;
  typedef value_type*                         pointer;
  typedef value_type&                         reference;
  typedef std::ptrdiff_t                      difference_type;
  
#endif
  
  Halfedge_around_face_iterator()
    : pos(), g(0)
  {}

  Halfedge_around_face_iterator(halfedge_descriptor hd, const Graph& g, int n=0)
    : anchor(hd), pos(hd), g(&g), winding((hd == halfedge_descriptor())?1:n)
  {}
#ifndef DOXYGEN_RUNNING
  reference         operator *  ( )       { return  pos; }
  const value_type& operator *  ( ) const { return  pos; }
  pointer           operator -> ( )       { return &pos; }
  const value_type* operator -> ( ) const { return &pos; }

  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Halfedge_around_face_iterator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base() == NULL)) ?
      &Halfedge_around_face_iterator::this_type_does_not_support_comparisons : 0;
  }

  bool operator==( const Self& i) const {
    CGAL_assertion( anchor == anchor);
    return  ( g == i.g) && ( pos == i.pos) && ( winding == i.winding);
  }

  bool operator!=( const Self& i) const {
    return !(*this == i);
  }

  Self& operator++() 
  {
    CGAL_assertion(g != NULL);
    pos = next(pos,*g);
    if ( pos == anchor)
      ++winding;
    return *this;
  }

  Self operator++(int) 
  {
    CGAL_assertion(g != NULL);
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  Self& operator--() 
  {
    CGAL_assertion(g != NULL); 
    if ( pos == anchor)
      --winding;
  
    pos = prev(pos,*g);
    return *this;
  }

  Self operator--(int) 
  {
    CGAL_assertion(g != NULL);
    Self tmp = *this;
    --*this;
    return tmp;
  }
#endif
private:
  halfedge_descriptor anchor, pos;
  const Graph* g;
  int      winding;
};


template <typename Graph>
class Halfedge_around_target_circulator;

/**
 * \ingroup PkgBGLIterators
 * A bidirectional circulator  with value type `boost::graph_traits<Graph>::%halfedge_descriptor` over all halfedges having the same vertex as source.
 * Let `h` be a halfedge of graph `g`. For a `Halfedge_around_source_circulator` `havc` with `h = *havc;` 
 * the following holds: `next(opposite(h,g),g) == *++havc`. 
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalCirculator`
 */
template <typename Graph>
class Halfedge_around_source_circulator
#ifndef DOXYGEN_RUNNING
 : public boost::iterator_adaptor<
             Halfedge_around_source_circulator<Graph>                    // Derived
             , Halfedge_around_target_circulator<Graph>                  // Base
             , typename boost::graph_traits<Graph>::halfedge_descriptor  // Value
             , Bidirectional_circulator_tag                              // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::halfedge_descriptor  // Reference
             >
#endif
{
private:
  internal::Opposite_halfedge<Graph> opp;

#ifndef DOXYGEN_RUNNING
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
#endif

public:

#ifndef DOXYGEN_RUNNING
  typedef std::size_t size_type;
#endif

  Halfedge_around_source_circulator()
  {}

  Halfedge_around_source_circulator(halfedge_descriptor hd, const Graph& g)
    : Halfedge_around_source_circulator::iterator_adaptor_(Halfedge_around_target_circulator<Graph>(opposite(hd,g),g)), opp(g)
  {}

  Halfedge_around_source_circulator(vertex_descriptor vd, const Graph& g)
    : Halfedge_around_source_circulator::iterator_adaptor_(Halfedge_around_target_circulator<Graph>(halfedge(vd,g),g)), opp(g)
  {}
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Halfedge_around_source_circulator::*bool_type)() const;
  
  void this_type_does_not_support_comparisons() const {}
  
  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &Halfedge_around_source_circulator::this_type_does_not_support_comparisons : 0;
  }
  
  bool operator== (void*) const
  {
    return this->base_reference() == NULL;
  }
  
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::halfedge_descriptor dereference() const { return opp(*this->base_reference()); }
}; 


/**
 * \ingroup PkgBGLIterators
 * A bidirectional circulator  with value type `boost::graph_traits<Graph>::%face_descriptor` over all faces incident to the same vertex.
 * It circulates over the same halfedges as the `Halfedge_around_target_circulator`.
 *
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalCirculator`
 */

template <typename Graph>
class Face_around_target_circulator
#ifndef DOXYGEN_RUNNING
  : public boost::iterator_adaptor<
             Face_around_target_circulator<Graph>                    // Derived
             , Halfedge_around_target_circulator<Graph>                  // Base
             , typename boost::graph_traits<Graph>::face_descriptor  // Value
             , Bidirectional_circulator_tag                              // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::face_descriptor  // Reference
             >
#endif
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  internal::Face<Graph> fct;

public:

  Face_around_target_circulator()
  {}

  Face_around_target_circulator(halfedge_descriptor hd, const Graph& g)
    : Face_around_target_circulator::iterator_adaptor_(Halfedge_around_target_circulator<Graph>(hd,g)), fct(g)
  {}
#ifndef DOXYGEN_RUNNING
  typedef std::size_t size_type;

  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Face_around_target_circulator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &Face_around_target_circulator::this_type_does_not_support_comparisons : 0;
  }

  bool operator== (void*) const
  {
    return this->base_reference() == NULL;
  }


private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::face_descriptor dereference() const { return fct(*this->base_reference()); }
#endif

}; 


/**
 * \ingroup PkgBGLIterators
 * A bidirectional circulator  with value type `boost::graph_traits<Graph>::%halfedge_descriptor` over all halfedges having the same vertex as target.
 * Let `h` be a halfedge of graph `g`. For a `Halfedge_around_target_circulator` `havc` with `h = *havc;` 
 * the following holds: `opposite(next(h,g),g) == *++havc`. 
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalCirculator`
 */



template <typename Graph>
class Halfedge_around_target_circulator {
  typedef Halfedge_around_target_circulator Self;

public:
#ifndef DOXYGEN_RUNNING
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef Bidirectional_circulator_tag        iterator_category;
  typedef halfedge_descriptor                 value_type;
  typedef value_type*                         pointer;
  typedef value_type&                         reference;
  typedef std::ptrdiff_t                      difference_type;
  typedef std::size_t                         size_type;
#endif 

  Halfedge_around_target_circulator()
    : g(0)
  {}

  Halfedge_around_target_circulator(halfedge_descriptor pos, const Graph& g)
    : pos(pos), g(&g)
  {}

  Halfedge_around_target_circulator(vertex_descriptor vd, const Graph& g)
    : pos(halfedge(vd,g)), g(&g)
  {}

#ifndef DOXYGEN_RUNNING
  reference         operator *  ( )       { return  pos; }
  const value_type& operator *  ( ) const { return  pos; }
  pointer           operator -> ( )       { return &pos; }
  const value_type* operator -> ( ) const { return &pos; }

  bool       operator == ( const Self& other) const { return  g == other.g && pos == other.pos; }
  bool       operator != ( const Self& other) const { return  g != other.g || pos != other.pos; }


  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Halfedge_around_target_circulator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (g == NULL)) ?
      &Halfedge_around_target_circulator::this_type_does_not_support_comparisons : 0;
  }


  bool operator== (void* ) const
  {
    return g == NULL;
  }

 
  Self& operator++() 
  {
    CGAL_assertion(g != NULL);
    pos = opposite(next(pos,*g),*g);
    return *this;
  }

  Self operator++(int) 
  {
    CGAL_assertion(g != NULL);
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  Self& operator--() 
  {
    CGAL_assertion(g != NULL);
    pos = prev(opposite(pos,*g),*g);
    return *this;
  }

  Self operator--(int) 
  {
    CGAL_assertion(g != NULL);
    Self tmp = *this;
    --*this;
    return tmp;
  }
#endif

private:
  halfedge_descriptor pos;
  const Graph* g;
};


/**
 * \ingroup PkgBGLIterators
 * A bidirectional circulator  with value type `boost::graph_traits<Graph>::%halfedge_descriptor` over all halfedges incident to the same face or border.
 * Let `h` be a halfedge of graph `g`. For a `Halfedge_around_face_circulator` `hafc` with  `h = *hafc` 
 * the following holds:  `next(h,g) == *++hafc`.
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalCirculator`
 */

template <typename Graph>
class Halfedge_around_face_circulator {
  typedef Halfedge_around_face_circulator Self;

public:
#ifndef DOXYGEN_RUNNING
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef Bidirectional_circulator_tag        iterator_category;
  typedef halfedge_descriptor                 value_type;
  typedef value_type*                         pointer;
  typedef value_type&                         reference;
  typedef std::ptrdiff_t                      difference_type;
  typedef std::size_t                         size_type;
#endif 
  
  Halfedge_around_face_circulator()
    : pos(), g(0)
  {}

  Halfedge_around_face_circulator(halfedge_descriptor pos, const Graph& g)
    : pos(pos), g(&g)
  {}
#ifndef DOXYGEN_RUNNING
  reference         operator *  ( )       { return  pos; }
  const value_type& operator *  ( ) const { return  pos; }
  pointer           operator -> ( )       { return &pos; }
  const value_type* operator -> ( ) const { return &pos; }

  bool       operator == ( const Self& other) const { return  g == other.g && pos == other.pos; }
  bool       operator != ( const Self& other) const { return  g != other.g || pos != other.pos; }
 

  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Halfedge_around_face_circulator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (g == NULL)) ?
      &Halfedge_around_face_circulator::this_type_does_not_support_comparisons : 0;
  }

  bool operator== (void* ) const
  {
    return g == NULL;
  }

  Self& operator++() 
  {
    CGAL_assertion(g != NULL);
    pos = next(pos,*g);
    return *this;
  }

  Self operator++(int) 
  {
    CGAL_assertion(g != NULL);
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  Self& operator--() 
  {
    CGAL_assertion(g != NULL);
    pos = prev(pos,*g);
    return *this;
  }

  Self operator--(int) 
  {
    CGAL_assertion(g != NULL);
    Self tmp = *this;
    --*this;
    return tmp;
  }
#endif

private:
  halfedge_descriptor pos;
  const Graph* g;
};


/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all halfedges with vertex `source(h,g)` as source.
 */
template<typename Graph>
Iterator_range<Halfedge_around_source_iterator<Graph> >
halfedges_around_source(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Halfedge_around_source_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}

/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all halfedges with vertex `v` as source.
 */
template<typename Graph>
Iterator_range<Halfedge_around_source_iterator<Graph> >
halfedges_around_source(typename boost::graph_traits<Graph>::vertex_descriptor v, const Graph& g)
{
  return halfedges_around_source(opposite(halfedge(v,g),g),g);
}

/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all halfedges with vertex `target(h,g)` as target. 
 */
template<typename Graph>
Iterator_range<Halfedge_around_target_iterator<Graph> >
halfedges_around_target(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Halfedge_around_target_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}

/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all halfedges with vertex `v` as target. 
 */
template<typename Graph>
Iterator_range<Halfedge_around_target_iterator<Graph> >
halfedges_around_target(typename boost::graph_traits<Graph>::vertex_descriptor v, const Graph& g)
{
  return halfedges_around_target(halfedge(v,g),g);
}

/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all halfedges incident to the same face or border as `h`. 
 */
template<typename Graph>
Iterator_range<Halfedge_around_face_iterator<Graph> >
halfedges_around_face(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Halfedge_around_face_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}

/**
 * \ingroup PkgBGLIterators
 * A bidirectional iterator with value type `boost::graph_traits<Graph>::%face_descriptor`.
 * It iterates over the same halfedges as the `Halfedge_around_face_iterator`,
 * and provides the face descriptor associated to the opposite halfedge.  The face descriptor
 * may be the null face, and it may be several times the same face descriptor.
 *
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalCirculator`
 */
template <typename Graph>
class Face_around_face_iterator
#ifndef DOXYGEN_RUNNING
  : public boost::iterator_adaptor<
            Face_around_face_iterator<Graph>                       // Derived
             , Halfedge_around_face_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::face_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::face_descriptor  // Reference
             >
#endif
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  internal::Opposite_face<Graph> fct;
public:

  Face_around_face_iterator()
  {}

  Face_around_face_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : Face_around_face_iterator::iterator_adaptor_(Halfedge_around_face_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), fct(g)
  {}
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::face_descriptor dereference() const { return fct(*this->base_reference()); }
}; 


/**
 * \ingroup PkgBGLIterators
 * A bidirectional circulator  with value type `boost::graph_traits<Graph>::%face_descriptor`.
 * It circulates over the same halfedges as the `Halfedge_around_face_circulator`,
 * and provides the face descriptor associated to the opposite halfedge.  The face descriptor
 * may be the null face, and it may be several times the same face descriptor.
 *
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalCirculator`
 */
template <typename Graph>
class Face_around_face_circulator
{};

/**
 * \ingroup PkgBGLIterators
 * A bidirectional iterator with value type `boost::graph_traits<Graph>::%face_descriptor`.
 * It iterates over the same halfedges as the `Halfedge_around_face_iterator`,
 * and provides the face descriptor associated to the opposite halfedge.  The face descriptor
 * may be the null face, and it may be several times the same face descriptor.
 *
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalIterator`
 */
template <typename Graph>
class Face_around_target_iterator
#ifndef DOXYGEN_RUNNING  
  : public boost::iterator_adaptor<
            Face_around_target_iterator<Graph>                       // Derived
             , Halfedge_around_target_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::face_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::face_descriptor  // Reference
             >
#endif
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  internal::Face<Graph> fct;

public:

  Face_around_target_iterator()
  {}

  Face_around_target_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : Face_around_target_iterator::iterator_adaptor_(Halfedge_around_target_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), fct(g)
  {}

private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::face_descriptor dereference() const { return fct(*this->base_reference()); }
}; 

/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all faces around  vertex `target(h,g)`. 
 */
template<typename Graph>
Iterator_range<Face_around_target_iterator<Graph> >
faces_around_target(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Face_around_target_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}

/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all edge-adjacent faces to the same face `face(h,g)`.
 */
template<typename Graph>
Iterator_range<Face_around_face_iterator<Graph> >
faces_around_face(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Face_around_face_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}

template <typename Graph>
class Vertex_around_face_circulator
#ifndef DOXYGEN_RUNNING 
  : public boost::iterator_adaptor<
             Vertex_around_face_circulator<Graph>                    // Derived
             , Halfedge_around_face_circulator<Graph>                  // Base
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Value
             , Bidirectional_circulator_tag                              // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Reference
             >
#endif
{
  internal::Target<Graph> fct;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

public:
#ifndef DOXYGEN_RUNNING
  typedef std::size_t size_type;
#endif

  Vertex_around_face_circulator()
  {}

  Vertex_around_face_circulator(halfedge_descriptor h, const Graph& g)
    : Vertex_around_face_circulator::iterator_adaptor_(Halfedge_around_face_circulator<Graph>(h,g)), fct(g)
  {}

#ifndef DOXYGEN_RUNNING
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Vertex_around_face_circulator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &Vertex_around_face_circulator::this_type_does_not_support_comparisons : 0;
  }

  bool operator== (void*) const
  {
    return this->base_reference()== NULL;
  }
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::vertex_descriptor dereference() const { return fct(*this->base_reference()); }
#endif
}; 

/**
 * \ingroup PkgBGLIterators
 * A bidirectional iterator  with value type `boost::graph_traits<Graph>::%vertex_descriptor`
 *  over all vertices incident to the same face or border.
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalIterator`
 */
template <typename Graph>
class Vertex_around_face_iterator
#ifndef DOXYGEN_RUNNING
   : public boost::iterator_adaptor<
            Vertex_around_face_iterator<Graph>                       // Derived
             , Halfedge_around_face_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Reference
             >
#endif
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  internal::Target<Graph> fct;
public:

  Vertex_around_face_iterator()
  {}

  Vertex_around_face_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : Vertex_around_face_iterator::iterator_adaptor_(Halfedge_around_face_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), fct(g)
  {}

#ifndef DOXYGEN_RUNNING
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Vertex_around_face_iterator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &Vertex_around_face_iterator::this_type_does_not_support_comparisons : 0;
  }

  bool operator== (void*) const
  {
    return this->base_reference()== NULL;
  }
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::vertex_descriptor dereference() const { return fct(*this->base_reference()); }
#endif
}; 


template <typename Graph>
class Opposite_edge_around_face_iterator
#ifndef DOXYGEN_RUNNING
  : public boost::iterator_adaptor<
            Opposite_edge_around_face_iterator<Graph>                       // Derived
             , Halfedge_around_face_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::edge_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::edge_descriptor  // Reference
             >
#endif
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  internal::Opposite_edge<Graph> fct;
public:

  Opposite_edge_around_face_iterator()
  {}

  Opposite_edge_around_face_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : Opposite_edge_around_face_iterator::iterator_adaptor_(Halfedge_around_face_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), fct(g)
  {}
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::edge_descriptor dereference() const { return fct(*this->base_reference()); }
}; 

template<typename Graph>
Iterator_range<Opposite_edge_around_face_iterator<Graph> >
opposite_edges_around_face(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Opposite_edge_around_face_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}

template <typename Graph>
class Edge_around_face_iterator
#ifndef DOXYGEN_RUNNING
  : public boost::iterator_adaptor<
            Edge_around_face_iterator<Graph>                       // Derived
             , Halfedge_around_face_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::edge_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::edge_descriptor  // Reference
             >
#endif
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  internal::Edge<Graph> fct;
public:

  Edge_around_face_iterator()
  {}

  Edge_around_face_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : Edge_around_face_iterator::iterator_adaptor_(Halfedge_around_face_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), fct(g)
  {}
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::edge_descriptor dereference() const { return fct(*this->base_reference()); }
}; 

template<typename Graph>
Iterator_range<Edge_around_face_iterator<Graph> >
edges_around_face(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Edge_around_face_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}



/**
 * \ingroup PkgBGLIterators
 * A bidirectional circulator with value type `boost::graph_traits<Graph>::%vertex_descriptor` over all vertices adjacent to the same vertex.
 * It circulates over the same halfedges as the `Halfedge_around_target_circulator`.
 *
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalCirculator`
 */
template <typename Graph>
class Vertex_around_target_circulator
#ifndef DOXYGEN_RUNNING
  : public boost::iterator_adaptor<
            Vertex_around_target_circulator<Graph>                       // Derived
             , Halfedge_around_target_circulator<Graph>                // Base
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Value
             , Bidirectional_circulator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Reference
             >
#endif
{
  internal::Source<Graph> fct;
public:
#ifndef DOXYGEN_RUNNING
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef std::size_t size_type;
#endif

  Vertex_around_target_circulator()
  {}

  Vertex_around_target_circulator(halfedge_descriptor h, const Graph& g)
    : Vertex_around_target_circulator::iterator_adaptor_(Halfedge_around_target_circulator<Graph>(h,g)), fct(g)
  {}

#ifndef DOXYGEN_RUNNING  
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Vertex_around_target_circulator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &Vertex_around_target_circulator::this_type_does_not_support_comparisons : 0;
  }

  bool operator== (void*) const
  {
    return this->base_reference()== NULL;
  }

private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::vertex_descriptor dereference() const { return fct(*this->base_reference()); }
#endif
}; 



/**
 * \ingroup PkgBGLIterators
 * A bidirectional iterator with value type `boost::graph_traits<Graph>::%vertex_descriptor` over all vertices adjacent to the same vertex.
 * It iterates over the same halfedges as the `Halfedge_around_target_iterator`.
 *
 * \tparam Graph must be a model of the concept `HalfedgeGraph`
 * \cgalModels `BidirectionalIterator`
 */
template <typename Graph>
class Vertex_around_target_iterator
#ifndef DOXYGEN_RUNNING
  : public boost::iterator_adaptor<
            Vertex_around_target_iterator<Graph>                       // Derived
             , Halfedge_around_target_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::vertex_descriptor  // Reference
             >
#endif
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  internal::Source<Graph> fct;

public:

  Vertex_around_target_iterator()
  {}

  Vertex_around_target_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : Vertex_around_target_iterator::iterator_adaptor_(Halfedge_around_target_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), fct(g)
  {}

#ifndef DOXYGEN_RUNNING
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Vertex_around_target_iterator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &Vertex_around_target_iterator::this_type_does_not_support_comparisons : 0;
  }
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::vertex_descriptor dereference() const { return fct(*this->base_reference()); }
#endif
}; 


template <typename Graph>
Iterator_range<Vertex_around_target_iterator<Graph> >
adjacent_vertices(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Vertex_around_face_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}


template <typename Graph>
Iterator_range<Vertex_around_target_iterator<Graph> >
adjacent_vertices(typename boost::graph_traits<Graph>::vertex_descriptor v, const Graph& g)
{
  typedef Vertex_around_face_iterator<Graph> I;
  return make_range(I(halfedge(v,g),g), I(halfedge(v,g),g,1));
}

/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all vertices adjacent to the vertex `target(h,g)`. 
 */
template <typename Graph>
Iterator_range<Vertex_around_target_iterator<Graph> >
vertices_around_target(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Vertex_around_target_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}

template <typename Graph>
Iterator_range<Vertex_around_target_iterator<Graph> >
vertices_around_target(typename boost::graph_traits<Graph>::vertex_descriptor v, const Graph& g)
{
  typedef Vertex_around_target_iterator<Graph> I;
  return make_range(I(halfedge(v,g),g), I(halfedge(v,g),g,1));
}
/**  
 * \ingroup PkgBGLIterators
 * returns an iterator range over all vertices adjacent to the face `face(h,g)`. 
 */
template <typename Graph>
Iterator_range<Vertex_around_face_iterator<Graph> >
vertices_around_face(typename boost::graph_traits<Graph>::halfedge_descriptor h, const Graph& g)
{
  typedef Vertex_around_face_iterator<Graph> I;
  return make_range(I(h,g), I(h,g,1));
}


template <class Graph>
class Out_edge_iterator
  : public boost::iterator_adaptor<
             Out_edge_iterator<Graph>                                // Derived
             , Halfedge_around_target_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::edge_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::edge_descriptor  // Reference
             >
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
private:
  internal::Opposite_edge<Graph> opp;
public:
  Out_edge_iterator()
  {}
  
  Out_edge_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : Out_edge_iterator::iterator_adaptor_(Halfedge_around_target_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), opp(g) {}
  
  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (Out_edge_iterator::*bool_type)() const;
  
  void this_type_does_not_support_comparisons() const {}
  
  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &Out_edge_iterator::this_type_does_not_support_comparisons : 0;
  }
  
  
private:
  friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::edge_descriptor dereference() const { return opp(*this->base_reference()); }
};




template <class Graph>
class In_edge_iterator
  : public boost::iterator_adaptor<
             In_edge_iterator<Graph>                                 // Derived
             , Halfedge_around_target_iterator<Graph>                // Base
             , typename boost::graph_traits<Graph>::edge_descriptor  // Value
             , std::bidirectional_iterator_tag                       // CategoryOrTraversal
             , typename boost::graph_traits<Graph>::edge_descriptor  // Reference
             >
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
private:
  internal::Edge<Graph> fct;
public:
  In_edge_iterator()
  {}
  
  In_edge_iterator(halfedge_descriptor h, const Graph& g, int n = 0)
    : In_edge_iterator::iterator_adaptor_(Halfedge_around_target_iterator<Graph>(h,g,(h==halfedge_descriptor())?1:n)), fct(g)
  {}

  // design patter: "safe bool"
  // will be replaced by explicit operator bool with C++11
  typedef void (In_edge_iterator::*bool_type)() const;

  void this_type_does_not_support_comparisons() const {}

  operator bool_type() const
  {
    return (! (this->base_reference() == NULL)) ?
      &In_edge_iterator::this_type_does_not_support_comparisons : 0;
  }

 private:
    friend class boost::iterator_core_access;
  typename  boost::graph_traits<Graph>::edge_descriptor dereference() const { return fct(*this->base_reference()); }
};






} // CGAL


#endif /* CGAL_BGL_ITERATORS_H */
