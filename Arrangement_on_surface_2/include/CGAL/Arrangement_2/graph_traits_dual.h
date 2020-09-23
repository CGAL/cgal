// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Ron Wein         <wein@post.tau.ac.il>
//             Ophir Setter     <ophirset@post.tau.ac.il>
//             Sebastien Loriot <sebastien.loriot@cgal.org>
//             Efi Fogel        <efifogel@gmail.com>

// This file contains the follwoing three parts:
// 1. The common base class template of the specialized
//      Dual<specialized-arrangement> class template.
//
// 2. The common base class template of the specialized
//      boost::graph_traits<Dual<specialized-arrangement> > class template.
//
// 3. Macro definitions of free Function templates required by
//    the various Boost Graph concepts. There is one macro per required function
//    template. Each macro accepts the name of a template class, an instance of
//    which represents an arrangement data structure, e.g., Arrangement_2. The
//    definitios of the free functions templates for a given arrangement data
//    strcture must be present when a dual of this data structure is defined.

#include <CGAL/license/Arrangement_on_surface_2.h>

#ifndef CGAL_GRAPH_TRAITS_DUAL_H
#define CGAL_GRAPH_TRAITS_DUAL_H

namespace CGAL {

// Forward declaration.
template <class Type> class Dual;

/*! \class
 * Generic implementation of the common base class template of the specialized
 * Dual<specialized-arrangement> class template.
 */
template <typename Arrangement_>
class Dual_arrangement_on_surface {
public:
  typedef Arrangement_                                  Arrangement;
  typedef typename Arrangement::Geometry_traits_2       Geometry_traits_2;
  typedef typename Arrangement::Topology_traits         Topology_traits;

  typedef typename Arrangement::Size                    Size;
  typedef typename Arrangement::Face_handle             Vertex_handle;
  typedef typename Arrangement::Halfedge_handle         Edge_handle;

  typedef typename Arrangement::Face_iterator           Vertex_iterator;
  typedef typename Arrangement::Halfedge_iterator       Edge_iterator;

protected:
  typedef typename Arrangement::Face_handle             Face_handle;
  typedef typename Arrangement::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Arrangement::Outer_ccb_iterator      Outer_ccb_iterator;
  typedef typename Arrangement::Inner_ccb_iterator      Inner_ccb_iterator;

  /*! \class
   * Iterator over the neighbors of a dual vertex (a face in the primal
   * arrangement).
   * These neighbors are the adjacent faces along the outer boundaries of the
   * face and its inner boundaries.
   */
  class Face_neighbor_iterator {
    typedef Face_neighbor_iterator               Self;

  public:
    typedef std::forward_iterator_tag            iterator_category;
    typedef Edge_handle                          value_type;
    typedef value_type                           reference;
    typedef value_type*                          pointer;
    typedef int                                  difference_type;

  private:
    Outer_ccb_iterator _outer_ccb_iter;
    Inner_ccb_iterator _inner_ccb_iter;
    Ccb_halfedge_circulator _ccb_curr;
    Ccb_halfedge_circulator _ccb_first;
    Face_handle _face;
    bool _out;
    Edge_handle _hh;
    bool _end;

  public:
    /*! Default constructor. */
    Face_neighbor_iterator() : _end (true) {}

    /*! Constructor.
     * \param face The face (dual vertex).
     * \param out_edges Do we need the outgoing or the ingoing halfedges.
     * \param start Should we start traversing the edges.
     *              If false, we construct a past-the-end iterator.
     */
    Face_neighbor_iterator(Face_handle face, bool out_edges, bool start) :
      _face(face),
      _out(out_edges),
      _end(! start)
    {
      CGAL_precondition(! face->is_fictitious());

      if (start) {
        _outer_ccb_iter = _face->outer_ccbs_begin();
        _inner_ccb_iter = _face->inner_ccbs_begin();

        if (_outer_ccb_iter != _face->outer_ccbs_end()) {
          // Start from the first outer CCB, if one exists.
          _ccb_curr = _ccb_first = *_outer_ccb_iter;
        }
        else if (_inner_ccb_iter != face->inner_ccbs_end()) {
          // Otherwise, start from the first inner CCB.
          _ccb_curr = _ccb_first = *_inner_ccb_iter;
        }
        else {
          // In this case there are no CCBs to traverse:
          _end = true;
          return;
        }

        _hh = this->_dereference();

        // In case the incident face of the twin halfedge is fictitious,
        // skip it and proceed to the next edge.
        if (_hh->is_fictitious()) ++(*this);
      }
      else { // end iterator.
        _outer_ccb_iter = _face->outer_ccbs_end();
        _inner_ccb_iter = _face->inner_ccbs_end();
      }
    }

    /*! Equality operators. */
    bool operator==(const Self& it) const { return (this->_equal(it)); }

    bool operator!= (const Self& it) const { return (! this->_equal(it)); }

    /*! Dereference operators. */
    reference operator*() const { return (_hh); }

    pointer operator->() const { return (&_hh); }

    /* Increment operators. */
    Self& operator++()
    {
      do {
        this->_increment();
        if (_end) return (*this);
        _hh = this->_dereference();
      } while (_hh->is_fictitious());
      return (*this);
    }

    Self operator++(int )
    {
      Self tmp = *this;
      do {
        this->_increment();
        if (_end) return (tmp);
        _hh = this->_dereference();
      } while (_hh->is_fictitious());
      return (tmp);
    }

  private:
    /*! Check two iterators for equality. */
    bool _equal(const Self& it) const
    {
      return (_out == it._out && _face == it._face && ((_end && it._end) ||
               (_outer_ccb_iter == it._outer_ccb_iter &&
                _inner_ccb_iter == it._inner_ccb_iter &&
                _ccb_curr == it._ccb_curr)));
    }

    /*! Derefernce the current circulator. */
    Edge_handle _dereference() const
    {
      if (_out) return (_ccb_curr);
      else return (_ccb_curr->twin());
    }

    // Increments of the iterator.
    void _increment()
    {
      CGAL_assertion(! _end);

      // If we have not traversed the entire CCB in full, move to the next
      // halfedge along the current CCB.
      ++_ccb_curr;

      if (_ccb_curr != _ccb_first) return;

      // In this case we have completed the current CCB and we have to move
      // to the next one.
      if (_outer_ccb_iter != _face->outer_ccbs_end()) {
        // Try to move to the next outer CCB.
        ++_outer_ccb_iter;
        if (_outer_ccb_iter != _face->outer_ccbs_end()) {
          _ccb_curr = _ccb_first = *_outer_ccb_iter;
          return;
        }

        // In this case we start traversing the inner CCBs.
        if (_inner_ccb_iter != _face->inner_ccbs_end()) {
          CGAL_assertion (_inner_ccb_iter == _face->inner_ccbs_begin());

          // Otherwise, start from the first inner CCB.
          _ccb_curr = _ccb_first = *_inner_ccb_iter;
          return;
        }
      }
      else if (_inner_ccb_iter != _face->inner_ccbs_end()) {

        // In this case we have already traversed all outer CCBs (and at least
        // one inner CCB), so we try to move to the next inner CCB.
        ++_inner_ccb_iter;
        if (_inner_ccb_iter != _face->inner_ccbs_end()) {
          // Otherwise, start from the first inner CCB.
          _ccb_curr = _ccb_first = *_inner_ccb_iter;
          return;
        }
      }

      // In this case we finished traversing all outer and inner CCBs:
      _end = true;
      return;
    }
  };

  // Data members:
  mutable Arrangement* p_arr;           // The primal arrangement.

public:
  typedef Face_neighbor_iterator            Incident_edge_iterator;

  /*! Default constructor. */
  Dual_arrangement_on_surface() : p_arr(nullptr) {}

  /*! Constructor from an arrangement. */
  Dual_arrangement_on_surface(const Arrangement& arr) :
    p_arr(const_cast<Arrangement*>(&arr))
  {}

  /*! Obtain the primal arrangement (const version). */
  const Arrangement* arrangement() const { return (p_arr); }

  /*! Obtain the primal arrangement (non-const version). */
  Arrangement* arrangement() { return (p_arr); }

  /*! Obtain the number of vertices (face of the primal arrangement). */
  Size number_of_vertices() const { return (p_arr->number_of_faces()); }

  /*! Obtain the begin iterator of the vertices of the dual arrangement
   * (faces of the primal arrangement).
   */
  Vertex_iterator vertices_begin() const { return (p_arr->faces_begin()); }

  /*! Obtain the pass-the-end iterator of the vertices of the dual arrangement
   * (faces of the primal arrangement).
   */
  Vertex_iterator vertices_end() const { return (p_arr->faces_end()); }

  /*! Obtain the number of edges. */
  Size number_of_edges () const { return (p_arr->number_of_halfedges()); }

  /*! Obtain the begin iterator of the edges of the dual arrangement. */
  Edge_iterator edges_begin() const { return (p_arr->halfedges_begin()); }

  /*! Obtain the pass-the-end iterator of the edges of the dual arrangement. */
  Edge_iterator edges_end() const { return (p_arr->halfedges_end()); }

  /*! Obtain the dual vertex-degree (number of edges forming the face boundary).
   */
  Size degree(Vertex_handle v) const
  {
    Incident_edge_iterator begin = Incident_edge_iterator(v, true, true);
    Incident_edge_iterator end = Incident_edge_iterator(v, false, true);
    Size deg = 0;

    while (begin != end) {
      deg++;
      ++begin;
    }

    return (deg);
  }

  /*! Traverse the outgoing edges of a given vertex. */
  Incident_edge_iterator out_edges_begin(Vertex_handle v) const
  { return (Incident_edge_iterator (v, true, true)); }

  Incident_edge_iterator out_edges_end(Vertex_handle v) const
  { return (Incident_edge_iterator (v, true, false)); }

  /*! Traverse the ingoing edges of a given vertex. */
  Incident_edge_iterator in_edges_begin(Vertex_handle v) const
  { return (Incident_edge_iterator (v, false, true)); }

  Incident_edge_iterator in_edges_end(Vertex_handle v) const
  { return (Incident_edge_iterator (v, false, false)); }
};

} //namespace CGAL

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/iterator/counting_iterator.hpp>

namespace CGAL {

/*! \class
 * The common base class template of the specialized
 *   boost::graph_traits<Dual<specialized-arrangement> > class template.
 * The latter serves as a dual adapter for the specialied arrangment, where the
 * valid arrangement faces correspond to graph verices, and two graph vertices
 * are connected if the two corrsponding faces are adjacent.
 * We consider the graph as directed. We also allow parallel edges, as two
 * faces may have more than one common edges.
 */
template <typename Arrangement_>
class Graph_traits_dual_arr_on_surface_impl {
public:
  typedef Arrangement_                                  Arrangement;
  typedef typename Arrangement::Geometry_traits_2       Geometry_traits_2;
  typedef typename Arrangement::Topology_traits         Topology_traits;
  typedef Dual_arrangement_on_surface<Arrangement>      Dual_arr_2;

private:
  typedef typename Dual_arr_2::Vertex_iterator          Vertex_iterator;
  typedef typename Dual_arr_2::Edge_iterator            Edge_iterator;
  typedef typename Dual_arr_2::Incident_edge_iterator   Incident_edge_iterator;

  /*! \struct
   * Define the arrangement traversal category, which indicates the arrangement
   * models the BidirectionalGraph concept and the VertexListGraph and
   * EdgeListGraph concepts.
   */
  struct Dual_arr_traversal_category :
    public virtual boost::bidirectional_graph_tag, // This tag refines the
                                                   // incidence_graph_tag.
    public virtual boost::vertex_list_graph_tag,   // Can iterate over vertices.
    public virtual boost::edge_list_graph_tag      // Can iterate over edges.
  {};

public:
  // Types required of the Graph concept:
  typedef typename Dual_arr_2::Vertex_handle            vertex_descriptor;
  typedef boost::directed_tag                           directed_category;
  typedef boost::allow_parallel_edge_tag                edge_parallel_category;

  typedef Dual_arr_traversal_category                   traversal_category;

  // Types required by the IncidenceGraph concept:
  typedef typename Dual_arr_2::Edge_handle              edge_descriptor;
  typedef Incident_edge_iterator                        out_edge_iterator;
  typedef typename Dual_arr_2::Size                     degree_size_type;

  // Types required by the BidirectionalGraph concept:
  typedef Incident_edge_iterator                        in_edge_iterator;

  // Types required by the VertexListGraph concept:
  typedef boost::counting_iterator<Vertex_iterator>     vertex_iterator;
  typedef typename Dual_arr_2::Size                     vertices_size_type;

  // Types required by the EdgeListGraph concept:
  typedef boost::counting_iterator<Edge_iterator>       edge_iterator;
  typedef typename Dual_arr_2::Size                     edges_size_type;

  // Types not required by any of these concepts:
  typedef void                                          adjacency_iterator;
};

}

// Macro definitions of free Function templates required by the various Boost
// Graph concepts. Each macro provides the required free function template with
// the specific interface and its implementation.
// We use macros (and not base functions similar to
// Graph_traits_dual_arr_on_surface_impl) for simplicity. The implementation
// is typically a one-liner. However, the interface is typically several lines
// of code. The alternative implementation (with base common functions) while
// being more safe (provided tight type checking) would consume many repeated
// lines of code.

// Functions required by the IncidenceGraph concept:
// -------------------------------------------------

/*! Obtain the out-degree of a vertex in a given dual arrangement.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \param Number of halfedges around the boundary of the primal face.
 */
#define CGAL_DUAL_ARRANGEMENT_2_OUT_DEGREE(T)                                 \
template <typename T1, typename T2>                                           \
typename boost::graph_traits<Dual<T<T1, T2> > >::degree_size_type             \
out_degree(typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_descriptor v,\
           const Dual<T<T1, T2> >& darr)                                      \
{ return darr.degree(v); }

/*! Return a range of the out-edges of a vertex given by its descriptor and the
 * dual arrangement it belongs to.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \return A pair of out-edges iterators.
 */
#define CGAL_DUAL_ARRANGEMENT_2_OUT_EDGES(T)                                  \
template <typename T1, typename T2>                                           \
std::pair<typename boost::graph_traits<Dual<T<T1, T2> > >::out_edge_iterator, \
          typename boost::graph_traits<Dual<T<T1, T2> > >::out_edge_iterator> \
out_edges(typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_descriptor v,\
          const Dual<T<T1, T2> >& darr)                                       \
{ return std::make_pair(darr.out_edges_begin(v), darr.out_edges_end(v)); }    \

/*! Obtain the source vertex of a dual arrangement edge.
 * \param e The edge.
 * \param darr The dual arrangement.
 * \return The incident face of e in the primal arrangement.
 */
#define CGAL_DUAL_ARRANGEMENT_2_SOURCE(T)                                     \
template <typename T1, typename T2>                                           \
typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_descriptor            \
source(typename boost::graph_traits<Dual<T<T1, T2> > >::edge_descriptor e,    \
       const Dual<T<T1, T2> >& /* darr */)                                    \
{ return e->face(); }

/*! Obtain the target vertex of a dual arrangement edge.
 * \param e The edge.
 * \param darr The dual arrangement.
 * \return The incident face of the twin of e in the primal arrangement.
 */
#define CGAL_DUAL_ARRANGEMENT_2_TARGET(T)                                     \
template <typename T1, typename T2>                                           \
typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_descriptor            \
target(typename boost::graph_traits<Dual<T<T1, T2> > >::edge_descriptor e,    \
       const Dual<T<T1, T2> >& /* darr */)                                    \
{ return e->twin()->face(); }

// Functions required by the BidirectionalGraph concept:
// -----------------------------------------------------

/*! Obtain the in-degree of a vertex in a given dual arrangement.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \param Number of halfedges around the boundary of the primal face.
 */
#define CGAL_DUAL_ARRANGEMENT_2_IN_DEGREE(T)                                  \
template <typename T1, typename T2>                                           \
typename boost::graph_traits<Dual<T<T1, T2> > >::degree_size_type             \
in_degree(typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_descriptor v,\
          const Dual<T<T1, T2> >& darr)                                       \
{ return darr.degree(v); }

/*! Return a range of the in-edges of a vertex given by its descriptor and the
 * dual arrangement it belongs to.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \return A pair of in-edges iterators.
 */
#define CGAL_DUAL_ARRANGEMENT_2_IN_EDGES(T)                                   \
template <typename T1, typename T2>                                           \
std::pair<typename boost::graph_traits<Dual<T<T1, T2> > >::in_edge_iterator,  \
          typename boost::graph_traits<Dual<T<T1, T2> > >::in_edge_iterator>  \
in_edges(typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_descriptor v,\
         const Dual<T<T1, T2> >& darr)                                        \
{ return std::make_pair(darr.in_edges_begin(v), darr.in_edges_end(v)); }

/*! Obtain the degree of a vertex in a given dual arrangement.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \param Number of ingoing and outgoing halfedges incident to v.
 */
#define CGAL_DUAL_ARRANGEMENT_2_DEGREE(T)                                     \
template <typename T1, typename T2>                                           \
typename boost::graph_traits<Dual<T<T1, T2> > >::degree_size_type             \
degree(typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_descriptor v,  \
       const Dual<T<T1, T2> >& darr)                                          \
{ return (2 * darr.degree(v)); }

// Functions required by the VertexListGraph concept:
// --------------------------------------------------

/*! Obtain the number of vertices in the given dual arrangement.
 * \param darr The dual arrangement.
 * \return Number of faces in the primal arrangement.
 */
#define CGAL_DUAL_ARRANGEMENT_2_NUM_VERTICES(T)                               \
template <typename T1, typename T2>                                           \
typename boost::graph_traits<Dual<T<T1, T2> > >::vertices_size_type           \
num_vertices(const Dual<T<T1, T2> >& darr)                                    \
{ return darr.number_of_vertices(); }

/*! Obtain the range of vertices of the given dual arrangement.
 * \param darr The dual arrangement.
 * \return A pair of vertex iterators.
 */
#define CGAL_DUAL_ARRANGEMENT_2_VERTICES(T)                                   \
template <typename T1, typename T2>                                           \
std::pair<typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_iterator,   \
          typename boost::graph_traits<Dual<T<T1, T2> > >::vertex_iterator>   \
vertices(const Dual<T<T1, T2> >& darr)                                        \
{ return std::make_pair(darr.vertices_begin(), darr.vertices_end()); }

// Functions required by the EdgeListGraph concept:
// ------------------------------------------------

/*! Obtain the number of edges in the given dual arrangement.
 * \param darr The dual arrangement.
 * \return Number of halfedges in the primal arrangement.
 */
#define CGAL_DUAL_ARRANGEMENT_2_NUM_EDGES(T)                                  \
template <typename T1, typename T2>                                           \
typename boost::graph_traits<Dual<T<T1, T2> > >::edges_size_type              \
num_edges(const Dual<T<T1, T2> >& darr)                                       \
{ return darr.number_of_edges(); }

/*! Obtain the range of edges of the given dual arrangement.
 * \param darr The dual arrangement.
 * \return A pair of edge iterators.
 */
#define CGAL_DUAL_ARRANGEMENT_2_EDGES(T)                                      \
template <typename T1, typename T2>                                           \
std::pair<typename boost::graph_traits<Dual<T<T1, T2> > >::edge_iterator,     \
          typename boost::graph_traits<Dual<T<T1, T2> > >::edge_iterator>     \
edges(const Dual<T<T1, T2> >& darr)                                           \
{ return std::make_pair(darr.edges_begin(), darr.edges_end()); }

#endif
