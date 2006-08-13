// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_BOOST_GRAPH_TRAITS_ARRANGEMENT_2_H
#define CGAL_BOOST_GRAPH_TRAITS_ARRANGEMENT_2_H

/*! \file
 * Definition of the specialized boost::graph_traits<Arrangement_2> class.
 */

#include <boost/graph/graph_concepts.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/Arrangement_2.h>

namespace boost {

/*! \class
 * Specialization of the BGL graph-traits template, which serve as a (primal)
 * adapter for Arrangment_2, where the arrangement vertices correspond to graph
 * verices and arrangement halfedges correspond to arrangement edges.
 * As halfedges are directed, we consider the graph as directed. We also allow
 * parallel edges, as two or more different arrangement edges (halfedge pairs)
 * may connect two adjacent vertices.
 */
template <class Traits_, class Dcel_>
class graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >
{
public:
  
  typedef Traits_                                    Traits_2;
  typedef Dcel_                                      Dcel;
  typedef CGAL::Arrangement_2<Traits_2, Dcel>        Arrangement_2;

private:

  typedef typename Arrangement_2::Vertex_handle      Vertex_handle;
  typedef typename Arrangement_2::Vertex_iterator    Vertex_iterator;
  typedef typename Arrangement_2::Halfedge_handle    Halfedge_handle;
  typedef typename Arrangement_2::Halfedge_iterator  Halfedge_iterator;
  typedef typename Arrangement_2::Halfedge_around_vertex_circulator
                                             Halfedge_around_vertex_circulator;

  /*! \struct
   * Define the arrangement traversal category, which indicates the arrangement
   * models the BidirectionalGraph concept and the VertexListGraph and
   * EdgeListGraph concepts.
   */
  struct Arr_traversal_category : 
    public virtual boost::bidirectional_graph_tag,   // This tag refines the
                                                     // incidence_graph_tag.
    public virtual boost::vertex_list_graph_tag,  // Can iterate over vertices.
    public virtual boost::edge_list_graph_tag     // Can iterate over edges.
  {};

  /*! \class
   * Iteratator over all outgoing halfedges around a given vertex.
   * This is by adapting the Halfegde_around_vertex_circulator type to an
   * iterator. Moreover, as the circulator goes over all ingoing halfedges
   * of the vertex, the iterator adapter may return their twin halfedges, if
   * we need the outgoing halfedges.  
   */
  class Halfedge_around_vertex_iterator
  {
  public:

    // Type definitions:
    typedef Halfedge_around_vertex_iterator       Self;

    typedef std::forward_iterator_tag             iterator_category;
    typedef Halfedge_handle                       value_type;
    typedef value_type                            reference;
    typedef value_type*                           pointer;
    typedef int                                   difference_type;

  protected:

    Halfedge_around_vertex_circulator  _circ;     // The circulator.
    bool                               _out;      // Do we need the out edges.
    int                                _counter;  // Counter for the edges.
    Halfedge_handle                    _hh;       // The current halfedge.

public:

    /*! Default constructor. */
    Halfedge_around_vertex_iterator () :
      _counter(-1)
    {}

    /*!
     * Constructor. 
     * \param circ A ciruclator for the halfedges around a vertex.
     * \param out_edges Do we need the outgoing or the ingoing halfedges.
     * \param counter A counter associated with the iterator.
     */
    Halfedge_around_vertex_iterator (Halfedge_around_vertex_circulator circ,
                                     bool out_edges,
                                     int counter) :
      _circ (circ),
      _out (out_edges),
      _counter (counter)
    {
      if (out_edges)
        _hh = circ->twin();
      else
        _hh = circ;
    }

    /*! Equality operators. */
    bool operator== (const Self& it) const
    {
      return (_circ == it._circ && _out == it._out && _counter == it._counter);
    }
    
    bool operator!= (const Self& it) const
    {
      return (_circ != it._circ || _out != it._out || _counter != it._counter);
    }
    
    /*! Dereference operators. */
    reference operator* () const
    {
      return (_hh);
    }

    pointer operator-> () const
    {
      return (&_hh);
    }
    
    /* Increment operators. */
    Self& operator++()
    {
      ++_circ;
      ++_counter;

      if (_out)
        _hh = _circ->twin();
      else
        _hh = _circ;

      return (*this);
    }

    Self operator++ (int )
    {
      Self tmp = *this;
      ++_circ;
      ++_counter;

      if (_out)
        _hh = _circ->twin();
      else
        _hh = _circ;

      return (tmp);
    }
  };

  // Data members:
  Arrangement_2    *p_arr;
  
public:

  // Types required of the Graph concept:
  typedef typename Arrangement_2::Vertex_handle         vertex_descriptor;
  typedef boost::directed_tag                           directed_category;
  typedef boost::allow_parallel_edge_tag                edge_parallel_category;

  typedef Arr_traversal_category                        traversal_category;

  // Types required by the IncidenceGraph concept:
  typedef typename Arrangement_2::Halfedge_handle       edge_descriptor;
  typedef Halfedge_around_vertex_iterator               out_edge_iterator;
  typedef typename Arrangement_2::Size                  degree_size_type;

  // Types required by the BidirectionalGraph concept:
  typedef Halfedge_around_vertex_iterator               in_edge_iterator;

  // Types required by the VertexListGraph concept:
  typedef boost::counting_iterator<Vertex_iterator>     vertex_iterator;
  typedef typename Arrangement_2::Size                  vertices_size_type;

  // Types required by the EdgeListGraph concept:
  typedef boost::counting_iterator<Halfedge_iterator>   edge_iterator;
  typedef typename Arrangement_2::Size                  edges_size_type;

  // Types not required by any of these concepts:
  typedef void                                          adjacency_iterator;

  /*! Constructor. */
  graph_traits (const Arrangement_2& arr) :
    p_arr (const_cast<Arrangement_2 *> (&arr))
  {}

  /*! Traverse the vertices. */
  vertex_iterator vertices_begin()
  {
    return (p_arr->vertices_begin());
  }

  vertex_iterator vertices_end()
  {
    return (p_arr->vertices_end());
  }

  /*! Traverse the edges. */
  edge_iterator edges_begin()
  {
    return (p_arr->halfedges_begin());
  }

  edge_iterator edges_end()
  {
    return (p_arr->halfedges_end());
  }

  /*! Traverse the outgoing halfedges of a given vertex. */
  out_edge_iterator out_edges_begin (vertex_descriptor v)
  {
    if (v->is_isolated())
      return out_edge_iterator();

    return (out_edge_iterator (v->incident_halfedges(), true, 0));
  }

  out_edge_iterator out_edges_end (vertex_descriptor v)
  {
    if (v->is_isolated())
      return out_edge_iterator ();

    return (out_edge_iterator (v->incident_halfedges(), true, v->degree()));
  }

  /*! Traverse the ingoing halfedges of a given vertex. */
  in_edge_iterator in_edges_begin (vertex_descriptor v)
  {
    if (v->is_isolated())
      return in_edge_iterator();

    return (in_edge_iterator (v->incident_halfedges(), false, 0));
  }

  in_edge_iterator in_edges_end (vertex_descriptor v)
  {
    if (v->is_isolated())
      return in_edge_iterator ();

    return (in_edge_iterator (v->incident_halfedges(), v->degree()));
  }

};

// Functions required by the IncidenceGraph concept:
// -------------------------------------------------

/*!
 * Get the out-degree of a vertex in a given arrangement.
 * \param v The vertex.
 * \param arr The arrangement.
 * \param Number of outgoing halfedges from v.
 */
template <class Traits_, class Dcel_>
typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::degree_size_type
out_degree (typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                           vertex_descriptor v,
            const CGAL::Arrangement_2<Traits_, Dcel_>& /* arr */)
{
  return (v->degree());
}

/*!
 * Return a range of the out-edges of a vertex given by its descriptor and the
 * arrangement it belongs to.
 * \param v The vertex.
 * \param arr The arrangement.
 * \return A pair of out-edges iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                             out_edge_iterator,
          typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                             out_edge_iterator>
out_edges (typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                           vertex_descriptor v,
           const CGAL::Arrangement_2<Traits_, Dcel_>& arr)
{
  graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >  gt_arr (arr);

  return (std::make_pair (gt_arr.out_edges_begin (v),
                          gt_arr.out_edges_end (v)));
}

/*!
 * Get the source vertex of an arrangement edge.
 * \param e The edge.
 * \param arr The arrangement.
 * \return The source vertex of e.
 */
template <class Traits_, class Dcel_>
typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::vertex_descriptor
source (typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                           edge_descriptor e,
        const CGAL::Arrangement_2<Traits_, Dcel_>& /* arr */)
{
  return (e->source());
}

/*!
 * Get the target vertex of an arrangement edge.
 * \param e The edge.
 * \param arr The arrangement.
 * \return The source vertex of e.
 */
template <class Traits_, class Dcel_>
typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::vertex_descriptor
target (typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                           edge_descriptor e,
        const CGAL::Arrangement_2<Traits_, Dcel_>& /* arr */)
{
  return (e->target());
}

// Functions required by the BidirectionalGraph concept:
// -----------------------------------------------------

/*!
 * Get the in-degree of a vertex in a given arrangement.
 * \param v The vertex.
 * \param arr The arrangement.
 * \param Number of ingoing halfedges from v.
 */
template <class Traits_, class Dcel_>
typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::degree_size_type
in_degree (typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                           vertex_descriptor v,
           const CGAL::Arrangement_2<Traits_, Dcel_>& /* arr */)
{
  return (v->degree());
}

/*!
 * Return a range of the in-edges of a vertex given by its descriptor and the
 * arrangement it belongs to.
 * \param v The vertex.
 * \param arr The arrangement.
 * \return A pair of in-edges iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                             in_edge_iterator,
          typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                             in_edge_iterator>
in_edges (typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                           vertex_descriptor v,
          const CGAL::Arrangement_2<Traits_, Dcel_>& arr)
{
  graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >  gt_arr (arr);

  return (std::make_pair (gt_arr.in_edges_begin (v),
                          gt_arr.in_edges_end (v)));
}

/*!
 * Get the degree of a vertex in a given arrangement.
 * \param v The vertex.
 * \param arr The arrangement.
 * \param Number of ingoing and outgoing halfedges incident to v.
 */
template <class Traits_, class Dcel_>
typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::degree_size_type
degree (typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                           vertex_descriptor v,
        const CGAL::Arrangement_2<Traits_, Dcel_>& /* arr */)
{
  return (2 * v->degree());
}

// Functions required by the VertexListGraph concept:
// --------------------------------------------------

/*!
 * Get the number of vertices in the given arrangement. 
 * \param arr The arrangement.
 * \return Number of vertices.
 */
template <class Traits_, class Dcel_>
typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::vertices_size_type
num_vertices (const CGAL::Arrangement_2<Traits_, Dcel_>& arr)
{
  return (arr.number_of_vertices()); 
}

/*!
 * Get the range of vertices of the given arrangement.
 * \param arr The arrangement.
 * \return A pair of vertex iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                               vertex_iterator,
          typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                               vertex_iterator>
vertices (const CGAL::Arrangement_2<Traits_, Dcel_>& arr)
{
  graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >  gt_arr (arr);

  return (std::make_pair (gt_arr.vertices_begin(),
                          gt_arr.vertices_end()));
}

// Functions required by the EdgeListGraph concept:
// ------------------------------------------------

/*!
 * Get the number of halfedges in the given arrangement. 
 * \param arr The arrangement.
 * \return Number of halfedges (graph edges).
 */
template <class Traits_, class Dcel_>
typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::edges_size_type
num_edges (const CGAL::Arrangement_2<Traits_, Dcel_>& arr)
{
  return (arr.number_of_halfedges()); 
}

/*!
 * Get the range of halfedges of the given arrangement.
 * \param arr The arrangement.
 * \return A pair of halfedge iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                               edge_iterator,
          typename graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >::
                                                               edge_iterator>
edges (const CGAL::Arrangement_2<Traits_, Dcel_>& arr)
{
  graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> >  gt_arr (arr);

  return (std::make_pair (gt_arr.edges_begin(),
                          gt_arr.edges_end()));
}

}; // namespace boost

#endif
