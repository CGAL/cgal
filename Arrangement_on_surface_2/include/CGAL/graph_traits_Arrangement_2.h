// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_BOOST_GRAPH_TRAITS_ARRANGEMENT_2_H
#define CGAL_BOOST_GRAPH_TRAITS_ARRANGEMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>


/*! \file
 * Definition of the specialized boost::graph_traits<Arrangement_2> class.
 */

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_accessor.h>

namespace boost {

/*! \class
 * Specialization of the BGL graph-traits template, which serves as a (primal)
 * adapter for Arrangment_on_surface_2, where the valid arrangement vertices
 * correspond to graph verices and arrangement halfedges correspond to
 * arrangement edges.
 * Note that non-fictitious vertices at infinity are also considered as graph
 * vertices, as they have incident non-fictitious edges.
 * As halfedges are directed, we consider the graph as directed. We also allow
 * parallel edges, as two or more different arrangement edges (halfedge pairs)
 * may connect two adjacent vertices.
 */
template <class GeomTraits, class TopTraits>
class graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
{
public:
  
  typedef GeomTraits                                 Geometry_traits_2;
  typedef TopTraits                                  Topology_traits;
  typedef CGAL::Arrangement_on_surface_2<Geometry_traits_2, Topology_traits>
                                                     Arrangement_on_surface_2;

private:

  typedef CGAL::Arr_accessor<Arrangement_on_surface_2>       Arr_accessor;
  typedef typename Arrangement_on_surface_2::Vertex_handle   Vertex_handle;
  typedef typename Arr_accessor::Valid_vertex_iterator       Vertex_iterator;
  typedef typename Arrangement_on_surface_2::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement_on_surface_2::Halfedge_iterator
                                                             Halfedge_iterator;
  typedef typename Arrangement_on_surface_2::Halfedge_around_vertex_circulator
                                             Halfedge_around_vertex_circulator;

  /*! \struct
   * Define the arrangement traversal category, which indicates the arrangement
   * models the BidirectionalGraph concept as well as the VertexListGraph and
   * EdgeListGraph concepts.
   */
  struct Arr_traversal_category : 
    public virtual boost::bidirectional_graph_tag,   // This tag refines the
                                                     // incidence_graph_tag.
    public virtual boost::vertex_list_graph_tag,  // Can iterate over vertices.
    public virtual boost::edge_list_graph_tag     // Can iterate over edges.
  {};

  /*! \class
   * Iteratator over all outgoing halfedges around a given vertex., skipping
   * fictitious halfedges.
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
    int                                _cend;     // The end counter.
    Halfedge_handle                    _hh;       // The current halfedge.

public:

    /*! Default constructor. */
    Halfedge_around_vertex_iterator () :
      _counter(-1),
      _cend(-1)
    {}

    /*!
     * Constructor. 
     * \param circ A ciruclator for the halfedges around a vertex.
     * \param out_edges Do we need the outgoing or the ingoing halfedges.
     * \param counter A counter associated with the iterator.
     * \param cend The past-the-end counter value.
     */
    Halfedge_around_vertex_iterator (Halfedge_around_vertex_circulator circ,
                                     bool out_edges,
                                     int counter,
                                     int cend) :
      _circ (circ),
      _out (out_edges),
      _counter (counter),
      _cend (cend)
    {
      if (_circ->is_fictitious() && _counter < _cend)
        ++(*this);

      if (out_edges)
        _hh = _circ->twin();
      else
        _hh = _circ;
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
      return _hh;
    }

    pointer operator-> () const
    {
      return &_hh;
    }
    
    /* Increment operators. */
    Self& operator++()
    {
      do
      {
        ++_circ;
        ++_counter;
        
        if (_out)
          _hh = _circ->twin();
        else
          _hh = _circ;

      } while (_circ->is_fictitious() && _counter < _cend);

      return (*this);
    }

    Self operator++ (int )
    {
      Self tmp = *this;
      
      do
      {
        ++_circ;
        ++_counter;
        
        if (_out)
          _hh = _circ->twin();
        else
          _hh = _circ;

      } while (_circ->is_fictitious() && _counter < _cend);

      return tmp;
    }
  };

  // Data members:
  Arrangement_on_surface_2    *p_arr;
  Arr_accessor                 arr_access;

public:

  // Types required of the Graph concept:
  typedef typename Arrangement_on_surface_2::Vertex_handle
                                                        vertex_descriptor;
  typedef boost::directed_tag                           directed_category;
  typedef boost::allow_parallel_edge_tag                edge_parallel_category;

  typedef Arr_traversal_category                        traversal_category;

  // Types required by the IncidenceGraph concept:
  typedef typename Arrangement_on_surface_2::Halfedge_handle
                                                        edge_descriptor;
  typedef Halfedge_around_vertex_iterator               out_edge_iterator;
  typedef typename Arrangement_on_surface_2::Size       degree_size_type;

  // Types required by the BidirectionalGraph concept:
  typedef Halfedge_around_vertex_iterator               in_edge_iterator;

  // Types required by the VertexListGraph concept:
  typedef boost::counting_iterator<Vertex_iterator>     vertex_iterator;
  typedef typename Arrangement_on_surface_2::Size       vertices_size_type;

  // Types required by the EdgeListGraph concept:
  typedef boost::counting_iterator<Halfedge_iterator>   edge_iterator;
  typedef typename Arrangement_on_surface_2::Size       edges_size_type;

  // Types not required by any of these concepts:
  typedef void                                          adjacency_iterator;

  /*! Constructor. */
  graph_traits (const Arrangement_on_surface_2& arr) :
    p_arr (const_cast<Arrangement_on_surface_2 *> (&arr)),
    arr_access (const_cast<Arrangement_on_surface_2&> (arr))
  {}

  /*! Nulls */
  static vertex_descriptor null_vertex() { return vertex_descriptor(); }

  /*! Traverse the vertices. */
  vertices_size_type number_of_vertices()
  {
    return arr_access.number_of_valid_vertices();
  }

  vertex_iterator vertices_begin()
  {
    return arr_access.valid_vertices_begin();
  }

  vertex_iterator vertices_end()
  {
    return arr_access.valid_vertices_end();
  }

  /*! Traverse the edges. */
  edge_iterator edges_begin()
  {
    return p_arr->halfedges_begin();
  }

  edge_iterator edges_end()
  {
    return p_arr->halfedges_end();
  }

  /*! Get the vertex degree (in degree or out degree). */
  degree_size_type degree (vertex_descriptor v)
  {
    if (v->is_isolated())
      return 0;

    Halfedge_around_vertex_circulator   first = v->incident_halfedges();
    Halfedge_around_vertex_circulator   circ = first;
    degree_size_type                    deg = 0;
    
    do
    {
      if (! circ->is_fictitious())
        deg++;

      ++circ;
    } while (circ != first);

    return deg;
  }

  /*! Traverse the outgoing halfedges of a given vertex. */
  out_edge_iterator out_edges_begin (vertex_descriptor v)
  {
    if (v->is_isolated())
      return out_edge_iterator();

    return out_edge_iterator (v->incident_halfedges(), true, 0, static_cast<int>(v->degree()));
  }

  out_edge_iterator out_edges_end (vertex_descriptor v)
  {
    if (v->is_isolated())
      return out_edge_iterator ();

    const int deg = static_cast<int>(v->degree());
    return out_edge_iterator (v->incident_halfedges(), true, deg, deg);
  }

  /*! Traverse the ingoing halfedges of a given vertex. */
  in_edge_iterator in_edges_begin (vertex_descriptor v)
  {
    if (v->is_isolated())
      return in_edge_iterator();

    const int deg = static_cast<int>(v->degree());
    return in_edge_iterator (v->incident_halfedges(), false, 0, deg);
  }

  in_edge_iterator in_edges_end (vertex_descriptor v)
  {
    if (v->is_isolated())
      return in_edge_iterator ();

    const int deg = static_cast<int>(v->degree());
    return in_edge_iterator (v->incident_halfedges(), false, deg, deg);
  }
};

/*! \class
 * Specialization of the BGL graph-traits template, which serves as a (primal)
 * adapter for Arrangment_2, where the arrangement vertices correspond to graph
 * verices and arrangement halfedges correspond to arrangement edges.
 */
template <class Traits_, class Dcel_>
class graph_traits<CGAL::Arrangement_2<Traits_, Dcel_> > :
    public graph_traits<CGAL::Arrangement_on_surface_2<
      typename CGAL::Arrangement_2<Traits_, Dcel_>::Geometry_traits_2,
      typename CGAL::Arrangement_2<Traits_, Dcel_>::Topology_traits> >
{
  typedef Traits_                                                     Traits_2;
  typedef Dcel_                                                       Dcel;
  typedef graph_traits<CGAL::Arrangement_on_surface_2<
    typename CGAL::Arrangement_2<Traits_2, Dcel>::Geometry_traits_2,
    typename CGAL::Arrangement_2<Traits_2, Dcel>::Topology_traits> >  Base;

public:

  /*! Constructor. */
  graph_traits (const CGAL::Arrangement_2<Traits_2, Dcel>& arr) :
    Base (arr)
  {}
};

} // namespace boost

namespace CGAL {

// Functions required by the IncidenceGraph concept:
// -------------------------------------------------

/*!
 * Get the out-degree of a vertex in a given arrangement.
 * \param v The vertex.
 * \param arr The arrangement.
 * \param Number of outgoing halfedges from v.
 */
template <class GeomTraits, class TopTraits>
typename
boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >::
                                                              degree_size_type
out_degree (typename
            boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                               TopTraits> >::
                                                           vertex_descriptor v,
            const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return gt_arr.degree (v);
}

/*!
 * Return a range of the out-edges of a vertex given by its descriptor and the
 * arrangement it belongs to.
 * \param v The vertex.
 * \param arr The arrangement.
 * \return A pair of out-edges iterators.
 */
template <class GeomTraits, class TopTraits>
std::pair<typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                         out_edge_iterator,
          typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                         out_edge_iterator>
out_edges (typename
           boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                              TopTraits> >::
                                                          vertex_descriptor v,
           const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return std::make_pair (gt_arr.out_edges_begin (v), gt_arr.out_edges_end (v));
}

/*!
 * Get the source vertex of an arrangement edge.
 * \param e The edge.
 * \param arr The arrangement.
 * \return The source vertex of e.
 */
template <class GeomTraits, class TopTraits>
typename
boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >::
                                                           vertex_descriptor
source (typename
        boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                           TopTraits> >::
                                                           edge_descriptor e,
        const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& /* arr */)
{
  return e->source();
}

/*!
 * Get the target vertex of an arrangement edge.
 * \param e The edge.
 * \param arr The arrangement.
 * \return The source vertex of e.
 */
template <class GeomTraits, class TopTraits>
typename
boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >::
                                                           vertex_descriptor
target (typename
        boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                           TopTraits> >::
                                                           edge_descriptor e,
        const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& /* arr */)
{
  return e->target();
}

// Functions required by the BidirectionalGraph concept:
// -----------------------------------------------------

/*!
 * Get the in-degree of a vertex in a given arrangement.
 * \param v The vertex.
 * \param arr The arrangement.
 * \param Number of ingoing halfedges to v.
 */
template <class GeomTraits, class TopTraits>
typename
boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >::
                                                           degree_size_type
in_degree (typename
           boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                              TopTraits> >::
                                                           vertex_descriptor v,
           const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return gt_arr.degree (v);
}

/*!
 * Return a range of the in-edges of a vertex given by its descriptor and the
 * arrangement it belongs to.
 * \param v The vertex.
 * \param arr The arrangement.
 * \return A pair of in-edges iterators.
 */
template <class GeomTraits, class TopTraits>
std::pair<typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                             in_edge_iterator,
          typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                             in_edge_iterator>
in_edges (typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                           vertex_descriptor v,
          const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return std::make_pair (gt_arr.in_edges_begin (v), gt_arr.in_edges_end (v));
}

/*!
 * Get the degree of a vertex in a given arrangement.
 * \param v The vertex.
 * \param arr The arrangement.
 * \param Number of ingoing and outgoing halfedges incident to v.
 */
template <class GeomTraits, class TopTraits>
typename
boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >::
                                                             degree_size_type
degree (typename
        boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                           TopTraits> >::
                                                           vertex_descriptor v,
        const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return (2 * gt_arr.degree (v));
}

// Functions required by the VertexListGraph concept:
// --------------------------------------------------

/*!
 * Get the number of vertices in the given arrangement. 
 * \param arr The arrangement.
 * \return Number of vertices.
 */
template <class GeomTraits, class TopTraits>
typename
boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >::
                                                            vertices_size_type
num_vertices (const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return gt_arr.number_of_vertices(); 
}

/*!
 * Get the range of vertices of the given arrangement.
 * \param arr The arrangement.
 * \return A pair of vertex iterators.
 */
template <class GeomTraits, class TopTraits>
std::pair<typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                             vertex_iterator,
          typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                             vertex_iterator>
vertices (const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return std::make_pair (gt_arr.vertices_begin(), gt_arr.vertices_end());
}

// Functions required by the EdgeListGraph concept:
// ------------------------------------------------

/*!
 * Get the number of halfedges in the given arrangement. 
 * \param arr The arrangement.
 * \return Number of halfedges (graph edges).
 */
template <class GeomTraits, class TopTraits>
typename
boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >::
                                                               edges_size_type
num_edges (const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  return arr.number_of_halfedges(); 
}

/*!
 * Get the range of halfedges of the given arrangement.
 * \param arr The arrangement.
 * \return A pair of halfedge iterators.
 */
template <class GeomTraits, class TopTraits>
std::pair<typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                             edge_iterator,
          typename
          boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits,
                                                             TopTraits> >::
                                                             edge_iterator>
edges (const CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits>& arr)
{
  boost::graph_traits<CGAL::Arrangement_on_surface_2<GeomTraits, TopTraits> >
    gt_arr (arr);

  return std::make_pair (gt_arr.edges_begin(), gt_arr.edges_end());
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
