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
// Author(s)     : Ron Wein     <wein@post.tau.ac.il>
//                 Ophir Setter <ophirset@post.tau.ac.il>
//                 Efi Fogel    <efif@post.tau.ac.il>

#ifndef CGAL_BOOST_GRAPH_TRAITS_DUAL_ARRANGEMENT_2_H
#define CGAL_BOOST_GRAPH_TRAITS_DUAL_ARRANGEMENT_2_H

/*! \file
 * Definition of the specialized Dual<Arrangement_2> class,
 * and the specialized graph_traits<Dual<Arrangement_2> >class.
 */

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/Arrangement_2.h>

CGAL_BEGIN_NAMESPACE

// Forward declaration.
template <class Type> class Dual;

/*! \class
 * Specilaization of the Dual<> template for Arrangement_2.
 */
template <class Traits_, class Dcel_>
class Dual<Arrangement_2<Traits_, Dcel_> >
{
public:
  
  typedef Traits_                              Traits_2;
  typedef Dcel_                                Dcel;
  typedef CGAL::Arrangement_2<Traits_2, Dcel>  Arrangement_2;

  typedef typename Arrangement_2::Size                   Size;
  typedef typename Arrangement_2::Face_handle            Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle        Edge_handle;

  typedef typename Arrangement_2::Face_iterator          Vertex_iterator;
  typedef typename Arrangement_2::Halfedge_iterator      Edge_iterator;

protected:

  typedef typename Arrangement_2::Face_handle            Face_handle;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
                                                      Ccb_halfedge_circulator;
  typedef typename Arrangement_2::Inner_ccb_iterator         Hole_iterator;

  /*! \class
   * Iterator over the neighbors of a dual vertex (a face in the primal
   * arrangement).
   * These neighbors are the adjacent faces along the outer boundary of the
   * face (if it is bounded), and its holes.
   */
  class Face_neighbor_iterator
  {
    typedef Face_neighbor_iterator               Self;

  public:

    typedef std::forward_iterator_tag            iterator_category;
    typedef Edge_handle                          value_type;
    typedef value_type                           reference;
    typedef value_type*                          pointer;
    typedef int                                  difference_type;

  private:

    bool                     _ccb_incremented;
    Ccb_halfedge_circulator  _outer_ccb_circ;
    Hole_iterator            _hole_iter;
    Ccb_halfedge_circulator  _curr_hole_circ;
    bool                     _curr_hole_incremented;
    Face_handle              _face;
    bool                     _out;
    Edge_handle              _hh;
    bool                     _end;

  public:

    /*! Default constructor. */
    Face_neighbor_iterator() :
      _end(true)
    {}

    /*!
     * Constructor.
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

      _outer_ccb_circ = face->outer_ccb();
      _ccb_incremented = !start;

      if (start)
      {
        _hole_iter = face->holes_begin();
        _curr_hole_incremented = true;
        
        if (_hole_iter != face->holes_end())
        {
          _curr_hole_circ = *_hole_iter;
          _curr_hole_incremented = false;
        }

        _hh = this->_dereference();

        // In case the incident face of the twin halfedge is fictitious,
        // skip it and proceed to the next edge.
        if (_hh->is_fictitious())
          ++(*this);
      }
      else // end iterator.
      {
        _hole_iter = face->holes_end();
      }
    }  

    /*! Equality operators. */
    bool operator==(const Self& it) const
    {
      return this->_equal(it);
    }
    
    bool operator!=(const Self& it) const
    {
      return (! this->_equal(it));
    }
    
    /*! Dereference operators. */
    reference operator*() const
    {
      return _hh;
    }

    pointer operator->() const
    {
      return &_hh;
    }
    
    /* Increment operators. */
    Self& operator++()
    {
      do
      {
        this->_increment();
        if (_end)
          return (*this);

        _hh = this->_dereference();

      } while (_hh->is_fictitious());

      return (*this);
    }

    Self operator++(int)
    {
      Self tmp = *this;

      do
      {
        this->_increment();
        if (_end)
          return tmp;

        _hh = this->_dereference();

      } while (_hh->is_fictitious());

      return tmp;
    }

  private:

    /*! Check two iterators for equality. */
    bool _equal(const Self& it) const
    {
      return (_out == it._out && _face == it._face &&
              ((_end && it._end) ||
               (_ccb_incremented == it._ccb_incremented &&
                _outer_ccb_circ == it._outer_ccb_circ && 
                _hole_iter == it._hole_iter)));
    }

    /*! Derefernce the current circulator. */
    Edge_handle _dereference() const
    {
      if (! _ccb_incremented ||
          _outer_ccb_circ != _face->outer_ccb())
      {
        if (_out)
          return _outer_ccb_circ;
        else
          return _outer_ccb_circ->twin();
      }
    
      if (_out)
        return _curr_hole_circ;
      else
        return _curr_hole_circ->twin();
    }

    // Increments of the iterator.
    void _increment()
    {
      CGAL_assertion(! _end);

      // If we have not traversed the entire outer CCB (namely this is the
      // first increment operation, or we still have not completed a full
      // cycle around the outer CCB), move to the next halfedge along the
      // outer CCB.
      if (! _ccb_incremented)
      {
        ++_outer_ccb_circ;
        _ccb_incremented = true;

        return;
      }

      if (_outer_ccb_circ != _face->outer_ccb())
      {
        ++_outer_ccb_circ;

        if ((_outer_ccb_circ == _face->outer_ccb()) &&
            (_hole_iter == _face->holes_end()))
          _end = true;
        return;
      }

      // Otherwise, we have to move along the current hole boundary.
      if (_hole_iter != _face->holes_end())
      {
        // If we have not traversed the entire current hole (namely this is the
        // first increment operation, or we still have not completed a full
        // cycle around the current hole), move to the next halfedge along the
        // hole.
        if (! _curr_hole_incremented)
        {
          ++_curr_hole_circ;
          _curr_hole_incremented = true;
          return;
        }

        if (_curr_hole_circ != *_hole_iter)
        {
          ++_curr_hole_circ;
          
          if (_curr_hole_circ != *_hole_iter)
            return;
        }

        // If we reached here, we have to proceed to the next hole.
        ++_hole_iter;
        if (_hole_iter != _face->holes_end())
        {
          _curr_hole_circ = *_hole_iter;
          _curr_hole_incremented = false;
        }
        else
        {
          // In this case we finished traversing all outer and inner CCBs:
          _end = true;
        }
      }
      else
      {
        // In this case we finished traversing all outer and inner CCBs:
        _end = true;
      }
      
      return;
    }

  };

  // Data members:
  mutable Arrangement_2    *p_arr;    // The primal arrangement.
  
public:

  typedef Face_neighbor_iterator            Incident_edge_iterator;

  /*! Default constructor. */
  Dual() :
    p_arr(NULL)
  {}

  /*! Constructor from an arrangement. */
  Dual(const Arrangement_2& arr) :
    p_arr(const_cast<Arrangement_2 *>(&arr))
  {}

  /*! Get the number of vertices (face of the primal arrangement). */
  Size number_of_vertices() const
  {
    return p_arr->number_of_faces();
  }

  /*! Traverse the vertices (faces of the primal arrangement). */
  Vertex_iterator vertices_begin() const
  {
    return p_arr->faces_begin();
  }

  Vertex_iterator vertices_end() const
  {
    return p_arr->faces_end();
  }

  /*! Get the number of edges. */
  Size number_of_edges() const
  {
    return p_arr->number_of_halfedges();
  }

  /*! Traverse the edges. */
  Edge_iterator edges_begin() const
  {
    return p_arr->halfedges_begin();
  }

  Edge_iterator edges_end() const
  {
    return p_arr->halfedges_end();
  }

  /*!
   * Get the dual vertex-degree (number of edges forming the face boundary).
   */
  Size degree(Vertex_handle v) const
  {
    Incident_edge_iterator   begin = Incident_edge_iterator(v, true, true);
    Incident_edge_iterator   end = Incident_edge_iterator(v, false, true);
    Size                     deg = 0;

    while (begin != end)
    {
      deg++;
      ++begin;
    }

    return deg;
  }

  /*! Traverse the outgoing edges of a given vertex. */
  Incident_edge_iterator out_edges_begin(Vertex_handle v) const
  {
    return Incident_edge_iterator(v, true, true);
  }

  Incident_edge_iterator out_edges_end(Vertex_handle v) const
  {
    return Incident_edge_iterator(v, true, false);
  }

  /*! Traverse the ingoing edges of a given vertex. */
  Incident_edge_iterator in_edges_begin(Vertex_handle v) const
  {
    return Incident_edge_iterator(v, false, true);
  }

  Incident_edge_iterator in_edges_end(Vertex_handle v) const
  {
    return Incident_edge_iterator(v, false, false);
  }
};

CGAL_END_NAMESPACE

#include <boost/graph/graph_concepts.hpp>
#include <boost/iterator/counting_iterator.hpp>

namespace boost {

/*! \class
 * Specialization of the BGL graph-traits template, which serve as a dual
 * adapter for Arrangment_2, where the arrangement faces correspond to graph
 * verices, and two graph vertices are connected if the two corrsponding
 * faces are adjacent.
 * We consider the graph as directed. We also allow parallel edges, as two
 * faces may have more than one common edges.
 */
template <class Traits_, class Dcel_>
class graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >
{
public:
  
  typedef Traits_                              Traits_2;
  typedef Dcel_                                Dcel;
  typedef CGAL::Arrangement_2<Traits_2, Dcel>  Arrangement_2;
  typedef CGAL::Dual<Arrangement_2>            Dual_arr_2;

private:

  typedef typename Dual_arr_2::Vertex_iterator         Vertex_iterator;
  typedef typename Dual_arr_2::Edge_iterator           Edge_iterator;
  typedef typename Dual_arr_2::Incident_edge_iterator  Incident_edge_iterator;
 
  /*! \struct
   * Define the arrangement traversal category, which indicates the arrangement
   * models the BidirectionalGraph concept and the VertexListGraph and
   * EdgeListGraph concepts.
   */
  struct Dual_arr_traversal_category : 
    public virtual boost::bidirectional_graph_tag,   // This tag refines the
                                                     // incidence_graph_tag.
    public virtual boost::vertex_list_graph_tag,  // Can iterate over vertices.
    public virtual boost::edge_list_graph_tag     // Can iterate over edges.
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

}; // namespace boost

CGAL_BEGIN_NAMESPACE

// Functions required by the IncidenceGraph concept:
// -------------------------------------------------

/*!
 * Get the out-degree of a vertex in a given dual arrangement.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \param Number of halfedges around the boundary of the primal face.
 */
template <class Traits_, class Dcel_>
typename
boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
degree_size_type
out_degree(typename boost::graph_traits<CGAL::Dual<CGAL::
                   Arrangement_2<Traits_, Dcel_> > >::vertex_descriptor v,
           const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return darr.degree(v);
}

/*!
 * Return a range of the out-edges of a vertex given by its descriptor and the
 * dual arrangement it belongs to.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \return A pair of out-edges iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::out_edge_iterator,
          typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::out_edge_iterator>
out_edges(typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::vertex_descriptor v,
          const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return std::make_pair(darr.out_edges_begin(v), darr.out_edges_end(v));
}

/*!
 * Get the source vertex of a dual arrangement edge.
 * \param e The edge.
 * \param darr The dual arrangement.
 * \return The incident face of e in the primal arrangement.
 */
template <class Traits_, class Dcel_>
typename
boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
vertex_descriptor
source(typename boost::graph_traits<CGAL::Dual<CGAL::
         Arrangement_2<Traits_, Dcel_> > >::edge_descriptor e,
       const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& /* darr */)
{
  return e->face();
}

/*!
 * Get the target vertex of a dual arrangement edge.
 * \param e The edge.
 * \param darr The dual arrangement.
 * \return The incident face of e's twin in the primal arrangement.
 */
template <class Traits_, class Dcel_>
typename
boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
vertex_descriptor
target(typename boost::graph_traits<CGAL::Dual<CGAL::
         Arrangement_2<Traits_, Dcel_> > >::edge_descriptor e,
       const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& /* darr */)
{
  return e->twin()->face();
}

// Functions required by the BidirectionalGraph concept:
// -----------------------------------------------------

/*!
 * Get the in-degree of a vertex in a given dual arrangement.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \param Number of halfedges around the boundary of the primal face.
 */
template <class Traits_, class Dcel_>
typename
boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
degree_size_type
in_degree(typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::vertex_descriptor v,
          const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return darr.degree(v);
}

/*!
 * Return a range of the in-edges of a vertex given by its descriptor and the
 * dual arrangement it belongs to.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \return A pair of in-edges iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::in_edge_iterator,
          typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::in_edge_iterator>
in_edges(typename boost::graph_traits<CGAL::Dual<CGAL::
           Arrangement_2<Traits_, Dcel_> > >::vertex_descriptor v,
         const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return std::make_pair(darr.in_edges_begin(v), darr.in_edges_end(v));
}

/*!
 * Get the degree of a vertex in a given dual arrangement.
 * \param v The vertex.
 * \param darr The dual arrangement.
 * \param Number of ingoing and outgoing halfedges incident to v.
 */
template <class Traits_, class Dcel_>
typename
boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
degree_size_type
degree(typename boost::graph_traits<CGAL::Dual<CGAL::
         Arrangement_2<Traits_, Dcel_> > >::vertex_descriptor v,
       const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return (2 * darr.degree(v));
}

// Functions required by the VertexListGraph concept:
// --------------------------------------------------

/*!
 * Get the number of vertices in the given dual arrangement. 
 * \param darr The dual arrangement.
 * \return Number of faces in the primal arrangement.
 */
template <class Traits_, class Dcel_>
typename
boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
vertices_size_type
num_vertices(const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return darr.number_of_vertices();
}

/*!
 * Get the range of vertices of the given dual arrangement.
 * \param darr The dual arrangement.
 * \return A pair of vertex iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::vertex_iterator,
          typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::vertex_iterator>
vertices(const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  typedef typename
    boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
                                              vertex_iterator vertex_iterator;
  return std::pair<vertex_iterator, vertex_iterator>(darr.vertices_begin(),
                                                     darr.vertices_end());
}

// Functions required by the EdgeListGraph concept:
// ------------------------------------------------

/*!
 * Get the number of edges in the given dual arrangement. 
 * \param darr The dual arrangement.
 * \return Number of halfedges in the primal arrangement.
 */
template <class Traits_, class Dcel_>
typename
boost::graph_traits<CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> > >::
edges_size_type
num_edges(const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return darr.number_of_edges(); 
}

/*!
 * Get the range of edges of the given dual arrangement.
 * \param darr The dual arrangement.
 * \return A pair of edge iterators.
 */
template <class Traits_, class Dcel_>
std::pair<typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::edge_iterator,
          typename boost::graph_traits<CGAL::Dual<CGAL::
            Arrangement_2<Traits_, Dcel_> > >::edge_iterator>
edges(const CGAL::Dual<CGAL::Arrangement_2<Traits_, Dcel_> >& darr)
{
  return std::make_pair(darr.edges_begin(), darr.edges_end());
}

CGAL_END_NAMESPACE

#endif
