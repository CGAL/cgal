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
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Ophir Setter <ophirset@post.tau.ac.il>

/*! \file
 * This file contains the BGL dual adaptor for Arrangment_2, where the
 * arrangement faces are considered as the graph verticies. The file also
 * includes a generator function for the adaptor.
 */

#ifndef ARR_BGL_DUAL_ADAPTOR_H
#define ARR_BGL_DUAL_ADAPTOR_H

#include <boost/graph/graph_concepts.hpp>
#include <boost/iterator/counting_iterator.hpp>

CGAL_BEGIN_NAMESPACE

/************************************************************************
* Template functions and classes
*************************************************************************/
template <typename T_Arrangement_2>
class Arr_bgl_dual_adator
{
public:
  typedef T_Arrangement_2                       Arr_2;

  /////////////////////////////////////////////////////////////////////////
  // Requirments for GraphConcept
  /////////////////////////////////////////////////////////////////////////
  // Face_handle will be used as a vertex_descriptor since we'd like the
  // faces of the planar map to be the vertecies in the graph.
  typedef typename Arr_2::Face_handle           vertex_descriptor;

  // Halfedge_handle will be used as an edge_descriptor since it is the 
  // border between one face and another face.
  typedef typename Arr_2::Halfedge_handle       edge_descriptor;

  // Since in CGAL planar map there is no importance to directed/not directed
  // our graph will be directed and will allow parallel edges (This fact will
  // make our life easier).
  typedef boost::directed_tag                   directed_category;
  typedef boost::allow_parallel_edge_tag        edge_parallel_category;

  // The struct is needed for the traversal_category since our adaptor supports:
  // 1) Incidence Graph
  // 2) Vertex List Graph
  struct planar_map_traversal_category : 
    public virtual boost::incidence_graph_tag,
    public virtual boost::vertex_list_graph_tag { };
  typedef planar_map_traversal_category         traversal_category;

  /////////////////////////////////////////////////////////////////////////
  // Requirments for IncidenceGraphConcept
  /////////////////////////////////////////////////////////////////////////
  // The neighbors of a face are the faces on the borders of the outer boundary,
  // and its holes. That means for that the out_edge_iterator we need an
  // iterator that passes on all the outer boundary half-edges and then on one
  // half edge from each hole. This is exactly what the following class does.
  // For explaination on iterator_facade, check the BOOST documentation (It
  // helps us to create iterator classes).
  class Face_neighbors_iterator :
    public boost::iterator_facade<Face_neighbors_iterator, 
                                  typename Arr_2::Ccb_halfedge_circulator,
                                  boost::forward_traversal_tag,
                                  typename Arr_2::Ccb_halfedge_circulator>
  {
  public:
    typedef typename Arr_2::Face_handle                 Face_handle;
    typedef typename Arr_2::Ccb_halfedge_circulator     Ccb_halfedge_circulator;
    typedef typename Arr_2::Holes_iterator              Holes_iterator;

    Face_neighbors_iterator( ) {};

    /*!
     * Specialize the get property_map function, so that BOOST is able to use
     * the above classes.
     * @param face The face that the iterator reffers to.
     * @param bStart Whether the iterator is "begin" or "end".
     */
    Face_neighbors_iterator( Face_handle face, bool bStart = true )
    {
      m_face = face;
      m_has_outer_ccb = !(*face).is_unbounded();
      m_ccb_incremented = true;
      if( m_has_outer_ccb )
      {
        m_outer_ccb_circ = (*face).outer_ccb();
        m_ccb_incremented = !bStart;
      }

      if( bStart  )
      {
        m_hole_iter = (*face).holes_begin();
        m_current_hole_incremented = true;
        if( m_hole_iter!=(*face).holes_end() )
        {
          m_current_hole_circ = *m_hole_iter;
          m_current_hole_incremented = false;
        }
      }
      else // end iterator.
      {
        m_hole_iter = (*face).holes_end();
      }
    }  

  private:
    // In the order of the iterator_facade be able to access our private member
    // functions.
    friend class boost::iterator_core_access;

    bool equal(const Face_neighbors_iterator & it) const
    {
      return (m_has_outer_ccb == it.m_has_outer_ccb) &&
        (m_ccb_incremented == it.m_ccb_incremented) &&
        (m_outer_ccb_circ == it.m_outer_ccb_circ) && 
        (m_hole_iter == it.m_hole_iter) &&
        (m_face == it.m_face);
    }

    // Function name   :  dereference - Derefernces the iterator.
    Ccb_halfedge_circulator dereference() const
    {
      if( m_has_outer_ccb )
      {
        if( !m_ccb_incremented )
        {
          return m_outer_ccb_circ;
        }

        if( m_outer_ccb_circ != (*m_face).outer_ccb() )
        {
          return m_outer_ccb_circ;
        }
      }

      return m_current_hole_circ;
    }

    // Function name   :  increment - Increments of the iterator.
    void increment()
    {
      if( m_has_outer_ccb )
      {
        if( !m_ccb_incremented )
        {
          m_ccb_incremented = true;
          ++m_outer_ccb_circ;
          return;
        }

        if( m_outer_ccb_circ != (*m_face).outer_ccb() )
        {
          ++m_outer_ccb_circ;
          return;
        }
      }

      if( m_hole_iter!=(*m_face).holes_end() )
      {
        if( !m_current_hole_incremented )
        {
          m_current_hole_incremented = true;
          ++m_current_hole_circ;
          return;
        }

        if( m_current_hole_circ != (*m_hole_iter) )
        {
          ++m_current_hole_circ;
          return;
        }

        ++m_hole_iter;
        if( m_hole_iter!=(*m_face).holes_end() )
        {
          m_current_hole_circ = (*m_hole_iter);
          m_current_hole_incremented = false;
        }
      }
    }

  private:
    bool m_has_outer_ccb;
    bool m_ccb_incremented;
    Ccb_halfedge_circulator m_outer_ccb_circ;
    Holes_iterator m_hole_iter;
    Ccb_halfedge_circulator m_current_hole_circ;
    bool m_current_hole_incremented;
    Face_handle m_face;
  };
  typedef Face_neighbors_iterator                       out_edge_iterator;

  // degree_size_type will be int.
  typedef int                                           degree_size_type;

  /////////////////////////////////////////////////////////////////////////
  // Requirments for VertexListGraphConcept
  /////////////////////////////////////////////////////////////////////////
  // because vertex_descriptor is a Face_handle, and because iterators are 
  // assignable to handles, we can use the operator * to return the actual 
  // iterator, then the Vertex_iterator of the planar map can be used as a 
  // boost vertex iterator. 
  // This is exactly what boost::counting_iterator does (return the actual
  // object in operator *)
  typedef typename Arr_2::Face_iterator                 Face_iterator;
  typedef boost::counting_iterator<Face_iterator>       vertex_iterator;

  // vertices_size_type will be int.
  typedef int                           vertices_size_type;

  /////////////////////////////////////////////////////////////////////////
  // void for non interesting types
  /////////////////////////////////////////////////////////////////////////
  typedef void                          adjacency_iterator;
  typedef void                          in_edge_iterator;
  typedef void                          edge_iterator;
  typedef void                          edges_size_type;

  /////////////////////////////////////////////////////////////////////////
  // Constructor and delegation functions
  /////////////////////////////////////////////////////////////////////////
  Arr_bgl_dual_adator(Arr_2 & in_planar_map) : m_planar_map(in_planar_map) {}

  Face_iterator faces_begin() {return m_planar_map.faces_begin();}
  Face_iterator faces_end()   {return m_planar_map.faces_end();}

  typedef typename Arr_2::Size          Size;
  Size number_of_faces() const {return m_planar_map.number_of_faces();}

private:
  Arr_2 & m_planar_map;
};


/////////////////////////////////////////////////////////////////////////
// Functions for IncidenceGraphConcept
/////////////////////////////////////////////////////////////////////////

/*!
 * The function returns a range of the out-edges of a vertex given by its
 * descriptor and the planar map graph it is part of. In this case the vertex
 * is a planar map face.
 * @param in_vertex The vertex whos out-edges will be returned.
 * @param in_graph The planar map graph of the vertex.
 * @return iterator range of out-edges iterators.
 */
template <typename T_Arrangement_2>
  std::pair<
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
        out_edge_iterator,
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
        out_edge_iterator>
  out_edges( 
    typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
        vertex_descriptor in_vertex, 
    const Arr_bgl_dual_adator<T_Arrangement_2> & in_graph)
{
  // typedefs for making the next lines smaller and READABLE.
  typedef typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
    out_edge_iterator
    out_edge_iterator;

  // return the Face_neighbors_iterator of the face.
  out_edge_iterator itBegin = out_edge_iterator(in_vertex, true);
  out_edge_iterator itEnd = out_edge_iterator(in_vertex, false);

  return std::make_pair( itBegin, itEnd );
}

/*!
 * The function returns the out-degree of a vertex given by its descriptor
 * and the planar map graph it is part of. The out-degree of a vertex is the
 * number of the out-edges of that vertex.
 * @param in_vertex The vertex whos out-edges will be returned.
 * @param in_graph The planar map graph of the vertex.
 * @return Number of out-edges.
 */
template <typename T_Arrangement_2>
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
        degree_size_type
  out_degree(
    typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
      vertex_descriptor in_vertex, 
    const Arr_bgl_dual_adator<T_Arrangement_2> & in_graph)
{
  typedef typename
    boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::out_edge_iterator
    out_edge_iterator;

  std::pair<out_edge_iterator,out_edge_iterator> edges =
    out_edges(in_vertex, in_graph);
  return std::distance(edges.first, edges.second); 
}

/*!
 * The function returns source vertex of an edge given by  descriptor
 * and the planar map graph the edge is part of.
 * @param in_edge The edge whos source will be returned.
 * @param in_graph The planar map graph of the edge.
 * @return The source.
 */
template <typename T_Arrangement_2>
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
        vertex_descriptor
  source(
    typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
      edge_descriptor in_edge,
    const Arr_bgl_dual_adator<T_Arrangement_2> & in_graph)
{
  // the source face is face of the halfedge (arbitrary).
  return (*in_edge).face();
}

/*!
 * The function returns the target vertex of an edge given by its descriptor
 * and the planar map graph the edge is part of.
 * @param in_edge The edge whos target will be returned.
 * @param in_graph The planar map graph of the edge.
 * @return The target.
 */
template <typename T_Arrangement_2>
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
    vertex_descriptor
  target(
    typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
      edge_descriptor in_edge,
    const Arr_bgl_dual_adator<T_Arrangement_2> & in_graph)
{
  // the target face is face of the opposite halfedge.
  return (*((*in_edge).twin())).face();
}



/////////////////////////////////////////////////////////////////////////
// Functions for VertexListGraphConcept
/////////////////////////////////////////////////////////////////////////

/*!
  The function returns for a given planar map graph, 
  an iterator range of iterator to all vertices. 
  @param in_graph The planar map graph.
  @return Iterator range of the planar map faces.
*/
template <typename T_Arrangement_2>
  std::pair<
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
    vertex_iterator,
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
    vertex_iterator>
  vertices(const Arr_bgl_dual_adator<T_Arrangement_2> & in_graph)
{
  typedef typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
    vertex_descriptor    vertex_descriptor;

  // we need the const cast since vertex_descriptor is non-const, but the given
  // graph needs to
  // be const.
  Arr_bgl_dual_adator<T_Arrangement_2> & graph =
    const_cast<Arr_bgl_dual_adator<T_Arrangement_2>&>(in_graph);
  return std::make_pair( vertex_descriptor(graph.faces_begin()), 
    vertex_descriptor(graph.faces_end()) );
}

/*!
 * The function returns the number of vertices of a given graph, 
 * @param in_graph The planar map graph.
 * @return Number of vertices.
 */
template <typename T_Arrangement_2>
  typename boost::graph_traits<Arr_bgl_dual_adator<T_Arrangement_2> >::
    vertices_size_type
  num_vertices(const Arr_bgl_dual_adator<T_Arrangement_2> & in_graph)
{
  return in_graph.number_of_faces(); 
}

/////////////////////////////////////////////////////////////////////////
// Generate function for the adaptor
/////////////////////////////////////////////////////////////////////////
template <typename T_Arrangement_2>
Arr_bgl_dual_adator<T_Arrangement_2>
make_arr_face_graph_adaptor(T_Arrangement_2 & in_planar_map)
{
  return Arr_bgl_dual_adator<T_Arrangement_2>(in_planar_map);
}

CGAL_END_NAMESPACE

#endif
