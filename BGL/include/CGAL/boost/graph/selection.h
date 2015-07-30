// Copyright (c) 2015  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_BOOST_GRAPH_KTH_SIMPLICIAL_NEIGHBORHOOD_H
#define CGAL_BOOST_GRAPH_KTH_SIMPLICIAL_NEIGHBORHOOD_H

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {


/// Operation on faces
namespace internal{
// extract edges in non-selected faces (boundary excluded but one)
template <class FaceRange, class FaceGraph, class IsFaceSelectedPMap, class OutputIterator>
OutputIterator
extract_selection_boundary(
  FaceRange& face_range,
  FaceGraph& graph,
  IsFaceSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<FaceGraph> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  BOOST_FOREACH(face_descriptor fd, face_range)
  {
    BOOST_FOREACH(  halfedge_descriptor h,
                    halfedges_around_face(halfedge(fd, graph), graph) )
    {
      halfedge_descriptor opp_hd = opposite(h, graph);
      face_descriptor opp_fd = face( opp_hd, graph );
      if (opp_fd!=GT::null_face())
      {
        if ( !get(is_selected, opp_fd) )
          *out++=opp_hd;
      }
      else{
        opp_hd=opposite( next( opp_hd, graph), graph );
        if ( !get( is_selected, face(opp_hd, graph) ) )
          *out++=opp_hd;
      }
    }
  }
  return out;
}
} //end of namespace internal

template <class FaceRange, class FaceGraph, class IsFaceSelectedPMap, class OutputIterator>
OutputIterator
dilate_face_selection(
  const FaceRange& selection,
  FaceGraph& graph,
  unsigned int k,
  IsFaceSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<FaceGraph> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  std::vector<face_descriptor> current_selection(selection.begin(), selection.end());
  for (unsigned int i=0; i<k; ++i)
  {
    //extract faces on the boundary of the selection
    std::vector<halfedge_descriptor> selection_boundary_halfedges;
    internal::extract_selection_boundary(current_selection, graph, is_selected,
                                         std::back_inserter(selection_boundary_halfedges));

    if (selection_boundary_halfedges.empty()) break;

    //collect faces around the target vertex of the selection boundary halfedges
    std::set<face_descriptor> new_selection_set;
    BOOST_FOREACH(halfedge_descriptor hd, selection_boundary_halfedges)
    {
      face_descriptor fd=face(hd, graph);
      while( !get(is_selected,fd) )
      {
        new_selection_set.insert(fd);
        hd=opposite( next(hd, graph), graph );
        fd=face(hd, graph);
        if ( face(hd, graph)==GT::null_face() ) break;
      }
    }

    // extract unique selection
    std::vector<face_descriptor> new_selection;
    BOOST_FOREACH(face_descriptor fd, new_selection_set)
    {
      *out++=fd;
      new_selection.push_back(fd);
      put( is_selected, fd, true );
    }
    current_selection.swap(new_selection);
  }
  return out;
}

template <class FaceRange, class FaceGraph, class IsFaceSelectedPMap, class OutputIterator>
OutputIterator
erode_face_selection(
  const FaceRange& selection,
  FaceGraph& graph,
  unsigned int k,
  IsFaceSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<FaceGraph> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  std::vector<face_descriptor> current_selection(selection.begin(), selection.end());
  for (unsigned int i=0; i<k; ++i)
  {
    //extract faces on the boundary of the selection
    std::vector<halfedge_descriptor> selection_boundary_halfedges;
    internal::extract_selection_boundary(current_selection, graph, is_selected,
                                         std::back_inserter(selection_boundary_halfedges));

    if (selection_boundary_halfedges.empty()) break;


    //collect faces around the target vertex of the selection boundary halfedges
    std::set<face_descriptor> elements_to_remove;
    BOOST_FOREACH(halfedge_descriptor hd, selection_boundary_halfedges)
    {
      hd = opposite(hd, graph);
      face_descriptor fd=face( hd, graph );
      while( face(hd, graph)!=GT::null_face() && get(is_selected,fd) )
      {
        elements_to_remove.insert(fd);
        hd=opposite( next(hd, graph), graph );
        fd=face(hd, graph);
      }
    }

    /// update is-selected attribute and output iterator
    BOOST_FOREACH(face_descriptor fd, elements_to_remove)
    {
      *out++=fd;
      put( is_selected, fd, false );
    }

    // update the set of currently selected faces
    std::vector<face_descriptor> new_selection;
    BOOST_FOREACH(face_descriptor fd, current_selection)
      if ( !elements_to_remove.count(fd) )
        new_selection.push_back(fd);
    current_selection.swap(new_selection);
  }
  return out;
}

// select all faces incident to the target vertex of halfedges in `hedges`
template <class HalfedgeRange, class FaceGraph, class OutputIterator>
OutputIterator
select_incident_faces(
  const HalfedgeRange& hedges,
  FaceGraph& graph,
  OutputIterator out)
{
  typedef boost::graph_traits<FaceGraph> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  //collect faces around the target vertex of the selection boundary halfedges
  std::set<face_descriptor> selection_set;
  BOOST_FOREACH(halfedge_descriptor hd, hedges)
  {
    halfedge_descriptor first = hd;
    face_descriptor fd=face(hd, graph);
    do
    {
      if ( face(hd, graph)!=GT::null_face() && selection_set.insert(fd).second)
        *out++=fd;
      hd=opposite( next(hd, graph), graph );
      fd=face(hd, graph);
    }while( hd!=first );
  }

  return out;
}

/// Operations on edges
template <class EdgeRange, class HalfedgeGraph, class IsEdgeSelectedPMap, class OutputIterator>
OutputIterator
dilate_edge_selection(
  const EdgeRange& selection,
  HalfedgeGraph& graph,
  unsigned int k,
  IsEdgeSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<HalfedgeGraph> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  std::vector<edge_descriptor> current_selection(selection.begin(), selection.end());
  for (unsigned int i=0; i<k; ++i)
  {
    if (current_selection.empty()) break;

    //collect adjacent edges not already selected
    std::set<edge_descriptor> new_selection_set;
    BOOST_FOREACH(edge_descriptor ed, current_selection)
    {
      halfedge_descriptor hdi=halfedge(ed,graph);
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_source( hdi, graph))
      {
        edge_descriptor ned=edge(hd, graph);
        if (!get(is_selected, ned)) new_selection_set.insert(ned);
      }
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target( hdi, graph))
      {
        edge_descriptor ned=edge(hd, graph);
        if (!get(is_selected, ned)) new_selection_set.insert(ned);
      }
    }

    // extract unique selection
    std::vector<edge_descriptor> new_selection;
    BOOST_FOREACH(edge_descriptor ed, new_selection_set)
    {
      *out++=ed;
      new_selection.push_back(ed);
      put( is_selected, ed, true );
    }
    current_selection.swap(new_selection);
  }
  return out;
}

template <class EdgeRange, class HalfedgeGraph, class IsEdgeSelectedPMap, class OutputIterator>
OutputIterator
erode_edge_selection(
  const EdgeRange& selection ,
  HalfedgeGraph& graph,
  unsigned int k,
  IsEdgeSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<HalfedgeGraph> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  // extract the set of vertices on the border
  std::set<vertex_descriptor> unique_vertex_set;
  BOOST_FOREACH(edge_descriptor ed, selection)
  {
    halfedge_descriptor hd=halfedge(ed,graph);
    BOOST_FOREACH(halfedge_descriptor nhd, halfedges_around_source( hd, graph))
    {
      edge_descriptor ned=edge(nhd, graph);
      if (!get(is_selected, ned)){
        unique_vertex_set.insert(source(hd,graph));
        break;
      }
    }
    BOOST_FOREACH(halfedge_descriptor nhd, halfedges_around_target( hd, graph))
    {
      edge_descriptor ned=edge(nhd, graph);
      if (!get(is_selected, ned)){
        unique_vertex_set.insert(target(hd,graph));
        break;
      }
    }
  }

  std::vector<vertex_descriptor> current_selection_border(unique_vertex_set.begin(), unique_vertex_set.end());
  for (unsigned int i=0; i<k; ++i)
  {
    if (current_selection_border.empty()) break;

    //collect incident edges selected
    std::set<edge_descriptor> edges_to_deselect;
    unique_vertex_set.clear();
    BOOST_FOREACH(vertex_descriptor vd, current_selection_border)
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target( halfedge(vd,graph), graph))
      {
        edge_descriptor ed = edge(hd, graph);
        if (get(is_selected, ed)){
          edges_to_deselect.insert(ed);
          unique_vertex_set.insert(source(hd, graph));
        }
      }

    // extract unique selection
    BOOST_FOREACH(edge_descriptor ed, edges_to_deselect)
    {
      *out++=ed;
      put( is_selected, ed, false );
    }

    current_selection_border.assign(unique_vertex_set.begin(), unique_vertex_set.end());
  }
  return out;
}

/// Operations on vertices
template <class VertexRange, class HalfedgeGraph, class IsVertexSelectedPMap, class OutputIterator>
OutputIterator
dilate_vertex_selection(
  const VertexRange& selection,
  HalfedgeGraph& graph,
  unsigned int k,
  IsVertexSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<HalfedgeGraph> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  std::vector<vertex_descriptor> current_selection(selection.begin(), selection.end());
  for (unsigned int i=0; i<k; ++i)
  {
    if (current_selection.empty()) break;

    //collect adjacent vertices not already selected
    std::set<vertex_descriptor> new_selection_set;
    BOOST_FOREACH(vertex_descriptor vd, current_selection)
      BOOST_FOREACH(vertex_descriptor nvd, vertices_around_target( halfedge(vd,graph), graph))
        if (!get(is_selected, nvd)) new_selection_set.insert(nvd);

    // extract unique selection
    std::vector<vertex_descriptor> new_selection;
    BOOST_FOREACH(vertex_descriptor vd, new_selection_set)
    {
      *out++=vd;
      new_selection.push_back(vd);
      put( is_selected, vd, true );
    }
    current_selection.swap(new_selection);
  }
  return out;
}

template <class VertexRange, class HalfedgeGraph, class IsVertexSelectedPMap, class OutputIterator>
OutputIterator
erode_vertex_selection(
  const VertexRange& selection,
  HalfedgeGraph& graph,
  unsigned int k,
  IsVertexSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<HalfedgeGraph> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  // collect vertices incident to a selected one
  std::set<vertex_descriptor> unique_vertex_set;
  BOOST_FOREACH(vertex_descriptor vd, selection)
    BOOST_FOREACH(vertex_descriptor nvd, vertices_around_target( halfedge(vd,graph), graph))
        if (!get(is_selected, nvd)) unique_vertex_set.insert(nvd);

  std::vector<vertex_descriptor> current_selection_border(unique_vertex_set.begin(), unique_vertex_set.end());
  for (unsigned int i=0; i<k; ++i)
  {
    if (current_selection_border.empty()) break;

    //collect adjacent vertices selected
    std::set<vertex_descriptor> vertices_to_deselect;
    BOOST_FOREACH(vertex_descriptor vd, current_selection_border)
      BOOST_FOREACH(vertex_descriptor nvd, vertices_around_target( halfedge(vd,graph), graph))
        if (get(is_selected, nvd)) vertices_to_deselect.insert(nvd);

    // extract unique selection
    std::vector<vertex_descriptor> new_selection_border;
    BOOST_FOREACH(vertex_descriptor vd, vertices_to_deselect)
    {
      *out++=vd;
      new_selection_border.push_back(vd);
      put( is_selected, vd, false );
    }
    current_selection_border.swap(new_selection_border);
  }
  return out;
}

} //end of namespace CGAL

#endif //CGAL_BOOST_GRAPH_KTH_SIMPLICIAL_NEIGHBORHOOD_H
