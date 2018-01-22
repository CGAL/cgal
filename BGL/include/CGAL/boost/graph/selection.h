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


// Operation on faces
namespace internal{
// extract edges in non-selected faces (boundary excluded but one)
template <class FaceRange, class FaceGraph, class IsFaceSelectedPMap, class OutputIterator>
OutputIterator
extract_selection_boundary(
  FaceRange& face_range,
  FaceGraph& fg,
  IsFaceSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<FaceGraph> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  BOOST_FOREACH(face_descriptor fd, face_range)
  {
    BOOST_FOREACH(  halfedge_descriptor h,
                    halfedges_around_face(halfedge(fd, fg), fg) )
    {
      halfedge_descriptor opp_hd = opposite(h, fg);
      face_descriptor opp_fd = face( opp_hd, fg );
      if (opp_fd!=GT::null_face())
      {
        if ( !get(is_selected, opp_fd) )
          *out++=opp_hd;
      }
      else{
        opp_hd=opposite( next( opp_hd, fg), fg );
        if ( !get( is_selected, face(opp_hd, fg) ) )
          *out++=opp_hd;
      }
    }
  }
  return out;
}
} //end of namespace internal


/*!
\ingroup PkgBGLSelectionFct
Augments a selection with faces of `fg` that are adjacent
to a face in `selection`. This process is applied `k` times considering
all faces added in the previous steps.
Two faces are said to be adjacent if they share a vertex or an edge.
Each new face added in the selection is added exactly once in `out`.
\tparam FaceRange a range of face descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
\tparam FaceGraph a model of `FaceGraph`.
\tparam IsFaceSelectedPMap a model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%face_descriptor`
        as key type and `bool` as value type.
\tparam OutputIterator an output iterator accepting face descriptors.
\param selection the initial selection of faces that will be expanded.
\param fg the graph containing the selected faces.
\param k the number of times the expansion procedure is iteratively applied.
\param is_selected indicates if a face is part of the selection. It is updated by the function
       to accomodate new faces added to the selection.
\param out new faces added to the selection are added exactly once in `out`.
*/
template <class FaceRange, class FaceGraph, class IsFaceSelectedPMap, class OutputIterator>
OutputIterator
expand_face_selection(
  const FaceRange& selection,
  FaceGraph& fg,
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
    internal::extract_selection_boundary(current_selection, fg, is_selected,
                                         std::back_inserter(selection_boundary_halfedges));

    if (selection_boundary_halfedges.empty()) break;

    //collect faces around the target vertex of the selection boundary halfedges
    std::set<face_descriptor> new_selection_set;
    BOOST_FOREACH(halfedge_descriptor hd, selection_boundary_halfedges)
    {
      face_descriptor fd=face(hd, fg);
      while( !get(is_selected,fd) )
      {
        new_selection_set.insert(fd);
        hd=opposite( next(hd, fg), fg );
        fd=face(hd, fg);
        if ( face(hd, fg)==GT::null_face() ) break;
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

/*!
\ingroup PkgBGLSelectionFct
Diminishes a selection of faces from faces adjacent to a non-selected face.
This process is applied `k` times considering all faces removed in the previous steps.
Two faces are said to be adjacent if they share a vertex or an edge.
Each face removed from the selection is added exactly once in `out`.
\tparam FaceRange a range of face descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
\tparam FaceGraph a model of `FaceGraph`.
\tparam IsFaceSelectedPMap a model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%face_descriptor`
        as key type and `bool` as value type.
\tparam OutputIterator an output iterator accepting face descriptors.
\param selection the initial selection of faces that will be expanded.
\param fg the graph containing the selected faces.
\param k the number of times the reduction procedure is iteratively applied.
\param is_selected indicates if a face is part of the selection. It is updated by the function
       to accomodate faces removed from the selection.
\param out faces removed from the selection are added exactly once in `out`.
*/
template <class FaceRange, class FaceGraph, class IsFaceSelectedPMap, class OutputIterator>
OutputIterator
reduce_face_selection(
  const FaceRange& selection,
  FaceGraph& fg,
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
    internal::extract_selection_boundary(current_selection, fg, is_selected,
                                         std::back_inserter(selection_boundary_halfedges));

    if (selection_boundary_halfedges.empty()) break;


    //collect faces around the target vertex of the selection boundary halfedges
    std::set<face_descriptor> elements_to_remove;
    BOOST_FOREACH(halfedge_descriptor hd, selection_boundary_halfedges)
    {
      hd = opposite(hd, fg);
      face_descriptor fd=face( hd, fg );
      while( face(hd, fg)!=GT::null_face() && get(is_selected,fd) )
      {
        elements_to_remove.insert(fd);
        hd=opposite( next(hd, fg), fg );
        fd=face(hd, fg);
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


/*!
\ingroup PkgBGLSelectionFct
Discovers and puts in `out` all faces incident to the target vertex
of a halfedge in `hedges`. Faces are put exactly once in `out`.
\tparam HalfedgeRange a range of halfedge descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
\tparam HalfedgeGraph a model of `HalfedgeGraph`.
\tparam OutputIterator an output iterator accepting face descriptors.
\param hedges the range a halfedge descriptors consider during the face selection.
\param fg the graph containing the input halfedges.
\param out faces added to the selection are added exactly once in `out`.
*/
template <class HalfedgeRange, class FaceGraph, class OutputIterator>
OutputIterator
select_incident_faces(
  const HalfedgeRange& hedges,
  FaceGraph& fg,
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
    face_descriptor fd=face(hd, fg);
    do
    {
      if ( face(hd, fg)!=GT::null_face() && selection_set.insert(fd).second)
        *out++=fd;
      hd=opposite( next(hd, fg), fg );
      fd=face(hd, fg);
    }while( hd!=first );
  }

  return out;
}

/*!
\ingroup PkgBGLSelectionFct
Augments a selection with edges of `fg` that are adjacent
to an edge in `selection`. This process is applied `k` times considering
all edges added in the previous steps.
Two edges are said to be adjacent if they are incident to the same face or vertex.
Each new edge added in the selection is added exactly once in `out`.
\tparam EdgeRange a range of edge descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
\tparam FaceGraph a model of `FaceGraph`.
\tparam IsEdgeSelectedPMap a model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%edge_descriptor`
        as key type and `bool` as value type.
\tparam OutputIterator an output iterator accepting edge descriptors.
\param selection the initial selection of edges that will be expanded.
\param fg the graph containing the selected edges.
\param k the number of times the expansion procedure is iteratively applied.
\param is_selected indicates if an edge is part of the selection. It is updated by the function
       to accomodate new edges added to the selection.
\param out new edges added to the selection are added exactly once in `out`.
*/
template <class EdgeRange, class HalfedgeGraph, class IsEdgeSelectedPMap, class OutputIterator>
OutputIterator
expand_edge_selection(
  const EdgeRange& selection,
  HalfedgeGraph& fg,
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
      halfedge_descriptor hdi=halfedge(ed,fg);
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_source( hdi, fg))
      {
        edge_descriptor ned=edge(hd, fg);
        if (!get(is_selected, ned)) new_selection_set.insert(ned);
      }
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target( hdi, fg))
      {
        edge_descriptor ned=edge(hd, fg);
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

/*!
\ingroup PkgBGLSelectionFct
Diminishes a selection of edges from edges adjacent to a non-selected edge.
This process is applied `k` times considering all edges removed in the previous steps.
Two edges are said to be adjacent if they are incident to the same face or vertex.
Each edge removed from the selection is added exactly once in `out`.
\tparam EdgeRange a range of edge descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
\tparam FaceGraph a model of `FaceGraph`.
\tparam IsEdgeSelectedPMap a model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%edge_descriptor`
        as key type and `bool` as value type.
\tparam OutputIterator an output iterator accepting edge descriptors.
\param selection the initial selection of edges that will be reduced.
\param fg the graph containing the selected edges.
\param k the number of times the reduction procedure is iteratively applied.
\param is_selected indicates if an edge is part of the selection. It is updated by the function
       to accomodate edges removed from the selection.
\param out edges removed from the selection are added exactly once in `out`.
*/
template <class EdgeRange, class HalfedgeGraph, class IsEdgeSelectedPMap, class OutputIterator>
OutputIterator
reduce_edge_selection(
  const EdgeRange& selection ,
  HalfedgeGraph& fg,
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
    halfedge_descriptor hd=halfedge(ed,fg);
    BOOST_FOREACH(halfedge_descriptor nhd, halfedges_around_source( hd, fg))
    {
      edge_descriptor ned=edge(nhd, fg);
      if (!get(is_selected, ned)){
        unique_vertex_set.insert(source(hd,fg));
        break;
      }
    }
    BOOST_FOREACH(halfedge_descriptor nhd, halfedges_around_target( hd, fg))
    {
      edge_descriptor ned=edge(nhd, fg);
      if (!get(is_selected, ned)){
        unique_vertex_set.insert(target(hd,fg));
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
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target( halfedge(vd,fg), fg))
      {
        edge_descriptor ed = edge(hd, fg);
        if (get(is_selected, ed)){
          edges_to_deselect.insert(ed);
          unique_vertex_set.insert(source(hd, fg));
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

/*!
\ingroup PkgBGLSelectionFct
Augments a selection with vertices of `fg` that are adjacent
to a vertex in `selection`. This process is applied `k` times considering
all vertices added in the previous steps.
Two vertices are said to be adjacent if they are part of the same face.
Each new vertex added in the selection is added exactly once in `out`.
\tparam VertexRange a range of vertex descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
\tparam FaceGraph a model of `FaceGraph`.
\tparam IsVertexSelectedPMap a model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
        as key type and `bool` as value type.
\tparam OutputIterator an output iterator accepting vertex descriptors.
\param selection the initial selection of vertices that will be expanded.
\param fg the graph containing the selected vertices.
\param k the number of times the expansion procedure is iteratively applied.
\param is_selected indicates if a vertex is part of the selection. It is updated by the function
       to accomodate new vertices added to the selection.
\param out new vertices added to the selection are added exactly once in `out`.
*/
template <class VertexRange, class HalfedgeGraph, class IsVertexSelectedPMap, class OutputIterator>
OutputIterator
expand_vertex_selection(
  const VertexRange& selection,
  HalfedgeGraph& fg,
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
      BOOST_FOREACH(vertex_descriptor nvd, vertices_around_target( halfedge(vd,fg), fg))
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

/*!
\ingroup PkgBGLSelectionFct
Diminishes a selection of vertices from vertices adjacent to a non-selected vertex.
This process is applied `k` times considering all vertices removed in the previous steps.
Two vertices are said to be adjacent if they are part of the same face.
Each vertex removed from the selection is added exactly once in `out`.
\tparam VertexRange a range of vertex descriptors, model of `Range`.
          Its iterator type is `InputIterator`.
\tparam FaceGraph a model of `FaceGraph`.
\tparam IsVertexSelectedPMap a model of `ReadWritePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
        as key type and `bool` as value type.
\tparam OutputIterator an output iterator accepting vertex descriptors.
\param selection the initial selection of vertices that will be reduced.
\param fg the graph containing the selected vertices.
\param k the number of times the reduction procedure is iteratively applied.
\param is_selected indicates if a vertex is part of the selection. It is updated by the function
       to accomodate vertices removed from the selection.
\param out vertices removed from the selection are added exactly once in `out`.
*/
template <class VertexRange, class HalfedgeGraph, class IsVertexSelectedPMap, class OutputIterator>
OutputIterator
reduce_vertex_selection(
  const VertexRange& selection,
  HalfedgeGraph& fg,
  unsigned int k,
  IsVertexSelectedPMap is_selected,
  OutputIterator out)
{
  typedef boost::graph_traits<HalfedgeGraph> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  // collect vertices incident to a selected one
  std::set<vertex_descriptor> unique_vertex_set;
  BOOST_FOREACH(vertex_descriptor vd, selection)
    BOOST_FOREACH(vertex_descriptor nvd, vertices_around_target( halfedge(vd,fg), fg))
        if (!get(is_selected, nvd)) unique_vertex_set.insert(nvd);

  std::vector<vertex_descriptor> current_selection_border(unique_vertex_set.begin(), unique_vertex_set.end());
  for (unsigned int i=0; i<k; ++i)
  {
    if (current_selection_border.empty()) break;

    //collect adjacent vertices selected
    std::set<vertex_descriptor> vertices_to_deselect;
    BOOST_FOREACH(vertex_descriptor vd, current_selection_border)
      BOOST_FOREACH(vertex_descriptor nvd, vertices_around_target( halfedge(vd,fg), fg))
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
