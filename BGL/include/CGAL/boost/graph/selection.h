// Copyright (c) 2015  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_BOOST_GRAPH_SELECTION_H
#define CGAL_BOOST_GRAPH_SELECTION_H

#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <boost/unordered_set.hpp>

#include <CGAL/boost/graph/Dual.h>
#include <boost/graph/filtered_graph.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <CGAL/boost/graph/alpha_expansion_graphcut.h>
#include <CGAL/squared_distance_3.h>

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

  for(face_descriptor fd : face_range)
  {
    for(halfedge_descriptor h :
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

template <typename GeomTraits,
          typename FaceGraph,
          typename IsSelectedMap,
          typename FaceIndexMap,
          typename VertexPointMap>
struct Regularization_graph
{
  typedef boost::graph_traits<FaceGraph> GT;
  typedef typename GT::face_descriptor fg_face_descriptor;
  typedef typename GT::face_iterator fg_face_iterator;
  typedef typename GT::halfedge_descriptor fg_halfedge_descriptor;
  typedef typename GT::edge_descriptor fg_edge_descriptor;
  typedef typename GT::edge_iterator fg_edge_iterator;
  typedef typename GT::vertex_descriptor fg_vertex_descriptor;

  typedef fg_face_descriptor vertex_descriptor;
  typedef fg_face_iterator vertex_iterator;
  typedef fg_edge_descriptor edge_descriptor;
  typedef boost::undirected_tag directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  typedef boost::edge_list_graph_tag traversal_category;

  struct Filter_border_edges
  {
    FaceGraph* fg;
    Filter_border_edges (FaceGraph& fg) : fg (&fg) { }
    bool operator() (const fg_edge_descriptor ed) const
    {
      return !is_border (ed, *fg);
    }
  };

  typedef boost::filter_iterator<Filter_border_edges, fg_edge_iterator> edge_iterator;

  struct Vertex_label_map
  {
    typedef vertex_descriptor key_type;
    typedef std::size_t value_type;
    typedef std::size_t& reference;
    typedef boost::lvalue_property_map_tag category;

    Regularization_graph* rg;

    Vertex_label_map (Regularization_graph* rg)
      : rg (rg) { }

    friend reference get (const Vertex_label_map& map, key_type k)
    {
      return (map.rg->labels)[get(map.rg->face_index_map,k)];
    }
    friend void put (const Vertex_label_map& map, key_type k, const value_type& v)
    {
      (map.rg->labels)[get(map.rg->face_index_map,k)] = v;
    }
  };

  struct Vertex_label_probability_map
  {
    typedef vertex_descriptor key_type;
    typedef std::vector<double> value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    const Regularization_graph* rg;

    Vertex_label_probability_map (const Regularization_graph* rg)
      : rg (rg)
    { }

    friend reference get (const Vertex_label_probability_map& pmap, key_type fd)
    {
      double value = (1. - pmap.rg->weight) * pmap.rg->area (fd) / pmap.rg->total_area;

      std::vector<double> out(2);
      if (get(pmap.rg->is_selected_map, fd))
      {
        if (pmap.rg->prevent_unselection)
          out[0] = (std::numeric_limits<double>::max)();
        else
          out[0] = value;
        out[1] = 0.;
      }
      else
      {
        out[0] = 0.;
        out[1] = value;
      }

      return out;
    }
  };

  struct Edge_cost_map
  {
    typedef edge_descriptor key_type;
    typedef double value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    const Regularization_graph* rg;

    Edge_cost_map (const Regularization_graph* rg)
      : rg (rg) { }

    friend reference get (const Edge_cost_map& pmap, key_type ed)
    {
      fg_vertex_descriptor esource = source(ed, pmap.rg->fg);
      fg_vertex_descriptor etarget = target(ed, pmap.rg->fg);

      // Cost
      double edge_length = std::sqrt(CGAL::squared_distance (get (pmap.rg->vertex_point_map, esource),
                                                             get (pmap.rg->vertex_point_map, etarget)));
      return pmap.rg->weight * edge_length / pmap.rg->total_length;
    }
  };

  FaceGraph& fg;
  IsSelectedMap is_selected_map;
  FaceIndexMap face_index_map;
  VertexPointMap vertex_point_map;
  double total_length;
  double total_area;
  double weight;
  bool prevent_unselection;
  std::vector<std::size_t> labels;

  Regularization_graph (FaceGraph& fg,
                        IsSelectedMap is_selected_map,
                        FaceIndexMap face_index_map,
                        VertexPointMap vertex_point_map,
                        double weight,
                        bool prevent_unselection)
    : fg (fg),
      is_selected_map (is_selected_map),
      face_index_map (face_index_map),
      vertex_point_map (vertex_point_map),
      total_length(0),
      total_area(0),
      weight (weight),
      prevent_unselection (prevent_unselection)
  {
    labels.reserve(num_faces(fg));
    std::size_t nb_selected = 0;
    for (fg_face_descriptor fd : faces(fg))
    {
      if (get(is_selected_map,fd))
      {
        labels.push_back(1);
        ++ nb_selected;
      }
      else
        labels.push_back(0);
    }

    // Compute normalization factors
    for (fg_edge_descriptor ed : edges(fg))
      total_length += length (ed);
    for (fg_face_descriptor fd : faces(fg))
      total_area += area (fd);
  }

  double length (fg_edge_descriptor ed) const
  {
    fg_vertex_descriptor esource = source(ed, fg);
    fg_vertex_descriptor etarget = target(ed, fg);
    return approximate_sqrt (typename GeomTraits::Compute_squared_distance_3()
                             (get (vertex_point_map, esource),
                              get (vertex_point_map, etarget)));
  }

  double area (fg_face_descriptor fd) const
  {
    fg_halfedge_descriptor hd = halfedge (fd, fg);
    fg_halfedge_descriptor nhd = next (hd, fg);

    return approximate_sqrt (typename GeomTraits::Compute_squared_area_3()
                             (get (vertex_point_map, source (hd, fg)),
                              get (vertex_point_map, target (hd, fg)),
                              get (vertex_point_map, target (nhd, fg))));
  }

  friend CGAL::Iterator_range<vertex_iterator>
  vertices (const Regularization_graph& graph)
  {
    return faces (graph.fg);
  }

  friend std::size_t num_vertices (const Regularization_graph& graph) { return num_faces(graph.fg); }

  friend CGAL::Iterator_range<edge_iterator>
  edges (const Regularization_graph& graph)
  {
    return CGAL::make_range (boost::make_filter_iterator
                             (Filter_border_edges(graph.fg),
                              begin(edges(graph.fg)), end(edges(graph.fg))),
                             boost::make_filter_iterator
                             (Filter_border_edges(graph.fg),
                              end(edges(graph.fg)), end(edges(graph.fg))));
  }

  friend vertex_descriptor source (edge_descriptor ed, const Regularization_graph& graph)
  {
    return face (halfedge (ed, graph.fg), graph.fg);
  }

  friend vertex_descriptor target (edge_descriptor ed, const Regularization_graph& graph)
  {
    return face (opposite(halfedge (ed, graph.fg), graph.fg), graph.fg);
  }

  Vertex_label_map vertex_label_map() { return Vertex_label_map(this); }
  Vertex_label_probability_map vertex_label_probability_map() const
  { return Vertex_label_probability_map(this); }
  Edge_cost_map edge_cost_map() const
  { return Edge_cost_map(this); }
};


} //end of namespace internal


/*!
\ingroup PkgBGLSelectionFct
augments a selection with faces of `fg` that are adjacent
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
       to accommodate new faces added to the selection.
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
    for(halfedge_descriptor hd : selection_boundary_halfedges)
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
    for(face_descriptor fd : new_selection_set)
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
diminishes a selection of faces from faces adjacent to a non-selected face.
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
       to accommodate faces removed from the selection.
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
    for(halfedge_descriptor hd : selection_boundary_halfedges)
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
    for(face_descriptor fd : elements_to_remove)
    {
      *out++=fd;
      put( is_selected, fd, false );
    }

    // update the set of currently selected faces
    std::vector<face_descriptor> new_selection;
    for(face_descriptor fd : current_selection)
      if ( !elements_to_remove.count(fd) )
        new_selection.push_back(fd);
    current_selection.swap(new_selection);
  }
  return out;
}

/*!
  \ingroup PkgBGLSelectionFct

  regularizes a selection in order to minimize the length of the
  border of the selection.

  The alpha expansion algorithm is used (see
  `CGAL::alpha_expansion_graphcut()`) using the length of the edge
  between two faces as the edge cost and the initial
  selected/unselected property of a face as the face cost.

  If `prevent_unselection` is set to `true`, the cost of unselecting a
  face is set to infinity, which forces the regularization to only
  select new faces and ensures that the regularization keeps all
  selected faces.

  \tparam TriangleMesh a model of `FaceGraph`

  \tparam IsSelectedMap a model of `ReadWritePropertyMap` with
  `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type and
  `bool` as value type

  \tparam NamedParameters a sequence of named parameters

  \param mesh the mesh containing the selected faces.

  \param is_selected indicates if a face is part of the selection. It
  is updated by the function to accommodate faces added or removed
  from the selection.

  \param weight sets the tradeoff between data fidelity and
  regularity, ranging from 0 (no regularization at all, selection is
  left unaltered) to 1 (maximum regularization, usually selects or
  unselects everything so that the length of the border of the
  selection is 0)

  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `tm`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `TriangleMesh`.}
    \cgalParamNEnd

    \cgalParamNBegin{face_index_map}
      \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
                     as key type and `std::size_t` as value type}
      \cgalParamDefault{an automatically indexed internal map}
    \cgalParamNEnd

    \cgalParamNBegin{prevent_unselection}
      \cgalParamDescription{Boolean used to indicate if selection can be only extended or if it can also be shrinked.}
      \cgalParamType{`bool`}
      \cgalParamDefault{`false`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd

    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd
  \cgalNamedParamsEnd
*/
template <typename TriangleMesh, typename IsSelectedMap, typename NamedParameters>
void
regularize_face_selection_borders(
  TriangleMesh& mesh,
  IsSelectedMap is_selected,
  double weight,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition (0.0 <= weight && weight < 1.0);

  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor mesh_face_descriptor;

  typedef typename GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::type FaceIndexMap;
  FaceIndexMap face_index_map = CGAL::get_initialized_face_index_map(mesh, np);

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vertex_point_map
    = choose_parameter(get_parameter(np, internal_np::vertex_point),
                       get_const_property_map(vertex_point, mesh));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Kernel;

  bool prevent_unselection = choose_parameter(get_parameter(np, internal_np::prevent_unselection),
                                              false);

  internal::Regularization_graph<Kernel, TriangleMesh, IsSelectedMap, FaceIndexMap,
                                 VertexPointMap>
    graph (mesh, is_selected,
           face_index_map,
           vertex_point_map,
           weight,
           prevent_unselection);

  alpha_expansion_graphcut (graph,
                            graph.edge_cost_map(),
                            graph.vertex_label_probability_map(),
                            graph.vertex_label_map(),
                            CGAL::parameters::vertex_index_map
                            (face_index_map));

  for (mesh_face_descriptor fd : faces(mesh))
    put(is_selected, fd, graph.labels[get(face_index_map,fd)]);
}

/// \cond SKIP_IN_MANUAL
// variant with default np
template <typename TriangleMesh, typename IsSelectedMap>
void
regularize_face_selection_borders(
  TriangleMesh& fg,
  IsSelectedMap is_selected,
  double weight)
{
  regularize_face_selection_borders (fg, is_selected, weight,
                                     CGAL::parameters::all_default());
}
/// \endcond

/// \cond SKIP_IN_MANUAL
// TODO: improve and document if useful
//
// Variant of regularization without graphcut but with brut-force
// local expansions. Can be interesting in some cases but too
// experimental/messy so far to be officially integrated.
template <class FaceGraph, class IsSelectedMap, class VertexPointMap>
void
regularize_face_selection_borders(
  FaceGraph& fg,
  IsSelectedMap is_selected,
  VertexPointMap vertex_point_map)
{
  typedef boost::graph_traits<FaceGraph> GT;
  typedef typename GT::face_descriptor fg_face_descriptor;
  typedef typename GT::halfedge_descriptor fg_halfedge_descriptor;
  typedef typename GT::edge_descriptor fg_edge_descriptor;
  typedef typename GT::vertex_descriptor fg_vertex_descriptor;

  // TODO: this is a quick and dirty version, the complexity is
  // crazy and it should be easy to do better (with priority queues,
  // for example)

  auto border_length =
    [&]() -> double
    {
      double out = 0.;
      for(fg_edge_descriptor ed : edges(fg))
      {
        fg_face_descriptor f0 = face (halfedge (ed, fg), fg);
        fg_face_descriptor f1 = face (opposite(halfedge (ed, fg), fg), fg);
        if (get(is_selected,f0) == get(is_selected,f1))
          continue;

        fg_vertex_descriptor esource = source(ed, fg);
        fg_vertex_descriptor etarget = target(ed, fg);

        out += std::sqrt(CGAL::squared_distance (get (vertex_point_map, esource),
                                                 get (vertex_point_map, etarget)));
      }
      return out;
    };

  // First: try edges
  while (true)
  {
    fg_edge_descriptor chosen;
    double length_before = border_length();
    double shortest_length = length_before;

    for (fg_edge_descriptor ed : edges(fg))
    {
      fg_face_descriptor selected = face (halfedge (ed, fg), fg);
      fg_face_descriptor unselected = face (opposite(halfedge (ed, fg), fg), fg);
      if (get(is_selected,selected) == get(is_selected,unselected))
        continue;

      if (get(is_selected, unselected))
        std::swap (selected, unselected);

      put(is_selected, unselected, true);
      double length_after = border_length();

      if (length_after < shortest_length)
      {
        chosen = ed;
        shortest_length = length_after;
      }

      // Cancel
      put(is_selected, unselected, false);
    }

    if (shortest_length == length_before)
      break;

    fg_face_descriptor selected = face (halfedge (chosen, fg), fg);
    fg_face_descriptor unselected = face (opposite(halfedge (chosen, fg), fg), fg);
    if (get(is_selected,selected) == get(is_selected,unselected))
      continue;

    if (get(is_selected, unselected))
      std::swap (selected, unselected);

    put(is_selected, unselected, true);
  }

  // Second: try 1-ring of vertices
  while (true)
  {
    fg_vertex_descriptor chosen;
    double length_before = border_length();
    double shortest_length = length_before;

    for (fg_vertex_descriptor vd : vertices(fg))
    {
      fg_halfedge_descriptor hd = halfedge(vd, fg);
      bool adjacent_to_selected = false, adjacent_to_nonselected = false;
      for (fg_face_descriptor fd : faces_around_target (hd, fg))
      {
        if (get(is_selected, fd))
          adjacent_to_selected = true;
        else
          adjacent_to_nonselected = true;

        if (adjacent_to_selected && adjacent_to_nonselected)
          break;
      }

      if (!(adjacent_to_selected && adjacent_to_nonselected))
        continue;

      std::vector<fg_face_descriptor> newly_selected;
      for (fg_face_descriptor fd : faces_around_target (hd, fg))
      {
        if (!get(is_selected, fd))
        {
          newly_selected.push_back (fd);
          put(is_selected, fd, true);
        }
      }
      double length_after = border_length();

      if (length_after < shortest_length)
      {
        chosen = vd;
        shortest_length = length_after;
      }

      // Cancel
      for (fg_face_descriptor fd : newly_selected)
        put(is_selected, fd, false);
    }

    if (shortest_length == length_before)
      break;

    fg_halfedge_descriptor hd = halfedge (chosen, fg);

    for (fg_face_descriptor fd : faces_around_target (hd, fg))
      put(is_selected, fd, true);
  }
}
/// \endcond


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
  for(halfedge_descriptor hd : hedges)
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
augments a selection with edges of `fg` that are adjacent
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
       to accommodate new edges added to the selection.
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
    for(edge_descriptor ed : current_selection)
    {
      halfedge_descriptor hdi=halfedge(ed,fg);
      for(halfedge_descriptor hd : halfedges_around_source( hdi, fg))
      {
        edge_descriptor ned=edge(hd, fg);
        if (!get(is_selected, ned)) new_selection_set.insert(ned);
      }
      for(halfedge_descriptor hd : halfedges_around_target( hdi, fg))
      {
        edge_descriptor ned=edge(hd, fg);
        if (!get(is_selected, ned)) new_selection_set.insert(ned);
      }
    }

    // extract unique selection
    std::vector<edge_descriptor> new_selection;
    for(edge_descriptor ed : new_selection_set)
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
diminishes a selection of edges from edges adjacent to a non-selected edge.
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
       to accommodate edges removed from the selection.
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
  for(edge_descriptor ed : selection)
  {
    halfedge_descriptor hd=halfedge(ed,fg);
    for(halfedge_descriptor nhd : halfedges_around_source( hd, fg))
    {
      edge_descriptor ned=edge(nhd, fg);
      if (!get(is_selected, ned)){
        unique_vertex_set.insert(source(hd,fg));
        break;
      }
    }
    for(halfedge_descriptor nhd : halfedges_around_target( hd, fg))
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
    for(vertex_descriptor vd : current_selection_border)
      for(halfedge_descriptor hd : halfedges_around_target( halfedge(vd,fg), fg))
      {
        edge_descriptor ed = edge(hd, fg);
        if (get(is_selected, ed)){
          edges_to_deselect.insert(ed);
          unique_vertex_set.insert(source(hd, fg));
        }
      }

    // extract unique selection
    for(edge_descriptor ed : edges_to_deselect)
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
augments a selection with vertices of `fg` that are adjacent
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
       to accommodate new vertices added to the selection.
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
    for(vertex_descriptor vd : current_selection)
      for(vertex_descriptor nvd : vertices_around_target( halfedge(vd,fg), fg))
        if (!get(is_selected, nvd)) new_selection_set.insert(nvd);

    // extract unique selection
    std::vector<vertex_descriptor> new_selection;
    for(vertex_descriptor vd : new_selection_set)
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
diminishes a selection of vertices from vertices adjacent to a non-selected vertex.
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
       to accommodate vertices removed from the selection.
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
  for(vertex_descriptor vd : selection)
    for(vertex_descriptor nvd : vertices_around_target( halfedge(vd,fg), fg))
        if (!get(is_selected, nvd)) unique_vertex_set.insert(nvd);

  std::vector<vertex_descriptor> current_selection_border(unique_vertex_set.begin(), unique_vertex_set.end());
  for (unsigned int i=0; i<k; ++i)
  {
    if (current_selection_border.empty()) break;

    //collect adjacent vertices selected
    std::set<vertex_descriptor> vertices_to_deselect;
    for(vertex_descriptor vd : current_selection_border)
      for(vertex_descriptor nvd : vertices_around_target( halfedge(vd,fg), fg))
        if (get(is_selected, nvd)) vertices_to_deselect.insert(nvd);

    // extract unique selection
    std::vector<vertex_descriptor> new_selection_border;
    for(vertex_descriptor vd : vertices_to_deselect)
    {
      *out++=vd;
      new_selection_border.push_back(vd);
      put( is_selected, vd, false );
    }
    current_selection_border.swap(new_selection_border);
  }
  return out;
}

/**
 * \ingroup PkgBGLSelectionFct
 *
 * Expands a selection of faces so that their removal does not create any non manifold vertex.
 * For each vertex that is incident to a selected face, we turn around that vertex and check
 * if there is exactly one set of consecutive selected faces. If not, additional faces around
 * that vertex are selected to match this condition.
 *
 * @tparam TriangleMesh a model of `FaceGraph` that is triangulated.
 * @tparam FaceRange a range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 * with an iterator type model of `ForwardIterator`.
 * @tparam  IsSelectedMap a model of `ReadWritePropertyMap` with
 * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `bool` as value.

 * @param tm the triangle mesh.
 * @param faces_to_be_deleted the range of selected faces.
 * @param is_selected a property map containing the selected-or-not status of each face of `tm`.
 * It must associate `true` to each face of `faces_to_be_deleted` and `false` to any other face of `tm`.
 * It will be modified if the face selection must be expanded.
 *
 **/
template<class TriangleMesh, class FaceRange, class IsSelectedMap>
void expand_face_selection_for_removal(const FaceRange& faces_to_be_deleted,
                                       TriangleMesh& tm,
                                       IsSelectedMap is_selected)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  boost::unordered_set<vertex_descriptor> vertices_queue;

  // collect vertices belonging to at least a triangle that will be removed
  for(face_descriptor fd : faces_to_be_deleted)
  {
    halfedge_descriptor h = halfedge(fd, tm);
    vertices_queue.insert( target(h, tm) );
    vertices_queue.insert( target(next(h, tm), tm) );
    vertices_queue.insert( target(prev(h, tm), tm) );
  }

  while (!vertices_queue.empty())
  {
    vertex_descriptor vd = *vertices_queue.begin();
    vertices_queue.erase( vertices_queue.begin() );

    // make sure hd is not a border halfedge
    halfedge_descriptor hd = halfedge(vd, tm);
    while(is_border(hd,tm) || ( !get(is_selected, face(hd, tm))) )
    {
      hd = opposite( next(hd, tm), tm);
      CGAL_assertion( hd != halfedge(vd, tm) );
    }

    // set hd to the last selected face of a connected component
    // of selected faces around a vertex
    halfedge_descriptor start = hd;
    halfedge_descriptor next_around_vertex = opposite( next(hd, tm), tm);
    while(is_border(next_around_vertex,tm) || get(is_selected, face(next_around_vertex, tm) ) )
    {
      hd = next_around_vertex;
      next_around_vertex = opposite( next(hd, tm), tm);
      if (hd==start) break;
    }
    if ( is_border(next_around_vertex,tm) || get(is_selected, face(next_around_vertex, tm) ) ) continue; //all incident faces will be removed

    while( true )
    {
      // collect non-selected faces
      std::vector<halfedge_descriptor> faces_traversed;
      bool non_selected_face_range_has_boundary = false; // handle non-manifold situations when crossing a border
      do
      {
        faces_traversed.push_back(next_around_vertex);
        next_around_vertex = opposite( next(next_around_vertex, tm), tm);
        if (is_border(next_around_vertex,tm))
        {
          next_around_vertex = opposite( next(next_around_vertex, tm), tm);
          if (!get(is_selected, face(next_around_vertex, tm) ))
          {
            non_selected_face_range_has_boundary=true; // always non-manifold after removal of the selection
            break;
          }
        }
        CGAL_assertion(!is_border(next_around_vertex,tm));
      }
      while( !get(is_selected, face(next_around_vertex, tm) ) );

      if (!non_selected_face_range_has_boundary)
      {
        // go over the connected components of faces to remove
        do{
          if (next_around_vertex==start)
            break;
          next_around_vertex = opposite( next(next_around_vertex, tm), tm);
        }
        while(is_border(next_around_vertex,tm) || get(is_selected, face(next_around_vertex, tm) ) );

        if (next_around_vertex==start)
          break;
      }
      // else we simply mark the range of traversed faces and start a new range after the border

      for(halfedge_descriptor f_hd : faces_traversed)
      {
        assert(target(f_hd, tm) == vd);
        put(is_selected, face(f_hd, tm), true);
        vertices_queue.insert( target( next(f_hd, tm), tm) );
        vertices_queue.insert( source(f_hd, tm) );
      }
    }
  }
}

//todo: take non-manifold vertices into account.
template<class PolygonMesh, class FaceRange>
bool is_selection_a_topological_disk(const FaceRange& face_selection,
                                           PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  boost::unordered_set<vertex_descriptor> sel_vertices;
  boost::unordered_set<edge_descriptor> sel_edges;
  for(face_descriptor f : face_selection)
  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(f, pm), pm))
    {
      sel_vertices.insert(target(h, pm));
      sel_edges.insert(edge(h,pm));
    }
  }
  return (sel_vertices.size() - sel_edges.size() + face_selection.size() == 1);
}
} //end of namespace CGAL

#endif //CGAL_BOOST_GRAPH_SELECTION_H
