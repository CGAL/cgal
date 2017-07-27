// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Maxime Gimeno

#ifndef CGAL_BOOST_GRAPH_FACE_FILTERED_GRAPH_H
#define CGAL_BOOST_GRAPH_FACE_FILTERED_GRAPH_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/assertions.h>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/has_range_iterator.hpp>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL
{

  /*!
   * \ingroup PkgBGLAdaptors
   *
   * The class `Face_filtered_graph` is an adaptor that creates a filtered view of a graph
   * by restricting it to a subset of faces. Contrary to
   * <a href="http://www.boost.org/doc/libs/release/libs/graph/doc/filtered_graph.html"><code>boost::filtered_graph</code></a>,
   * this class only requires a way to access the selected faces and will automatically select the
   * edges/halfedges and vertices present in the adapted graph. A vertex is selected if it is incident to at least one
   * selected face. A edge is selected if it is incident to at least a selected face. A halfedge is selected if its edge
   * is selected.
   *
   * Since this class is a model of the `FaceGraph` concept, there is a restriction on the set of selected faces:
   * the adapted graph must define a manifold mesh. In order to check that this condition is verified, you can
   * use the function `is_selection_valid()`.
   *
   * There are two different ways to initialize this class. You can directly provide the set of faces selected, or
   * if you have a face patch map, select the patches of faces. The latter option is convenient if you want to access
   * some connected components of a graph after having called `CGAL::Polygon_mesh_processing::connected_components()`.
   *
   * The documented interface of this class is limited on purpose and free functions of the concept
   * this class is a model of must be used to manipulate it.
   *
   * A BGL-like named parameter mechanism is used in the constructors of this class. %Default values are available but if you need
   * to set them, you can pass for `np` `CGAL::parameters::face_index_map(fim).halfedge_index_map(him).vertex_index_map(vim)`
   * where `fim`, `him`, and `vim` are the respective index maps. The order of the arguments is not important and any of them can be
   * missing if the default is fine.
   *
   * \tparam Graph must be a model of a `FaceListGraph`, `HalfedgeListGraph`, and \bgllink{VertexListGraph}.
   * \tparam FIMap a model of `ReadablePropertyMap` with `face_descriptor` as key and `graph_traits<Graph>::%faces_size_type` as value
   * \tparam VIMap a model of `ReadablePropertyMap` with `vertex_descriptor` as key and `graph_traits<Graph>::%vertices_size_type` as value
   * \tparam HIMap a model of `ReadablePropertyMap` with `halfedge_descriptor` as key and `graph_traits<Graph>::%halfedges_size_type` as value
   *
   * \cgalModels `FaceListGraph`
   * \cgalModels `HalfedgeListGraph`
   * \cgalModels \bgllink{VertexListGraph}
   */
template<typename Graph,
         typename FIMap = typename boost::property_map<Graph, CGAL::face_index_t>::type,
         typename VIMap = typename boost::property_map<Graph, boost::vertex_index_t>::type,
         typename HIMap = typename boost::property_map<Graph, CGAL::halfedge_index_t>::type>
struct Face_filtered_graph
{
  typedef boost::graph_traits<Graph>                  gt;
  /// Vertex descriptor type
  typedef typename boost::graph_traits<Graph>::vertex_descriptor              vertex_descriptor;
  /// Halfedge descriptor type
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor            halfedge_descriptor;
  /// Edge descriptor type
  typedef typename boost::graph_traits<Graph>::edge_descriptor                edge_descriptor;
  /// Face descriptor type
  typedef typename boost::graph_traits<Graph>::face_descriptor                face_descriptor;
  /// Size type
  #ifndef DOXYGEN_RUNNING
  typedef boost::dynamic_bitset<>::size_type size_type;
  #else
  typedef unspecified_type size_type;
  #endif

  // non documented types
  typedef typename boost::property_traits< FIMap >::value_type face_index_type;
  typedef typename boost::property_traits< VIMap >::value_type vertex_index_type;
  typedef typename boost::property_traits< HIMap >::value_type halfedge_index_type;
  typedef Face_filtered_graph<Graph, FIMap, VIMap, HIMap>   Self;

  /*!
   * \brief Constructor where the set of selected faces is specified as a range of patch ids.
   *
   * \tparam FacePatchIndexMap a model of `ReadablePropertyMap` with
      `face_descriptor` as key type and
      `graph_traits<Graph>::%faces_size_type` as value type.
   * \tparam FacePatchIndexRange a model of `ConstRange` with `boost::property_traits<FacePatchIndexMap>::%value_type` as value type.
   * \tparam NamedParameters a sequence of named parameters
   *
   * \param graph the underlying graph.
   * \param face_patch_index_map the property_map that assigns a patch index to each face, with
      `face_descriptor` as key type and `boost::graph_traits<Graph>::%faces_size_type` as value type.
   * \param selected_face_patch_indices a range of the face patch indices to select.
   * \param np optional sequence of named parameters among the ones listed below
   *
   * \cgalNamedParamsBegin
   *   \cgalParamBegin{face_index_map}
   *     a property map containing an index for each face initialized from 0 to `num_vertices(graph)`
   *   \cgalParamEnd
   *   \cgalParamBegin{vertex_index_map}
   *     a property map containing an index for each vertex initialized 0 to `num_vertices(vertex)`
   *   \cgalParamEnd
   *   \cgalParamBegin{halfedge_index_map}
   *     a property map containing an index for each halfedge initialized 0 to `num_halfedges(graph)`
   *   \cgalParamEnd
   * \cgalNamedParamsEnd
   */
  template <typename FacePatchIndexMap, class FacePatchIndexRange, class CGAL_BGL_NP_TEMPLATE_PARAMETERS>
  Face_filtered_graph(const Graph& graph,
                      const FacePatchIndexRange& selected_face_patch_indices,
                            FacePatchIndexMap face_patch_index_map,
                            const CGAL_BGL_NP_CLASS& np
                           #ifndef DOXYGEN_RUNNING
                              , typename boost::enable_if<
                                  typename boost::has_range_const_iterator<FacePatchIndexRange>::type
                                >::type* = 0
                           #endif
                             )
    : _graph(const_cast<Graph&>(graph))
    , fimap(boost::choose_param(get_param(np, internal_np::face_index), get_const_property_map(face_index, graph)))
    , vimap(boost::choose_param(get_param(np, boost::vertex_index), get_const_property_map(boost::vertex_index, graph)))
    , himap(boost::choose_param(get_param(np, internal_np::halfedge_index), get_const_property_map(halfedge_index, graph)))
  {
    set_selected_faces(selected_face_patch_indices, face_patch_index_map);
  }

  template <typename FacePatchIndexMap, class FacePatchIndexRange>
  Face_filtered_graph(const Graph& graph,
                      const FacePatchIndexRange& selected_face_patch_indices,
                      FacePatchIndexMap face_patch_index_map
                      , typename boost::enable_if<
                      typename boost::has_range_const_iterator<FacePatchIndexRange>::type
                      >::type* = 0
                      )
    : _graph(const_cast<Graph&>(graph))
    , fimap(get(CGAL::face_index, graph))
    , vimap(get(boost::vertex_index, graph))
    , himap(get(CGAL::halfedge_index, graph))
  {
    set_selected_faces(selected_face_patch_indices, face_patch_index_map);
  }
  /*!
   * \brief Constructor where the set of selected faces is specified as a patch id.
   *
   * \tparam FacePatchIndexMap a model of `ReadablePropertyMap` with
      `face_descriptor` as key type and
      `graph_traits<Graph>::%faces_size_type` as value type.
   * \tparam NamedParameters a sequence of named parameters
   * \param graph the underlying graph.
   * \param face_patch_index_map the property_map that assigns a patch index to each face, with
      `face_descriptor` as key type and
      `graph_traits<Graph>::%faces_size_type` as value type.
   * \param selected_face_patch_index the index of the face patch selected.
   * \param np optional sequence of named parameters among the ones listed below
   *
   * \cgalNamedParamsBegin
   *   \cgalParamBegin{face_index_map}
   *     a property map containing an index for each face initialized from 0 to `num_vertices(graph)`
   *   \cgalParamEnd
   *   \cgalParamBegin{vertex_index_map}
   *     a property map containing an index for each vertex initialized 0 to `num_vertices(vertex)`
   *   \cgalParamEnd
   *   \cgalParamBegin{halfedge_index_map}
   *     a property map containing an index for each halfedge initialized 0 to `num_halfedges(graph)`
   *   \cgalParamEnd
   * \cgalNamedParamsEnd
   */
  template <typename FacePatchIndexMap, class CGAL_BGL_NP_TEMPLATE_PARAMETERS>
  Face_filtered_graph(const Graph& graph,
                            typename boost::property_traits<FacePatchIndexMap>::value_type selected_face_patch_index,
                            FacePatchIndexMap face_patch_index_map,
                            const CGAL_BGL_NP_CLASS& np
                             )
    : _graph(const_cast<Graph&>(graph))
    , fimap(boost::choose_param(get_param(np, internal_np::face_index), get_const_property_map(face_index, graph)))
    , vimap(boost::choose_param(get_param(np, boost::vertex_index), get_const_property_map(boost::vertex_index, graph)))
    , himap(boost::choose_param(get_param(np, internal_np::halfedge_index), get_const_property_map(halfedge_index, graph)))
  {
    set_selected_faces(selected_face_patch_index, face_patch_index_map);
  }

  template <typename FacePatchIndexMap>
  Face_filtered_graph(const Graph& graph,
                            typename boost::property_traits<FacePatchIndexMap>::value_type pid,
                            FacePatchIndexMap face_patch_index_map)
    : _graph(const_cast<Graph&>(graph))
    , fimap(get(CGAL::face_index, graph))
    , vimap(get(boost::vertex_index, graph))
    , himap(get(CGAL::halfedge_index, graph))
  {
    set_selected_faces(pid, face_patch_index_map);
  }

  /*!
   * \brief Constructor where the set of selected faces is specified as a range of face descriptors.
   *
   * \tparam FaceRange a model of `ConstRange` with `face_descriptor` as value type.
   * \tparam NamedParameters a sequence of named parameters
   * \param graph the graph containing the wanted patch.
   * \param selected_faces the set of selected faces.
   * \param np optional sequence of named parameters among the ones listed below
   *
   * \cgalNamedParamsBegin
   *   \cgalParamBegin{face_index_map}
   *     a property map containing an index for each face initialized from 0 to `num_vertices(graph)`
   *   \cgalParamEnd
   *   \cgalParamBegin{vertex_index_map}
   *     a property map containing an index for each vertex initialized 0 to `num_vertices(vertex)`
   *   \cgalParamEnd
   *   \cgalParamBegin{halfedge_index_map}
   *     a property map containing an index for each halfedge initialized 0 to `num_halfedges(graph)`
   *   \cgalParamEnd
   * \cgalNamedParamsEnd
   */
  template <typename FaceRange, class CGAL_BGL_NP_TEMPLATE_PARAMETERS>
  Face_filtered_graph(const Graph& graph,
                      const FaceRange& selected_faces,
                      const CGAL_BGL_NP_CLASS& np)
    : _graph(const_cast<Graph&>(graph))
    , fimap(boost::choose_param(get_param(np, internal_np::face_index), get_const_property_map(face_index, graph)))
    , vimap(boost::choose_param(get_param(np, boost::vertex_index), get_const_property_map(boost::vertex_index, graph)))
    , himap(boost::choose_param(get_param(np, internal_np::halfedge_index), get_const_property_map(halfedge_index, graph)))
  {
    set_selected_faces(selected_faces);
  }

  template <typename FaceRange>
  Face_filtered_graph(const Graph& graph,
                      const FaceRange& selected_faces)
    : _graph(const_cast<Graph&>(graph))
    , fimap(get(CGAL::face_index, graph))
    , vimap(get(boost::vertex_index, graph))
    , himap(get(CGAL::halfedge_index, graph))
  {
    set_selected_faces(selected_faces);
  }

  ///returns a const reference to the underlying graph.
  const Graph& graph()const{ return _graph; }
  ///returns a reference to the underlying graph.
  Graph& graph(){ return _graph; }

  ///change the set of selected faces using a patch id
  template<class FacePatchIndexMap>
  void set_selected_faces(typename boost::property_traits<FacePatchIndexMap>::value_type face_patch_id,
                          FacePatchIndexMap face_patch_index_map)
  {
    face_indices.clear();
    vertex_indices.clear();
    halfedge_indices.clear();

    selected_faces.resize(num_faces(_graph));
    selected_vertices.resize(num_vertices(_graph));
    selected_halfedges.resize(num_halfedges(_graph));
    selected_faces.reset();
    selected_vertices.reset();
    selected_halfedges.reset();
    BOOST_FOREACH(face_descriptor fd, faces(_graph) )
    {
      if(get(face_patch_index_map, fd) == face_patch_id)
      {
        selected_faces.set(get(fimap, fd));
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, _graph), _graph))
        {
          selected_halfedges.set(get(himap, hd));
          selected_halfedges.set(get(himap, opposite(hd, _graph)));
          selected_vertices.set(get(vimap, target(hd, _graph)));
        }
      }
    }
  }
  /// change the set of selected faces using a range of patch ids
  template<class FacePatchIndexRange, class FacePatchIndexMap>
  void set_selected_faces(const FacePatchIndexRange& selected_face_patch_indices,
                          FacePatchIndexMap face_patch_index_map
  #ifndef DOXYGEN_RUNNING
                          , typename boost::enable_if<
                              typename boost::has_range_const_iterator<FacePatchIndexRange>::type
                            >::type* = 0
  #endif
  )
  {
    face_indices.clear();
    vertex_indices.clear();
    halfedge_indices.clear();

    selected_faces.resize(num_faces(_graph));
    selected_vertices.resize(num_vertices(_graph));
    selected_halfedges.resize(num_halfedges(_graph));
    selected_faces.reset();
    selected_vertices.reset();
    selected_halfedges.reset();
    typedef typename boost::property_traits<FacePatchIndexMap>::value_type Patch_index;
    boost::unordered_set<Patch_index> pids(boost::begin(selected_face_patch_indices),
                                           boost::end(selected_face_patch_indices));

    BOOST_FOREACH(face_descriptor fd, faces(_graph) )
    {
      if(pids.count(get(face_patch_index_map, fd)) != 0)
      {
        selected_faces.set(get(fimap, fd));
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, _graph), _graph))
        {
          selected_halfedges.set(get(himap, hd));
          selected_halfedges.set(get(himap, opposite(hd, _graph)));
          selected_vertices.set(get(vimap, target(hd, _graph)));
        }
      }
    }
  }
  /// change the set of selected faces using a range of face descriptors
  template<class FaceRange>
  void set_selected_faces(const FaceRange& selection)
  {
    face_indices.clear();
    vertex_indices.clear();
    halfedge_indices.clear();

    selected_faces.resize(num_faces(_graph));
    selected_vertices.resize(num_vertices(_graph));
    selected_halfedges.resize(num_halfedges(_graph));
    selected_faces.reset();
    selected_vertices.reset();
    selected_halfedges.reset();
    BOOST_FOREACH(face_descriptor fd, selection)
    {
      selected_faces.set(get(fimap, fd));
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, _graph), _graph))
      {
        selected_halfedges.set(get(himap, hd));
        selected_halfedges.set(get(himap, opposite(hd, _graph)));
        selected_vertices.set(get(vimap, target(hd, _graph)));
      }
    }
  }

  struct Is_simplex_valid
  {
    Is_simplex_valid(const Self* graph)
      :adapter(graph)
    {}

    Is_simplex_valid()
      :adapter(NULL)
    {}
    template<typename Simplex>
    bool operator()(Simplex s)
    {
      CGAL_assertion(adapter!=NULL);
      return (adapter->is_in_cc(s));
    }
    const Self* adapter;
  };

  bool is_in_cc(face_descriptor f) const
  {
    return selected_faces[get(fimap, f)];
  }

  bool is_in_cc(vertex_descriptor v) const
  {
    return selected_vertices[get(vimap, v)];
  }

  bool is_in_cc(halfedge_descriptor h) const
  {
    return selected_halfedges[get(himap, h)];
  }

  bool is_in_cc(edge_descriptor e) const
  {
    return selected_halfedges[get(himap, halfedge(e,_graph))];
  }
  ///returns the number of selected faces
  size_type number_of_faces()const
  {
    return selected_faces.count();
  }
///returns the number of selected vertices.
  size_type number_of_vertices()const
  {
    return selected_vertices.count();
  }
///returns the number of selected halfedges.
  size_type number_of_halfedges()const
  {
    return selected_halfedges.count();
  }

  Property_map_binder< FIMap, typename Pointer_property_map< typename boost::property_traits< FIMap >::value_type >::type >
  get_face_index_map() const
  {
    if (face_indices.empty())
    {
      face_index_type index = 0;
      face_indices.resize(num_faces(_graph));
      for (std::size_t i=selected_faces.find_first(); i < selected_faces.npos; i = selected_faces.find_next(i))
     {
        face_indices[i] = index++;
     }
    }
    return bind_property_maps(fimap, make_property_map(face_indices) );
  }

  Property_map_binder< VIMap, typename Pointer_property_map< typename boost::property_traits< VIMap >::value_type >::type >
  get_vertex_index_map() const
  {
    if (vertex_indices.empty())
    {
      vertex_index_type index = 0;
      vertex_indices.resize(num_vertices(_graph));
      for (std::size_t i=selected_vertices.find_first(); i < selected_vertices.npos; i = selected_vertices.find_next(i))
     {
        vertex_indices[i] = index++;
     }
    }
    return bind_property_maps(vimap, make_property_map(vertex_indices) );
  }

  Property_map_binder< HIMap, typename Pointer_property_map< typename boost::property_traits< HIMap >::value_type >::type >
  get_halfedge_index_map() const
  {
    if (halfedge_indices.empty())
    {
      halfedge_index_type index = 0;
      halfedge_indices.resize(num_halfedges(_graph));
      for (std::size_t i=selected_halfedges.find_first(); i < selected_halfedges.npos; i = selected_halfedges.find_next(i))
     {
        halfedge_indices[i] = index++;
     }
    }
    return bind_property_maps(himap, make_property_map(halfedge_indices) );
  }

  /// returns `true` if around any vertex of a selected face,
  /// there is at most one connected set of selected faces.
  bool is_selection_valid()
  {
    BOOST_FOREACH(vertex_descriptor vd, vertices(*this) )
    {
      face_descriptor first_selected = boost::graph_traits<Graph>::null_face();
      bool first_unselected_found(false),
           second_unselected_found(false);

      //find an unselected face, then find the first selected face.
      //Find another unselected face, the next selected face must be the first;
      //else this is not valid.
      halfedge_descriptor hd = halfedge(vd, _graph);
      face_descriptor first_tested = boost::graph_traits<Graph>::null_face();
      while(1) //will break if valid, return false if not valid
      {
        face_descriptor fd = face(hd, _graph);

        if(first_tested == boost::graph_traits<Graph>::null_face())
          first_tested = fd;
        else if(fd == first_tested )
        {
          //if there is no unselected face, break
          if(selected_faces[get(fimap, fd)] && !first_unselected_found)
          break;
          //if there is no selected face, break
          else if(!selected_faces[get(fimap, fd)] &&
             first_selected == boost::graph_traits<Graph>::null_face())
          break;
        }

        if(fd != boost::graph_traits<Graph>::null_face())
        {
          if(selected_faces[get(fimap, fd)])
          {
            if(first_unselected_found &&
               first_selected == boost::graph_traits<Graph>::null_face())
            {
              first_selected = fd;
            }
            else if(second_unselected_found)
            {
              if(fd == first_selected)
                break;
              else
                return false;
            }
          }
          else
          {
            if(first_selected == boost::graph_traits<Graph>::null_face())
              first_unselected_found = true;
            else
              second_unselected_found = true;
          }
        }
        hd = next(opposite(hd, _graph), _graph);
      }
    }
    return true;
  }

private:
  Graph& _graph;
  FIMap fimap;
  VIMap vimap;
  HIMap himap;
  boost::dynamic_bitset<> selected_faces;
  boost::dynamic_bitset<> selected_vertices;
  boost::dynamic_bitset<> selected_halfedges;
  mutable std::vector<face_index_type> face_indices;
  mutable std::vector<vertex_index_type> vertex_indices;
  mutable std::vector<halfedge_index_type> halfedge_indices;
};

} // namespace CGAL

namespace boost
{

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct graph_traits< CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >
{
  typedef CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap> G;
  typedef boost::graph_traits<Graph> BGTG;
  typedef typename BGTG::vertex_descriptor vertex_descriptor;
  typedef typename BGTG::halfedge_descriptor halfedge_descriptor;
  typedef typename BGTG::edge_descriptor edge_descriptor;
  typedef typename BGTG::face_descriptor face_descriptor;

  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::vertex_iterator>    vertex_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::halfedge_iterator>  halfedge_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::edge_iterator>      edge_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::face_iterator>      face_iterator;

  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::out_edge_iterator>  out_edge_iterator;
  typedef boost::filter_iterator<typename G::Is_simplex_valid, typename BGTG::in_edge_iterator>   in_edge_iterator;

  typedef typename BGTG::directed_category directed_category;
  typedef typename BGTG::edge_parallel_category edge_parallel_category;
  typedef typename BGTG::traversal_category traversal_category;
  typedef typename boost::dynamic_bitset<>::size_type vertices_size_type;
  typedef typename boost::dynamic_bitset<>::size_type edges_size_type;
  typedef typename boost::dynamic_bitset<>::size_type halfedges_size_type;
  typedef typename boost::dynamic_bitset<>::size_type faces_size_type;
  typedef typename BGTG::degree_size_type degree_size_type;

  static vertex_descriptor null_vertex()
  {
    return BGTG::null_vertex();
  }

  static halfedge_descriptor null_halfedge()
  {
    return BGTG::null_halfedge();
  }

  static edge_descriptor null_edge()
  {
    return edge_descriptor(BGTG::null_halfedge());
  }

  static face_descriptor null_face()
  {
    return BGTG::null_face();
  }
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct graph_traits< const CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >
    : public graph_traits< CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >
{};


} // namespace boost


namespace CGAL {

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap>::size_type
num_vertices(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.number_of_vertices();
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap>::size_type
num_edges(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.number_of_halfedges()/2;
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::degree_size_type
degree(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
       const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  CGAL_assertion(w.is_in_cc(v));
  typename boost::graph_traits<Graph>::degree_size_type v_deg = 0;
  typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h = halfedge(v, w);
  typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor hcirc = h;
  do
  {
    if(w.is_in_cc(hcirc))
      ++v_deg;
    hcirc = opposite(next(hcirc, w.graph()), w.graph());
  }while(hcirc != h);
  return v_deg;
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::degree_size_type
out_degree(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
           const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  CGAL_assertion(w.is_in_cc(v));
  return static_cast<typename boost::graph_traits<Graph>::degree_size_type>(
    std::distance(out_edges(v, w).first ,out_edges(v, w).second) );
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Graph>::degree_size_type
in_degree(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
          const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  CGAL_assertion(w.is_in_cc(v));
  return static_cast<typename boost::graph_traits<Graph>::degree_size_type>(
    std::distance(in_edges(v, w).first ,in_edges(v, w).second) );
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
source(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e,
       const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(e));
  return source(e, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
target(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e,
       const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(e));
  return target(e, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
std::pair<typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor, bool>
edge(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor u,
     typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
     const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(u) && w.is_in_cc(v));
  typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e = edge(u, v, w.graph()).first;
  bool res = w.is_in_cc(e);
  return std::make_pair(e, res);
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_iterator>
vertices(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<Graph >::vertex_iterator g_vertex_iterator;

  typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_vertex_iterator b,e;
  boost::tie(b,e) = vertices(w.graph());
  return make_range(vertex_iterator(predicate, b, e),
                    vertex_iterator(predicate, e, e));
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_iterator>
edges(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_iterator edge_iterator;
  typedef typename boost::graph_traits<Graph >::edge_iterator g_edge_iterator;

  typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_edge_iterator b,e;
  boost::tie(b,e) = edges(w.graph());
  return make_range(edge_iterator(predicate, b, e),
                    edge_iterator(predicate, e, e));
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::out_edge_iterator>
out_edges(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
          const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{

  typedef typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::out_edge_iterator out_edge_iterator;
  typedef typename boost::graph_traits<Graph >::out_edge_iterator g_out_edge_iterator;

  typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_out_edge_iterator b,e;
  boost::tie(b,e) = out_edges(v, w.graph());
  return make_range(out_edge_iterator(predicate, b, e),
                    out_edge_iterator(predicate, e, e));
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
CGAL::Iterator_range<typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::in_edge_iterator>
in_edges(typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
         const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{

  typedef typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::in_edge_iterator in_edge_iterator;
  typedef typename boost::graph_traits<Graph >::in_edge_iterator g_in_edge_iterator;

  typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  g_in_edge_iterator b,e;
  boost::tie(b,e) = in_edges(v, w.graph());
  return make_range(in_edge_iterator(predicate, b, e),
                    in_edge_iterator(predicate, e, e));
}

//
// HalfedgeGraph
//
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor
edge(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(h));
  return edge(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits<  Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::edge_descriptor e,
         const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(e));
  return halfedge(e, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
         const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(v));
  typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h = halfedge(v, w.graph());
  typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor hcirc = h;
  do
  {
    if(w.is_in_cc(hcirc))
      return hcirc;
    hcirc = opposite(next(hcirc, w.graph()), w.graph());
  }while(hcirc != h);
  return boost::graph_traits< CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::null_halfedge();
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
std::pair<typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor, bool>
halfedge(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor u,
         typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor v,
         const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(u) && w.is_in_cc(v));
  typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h = halfedge(u, v, w.graph()).first;
  return std::make_pair(h, w.is_in_cc(h));
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
opposite(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
         const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(h) );
     return opposite(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
source(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
       const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(h) );
  return source(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::vertex_descriptor
target(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
       const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(h) );
  return target(h, w.graph());
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
next(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(h));
  if(w.is_in_cc(next(h, w.graph())))
    return next(h, w.graph());

  //h is on the border of the selection
  CGAL_assertion( !w.is_in_cc(face(h, w.graph())) );
  typedef typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h_d;
  h_d candidate = next(h, w.graph());
  CGAL_assertion(!w.is_in_cc(candidate));
  do{
    candidate = next(opposite(candidate, w.graph()), w.graph());
    CGAL_assertion(candidate!=opposite(h,w.graph()));
  }while(!w.is_in_cc(candidate));
  return candidate;
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
prev(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{

  CGAL_assertion(w.is_in_cc(h));
  if(w.is_in_cc(prev(h, w.graph())))
    return prev(h, w.graph());

  //h is on the border of the selection
  CGAL_assertion( !w.is_in_cc(face(h, w.graph())) );
  typedef typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h_d;
  h_d candidate = prev(h, w.graph());
  CGAL_assertion(!w.is_in_cc(candidate));
  do{
    candidate = prev(opposite(candidate, w.graph()), w.graph());
    CGAL_assertion(candidate!=opposite(h,w.graph()));
  }while(!w.is_in_cc(candidate));
  return candidate;
}

//
// HalfedgeListGraph
//

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
std::pair<typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_iterator,
typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_iterator>
halfedges(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits<Graph >::halfedge_iterator g_halfedge_iterator;

  typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  std::pair<g_halfedge_iterator, g_halfedge_iterator> original_halfedges = halfedges(w.graph());

  return make_range(halfedge_iterator(predicate, original_halfedges.first, original_halfedges.second),
                    halfedge_iterator(predicate, original_halfedges.second, original_halfedges.second));
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap>::size_type
num_halfedges(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return w.number_of_halfedges();
}

// FaceGraph
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::face_descriptor
face(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor h,
     const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(h));
  if(face(h, w.graph()) == boost::graph_traits<Graph>::null_face()) // h is a border hafedge
    return boost::graph_traits< CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::null_face();
  else if(w.is_in_cc(face(h,w.graph())))
    return face(h,w.graph());
  else
    return boost::graph_traits< CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::null_face();
}

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::halfedge_descriptor
halfedge(typename boost::graph_traits< Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::face_descriptor f,
         const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  CGAL_assertion(w.is_in_cc(f));
  return halfedge(f,w.graph());
}


template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
Iterator_range<typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::face_iterator>
faces(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  typedef typename boost::graph_traits<Face_filtered_graph<Graph, FIMap, VIMap, HIMap> >::face_iterator face_iterator;
  typedef typename boost::graph_traits<Graph >::face_iterator g_face_iterator;

  typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap> ::Is_simplex_valid predicate(&w);
  std::pair<g_face_iterator, g_face_iterator> original_faces = faces(w.graph());

  return make_range(face_iterator(predicate, original_faces.first, original_faces.second),
                    face_iterator(predicate, original_faces.second, original_faces.second));
}



template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
typename Face_filtered_graph<Graph, FIMap, VIMap, HIMap>::size_type
num_faces(const Face_filtered_graph<Graph, FIMap, VIMap, HIMap> & w)
{
  return w.number_of_faces();
}

template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
typename boost::property_map<Graph, PropertyTag >::const_type
get(PropertyTag ptag, const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return get(ptag, w.graph());
}

template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
typename boost::property_map<Graph, PropertyTag >::type
get(PropertyTag ptag, Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return get(ptag, w.graph());
}

//specializations for indices
template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap>
typename boost::property_map<Face_filtered_graph<Graph, FIMap, VIMap, HIMap>, CGAL::face_index_t >::type
get(CGAL::face_index_t, const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.get_face_index_map();
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap>
typename boost::property_map<Face_filtered_graph<Graph, FIMap, VIMap, HIMap>, boost::vertex_index_t >::type
get(boost::vertex_index_t, const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.get_vertex_index_map();
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap>
typename boost::property_map<Face_filtered_graph<Graph, FIMap, VIMap, HIMap>, CGAL::halfedge_index_t >::type
get(CGAL::halfedge_index_t, const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w)
{
  return w.get_halfedge_index_map();
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type
get(PropertyTag ptag,
    const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w,
    const typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::key_type& k)
{
  return get(ptag, w.graph(), k);
}


template <class Graph,
          typename FIMap,
          typename VIMap,
          typename HIMap,
          class PropertyTag>
void
put(PropertyTag ptag, const Face_filtered_graph<Graph, FIMap, VIMap, HIMap>& w,
    const typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::key_type& k,
    typename boost::property_traits<typename boost::property_map<Graph,PropertyTag>::type>::value_type& v)
{
  put(ptag, w.graph(), k, v);
}

}//end namespace CGAL

namespace boost {
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap,
         typename PropertyTag>
struct property_map<CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap>,PropertyTag> {
  typedef typename boost::property_map<Graph, PropertyTag >::type type;
  typedef typename boost::property_map<Graph, PropertyTag >::const_type const_type;
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap,
         typename PropertyTag>
struct graph_has_property<CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap>, PropertyTag>
    : graph_has_property<Graph, PropertyTag> {};


//specializations for indices
template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct property_map<CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap>, CGAL::face_index_t>{
  typedef typename CGAL::Property_map_binder< FIMap,
  typename CGAL::Pointer_property_map< typename boost::property_traits< FIMap >::value_type >::type > type;
  typedef type const_type;
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>
struct property_map<CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap>, boost::vertex_index_t>{
  typedef typename CGAL::Property_map_binder< VIMap,
  typename CGAL::Pointer_property_map< typename boost::property_traits< VIMap >::value_type >::type > type;
  typedef type const_type;
};

template<typename Graph,
         typename FIMap,
         typename VIMap,
         typename HIMap>

struct property_map<CGAL::Face_filtered_graph<Graph, FIMap, VIMap, HIMap>, CGAL::halfedge_index_t>{
  typedef typename CGAL::Property_map_binder< HIMap,
  typename CGAL::Pointer_property_map< typename boost::property_traits< HIMap >::value_type >::type > type;
  typedef type const_type;
};
}// namespace boost
#endif // CGAL_BOOST_GRAPH_FACE_FILTERED_GRAPH_H
