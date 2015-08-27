// Copyright (c) 2011, 2015 GeometryFactory (France).
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
//
//
// Author(s)     : Sebastien Loriot and Andreas Fabri

#ifndef CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H
#define CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H

#include<set>
#include<vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/assertions.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/tuple.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/internal/corefinement/Polyhedron_constness_types.h>
#include <CGAL/Default.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif


namespace CGAL {
 namespace internal{
   namespace corefinement{

template <typename Polyhedron>
struct Compare_handle_ptr{
  typedef typename Polyhedron::Facet_const_handle  Facet_const_handle;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  
  bool operator()(Facet_const_handle f1,Facet_const_handle f2) const {
    return &(*f1) < &(*f2);
  }

  bool operator()(Vertex_const_handle v1,Vertex_const_handle v2) const {
    return &(*v1) < &(*v2);
  }
};

struct Dummy_true{
  template <typename T>
  bool operator()(T) const  {return true;}
};

template <typename Polyhedron,typename HDS=typename Polyhedron::HalfedgeDS>
class Build_polyhedron_subset : public ::CGAL::Modifier_base<HDS> {
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;       
  
  typedef typename HDS::Vertex::Point Point;
  std::list<Vertex_const_handle> points;
  std::list< std::vector<unsigned int> > facets;
    
  template <typename Facet_iterator>
  typename Polyhedron::Halfedge_const_handle get_facet_halfedge(Facet_iterator facet_it) const
  {
    return (*facet_it)->halfedge();
  }

  typename Polyhedron::Halfedge_const_handle get_facet_halfedge(typename Polyhedron::Facet_const_handle facet) const
  {
    return facet->halfedge();
  }
  
public:
  template <typename Facets_const_iterator>
  Build_polyhedron_subset(const Polyhedron&,Facets_const_iterator begin,Facets_const_iterator end) 
  {
    typedef std::map<Vertex_const_handle,unsigned int,Compare_handle_ptr<Polyhedron> > Vertices;
    Vertices vertices;
    unsigned int index=0;
    //get vertices and get face description relatively to the restricted set of vertices
    for (Facets_const_iterator it=begin;it!=end;++it)
    {
      Halfedge_const_handle start=get_facet_halfedge(it);
      Halfedge_const_handle curr=start;
      facets.push_back(std::vector<unsigned int>());
      std::vector<unsigned int>& indices = facets.back();
      do{
        bool is_new_vertex;
        typename Vertices::iterator it_vertex;
        ::CGAL::cpp11::tie(it_vertex,is_new_vertex)=vertices.insert(std::make_pair(curr->vertex(),index));
        if (is_new_vertex) {
          ++index;
          points.push_back(curr->vertex());
        }
        indices.push_back(it_vertex->second);
        curr=curr->next();
      }while(curr!=start);
    }
  }
  
  void operator()( HDS& hds) {
    ::CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface( points.size(), facets.size() );
    for (typename std::list<Vertex_const_handle>::iterator it=points.begin();it!=points.end();++it)
      B.add_vertex((*it)->point());
    for (typename std::list< std::vector<unsigned int> >::iterator
         it=facets.begin();it!=facets.end();++it)
    {
      B.begin_facet();
      for (std::vector<unsigned int>::iterator it_i=it->begin();it_i!=it->end();++it_i)
        B.add_vertex_to_facet(*it_i);
      B.end_facet();
    }
    B.end_surface();
  }
};



template <typename Polyhedron,typename Adjacency_criterium,typename Face_to_UF_handle_map,typename Result>
void extract_connected_components(
  Polyhedron& P,
  const Adjacency_criterium& adjacent,
  CGAL::Union_find<typename internal_IOP::Polyhedron_types_with_mpl<Polyhedron>::Facet_handle>& uf,
  Face_to_UF_handle_map& map_f2h,
  Result&  result
  )
{
  typedef typename internal_IOP::Polyhedron_types_with_mpl<Polyhedron>::Facet_handle Facet_handle;
  typedef typename internal_IOP::Polyhedron_types_with_mpl<Polyhedron>::Facet_iterator Facet_iterator;
  typedef typename internal_IOP::Polyhedron_types_with_mpl<Polyhedron>::Halfedge_handle Halfedge_handle;
  typedef ::CGAL::Union_find<Facet_handle> UF;
  typedef typename UF::handle UF_handle;
  typedef typename UF::iterator UF_iterator;

//init union-find: each facet is in its own set  
  for (Facet_iterator it=P.facets_begin();it!=P.facets_end();++it){
    map_f2h.insert(std::make_pair(it,uf.make_set(it)));
  }
//merge 2 facets if they share a common edge  
  for (Facet_iterator it=P.facets_begin();it!=P.facets_end();++it){
    Facet_handle facet=it;
    
    UF_handle current=map_f2h.find(it)->second;
    std::vector<Halfedge_handle> neighbors;
    Halfedge_handle hedge=facet->halfedge();
    do
    {
      neighbors.push_back( hedge->opposite() );
      hedge=hedge->next();
    }
    while(hedge!=facet->halfedge());

    std::size_t nb_edges=neighbors.size();
    for (std::size_t i=0;i<nb_edges;++i){
      if ( neighbors[i]->is_border() ) continue;
      UF_handle neigh=map_f2h.find(neighbors[i]->facet())->second;
      if ( adjacent(neighbors[i]) && !uf.same_set(current,neigh) ){
        uf.unify_sets(current,neigh);
      }
    }
  }
  
//recover merged sets
  for (UF_iterator it=uf.begin();it!=uf.end();++it){
    UF_handle master=uf.find(it);
    result[*master].push_back(*it);
  }
}

template <typename Polyhedron,typename Adjacency_criterium,typename Output_iterator>
void extract_connected_components(const Polyhedron& P,const Adjacency_criterium& adjacent,Output_iterator out)
{
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef ::CGAL::Union_find<Facet_const_handle> UF;
  typedef typename UF::handle UF_handle;
  typedef std::map<Facet_const_handle,std::list<Facet_const_handle>,Compare_handle_ptr<Polyhedron> > Result;
  typedef std::map<Facet_const_handle,UF_handle,Compare_handle_ptr<Polyhedron> > Facet_to_handle_map;
  
  UF uf;
  Facet_to_handle_map map_f2h;
  Result result;
  
  extract_connected_components(P,adjacent,uf,map_f2h,result);
  
  for (typename Result::iterator it=result.begin();it!=result.end();++it)
  {
    typedef std::list<Facet_const_handle> Facets;
    const Facets& facets=it->second;
    Polyhedron new_poly;
    Build_polyhedron_subset<Polyhedron> modifier(new_poly,facets.begin(),facets.end());
    new_poly.delegate(modifier);
    *out++=new_poly;
  }
}

template <typename Polyhedron, typename Adjacency_criterium, typename Face_marker>
void mark_connected_components(Polyhedron& P, const Adjacency_criterium& adjacent, Face_marker& face_marker)
{
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef ::CGAL::Union_find<Facet_handle> UF;
  typedef typename UF::handle UF_handle;
  typedef std::map<Facet_handle,std::list<Facet_handle>,Compare_handle_ptr<Polyhedron> > Result;
  typedef std::map<Facet_handle,UF_handle,Compare_handle_ptr<Polyhedron> > Facet_to_handle_map;
  
  UF uf;
  Facet_to_handle_map map_f2h;
  Result result;
  
  extract_connected_components(P,adjacent,uf,map_f2h,result);
  
  for (typename Result::iterator it=result.begin();it!=result.end();++it)
  {
    face_marker.start_new_connected_component();
    typedef std::list<Facet_handle> Facets;
    const Facets& facets=it->second;
    face_marker.mark(facets.begin(),facets.end());
  }
}

template <typename Polyhedron, typename Adjacency_criterium, typename Face_marker, typename OutputIterator>
OutputIterator
mark_connected_components(Polyhedron& P, const Adjacency_criterium& adjacent, Face_marker& face_marker, OutputIterator out)
{
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef ::CGAL::Union_find<Facet_handle> UF;
  typedef typename UF::handle UF_handle;
  typedef std::map<Facet_handle,std::list<Facet_handle>,Compare_handle_ptr<Polyhedron> > Result;
  typedef std::map<Facet_handle,UF_handle,Compare_handle_ptr<Polyhedron> > Facet_to_handle_map;

  UF uf;
  Facet_to_handle_map map_f2h;
  Result result;

  extract_connected_components(P,adjacent,uf,map_f2h,result);

  for (typename Result::iterator it=result.begin();it!=result.end();++it)
  {
    face_marker.start_new_connected_component();
    typedef std::list<Facet_handle> Facets;
    const Facets& facets=it->second;
    face_marker.mark(facets.begin(),facets.end());
    *out++=*facets.begin();
  }
  return out;
}

template <typename Polyhedron,typename Output_iterator>
void extract_connected_components(const Polyhedron& P,Output_iterator out)
{
  extract_connected_components(P,Dummy_true(),out);
}

template <typename Polyhedron, typename Polyhedron_facet_index_map>
std::size_t
init_facet_indices(
  const Polyhedron& P,
  Polyhedron_facet_index_map facet_index_map)
{
  //init facet indices
  std::size_t index=0;
  for (typename Polyhedron::Facet_const_iterator fit=P.facets_begin(),
                                                 fit_end=P.facets_end();
                                                 fit!=fit_end; ++fit)
  {
    put(facet_index_map, fit, index++);
  }
  return index;
}

// alternative method by propagation
template <typename Polyhedron, typename Adjacency_criterium, typename Polyhedron_facet_index_map>
std::size_t
mark_connected_components_v2(
  const Polyhedron& P,
  const Adjacency_criterium& adjacent,
  Polyhedron_facet_index_map facet_index_map,
  std::vector<std::size_t>& patch_ids,
  std::vector<std::size_t>& patch_sizes)
{
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_handle;

  std::size_t max_id=(std::numeric_limits<std::size_t>::max)();
  patch_ids.clear();
  patch_ids.resize(P.size_of_facets(), max_id);

  //traversal of the facets to discover connected components
  std::size_t patch_id=0;
  for (typename Polyhedron::Facet_const_iterator fit=P.facets_begin(),
                                                 fit_end=P.facets_end();
                                                 fit!=fit_end; ++fit)
  {
    std::size_t index=get(facet_index_map, fit);
    if ( patch_ids[index]==max_id )
    {
      patch_sizes.push_back(0);
      patch_ids[index]=patch_id;// set patch id
      ++(patch_sizes.back());
      std::vector<typename Polyhedron::Halfedge_const_handle> queue;
      if ( adjacent(fit->halfedge()) )
        queue.push_back( fit->halfedge()->opposite() );
      if ( adjacent(fit->halfedge()->next()) )
        queue.push_back( fit->halfedge()->next()->opposite() );
      if ( adjacent(fit->halfedge()->next()->next()) )
        queue.push_back( fit->halfedge()->next()->next()->opposite() );
      while (!queue.empty())
      {
        Halfedge_handle h=queue.back();
        queue.pop_back();
        index=get(facet_index_map, h->facet());
        if ( patch_ids[index]!=max_id ) continue;
        patch_ids[index]=patch_id;
        ++(patch_sizes.back());
        if ( adjacent(h->next()) )
          queue.push_back( h->next()->opposite() );
        if ( adjacent(h->next()->next()) )
          queue.push_back( h->next()->next()->opposite() );
      }
      ++patch_id;
    }
  }

  return patch_id;
}

} } // end of namespace internal::corefinement

namespace Polygon_mesh_processing{

  namespace internal {
    struct MoreSecond {
      typedef std::pair<std::size_t,std::size_t> T;
      bool operator()(const T& a, const T& b) const {
        return a.second > b.second;
      }
    };

    // A property map 
    template <typename G>
    struct No_constraint {
      friend bool get(No_constraint<G>, typename boost::graph_traits<G>::edge_descriptor)
      {
        return false;
      }
    };

    // A functor
    template <typename G, typename EdgeConstraintMap = No_constraint<G> >
    struct No_border {
      No_border()
      {}

      No_border(const G & g, EdgeConstraintMap ecm = EdgeConstraintMap())
        : g(&g), ecm(ecm)
      {}

      bool operator()(typename boost::graph_traits<G>::edge_descriptor e) const {
        if (!is_border(e, *g)){
          return !get(ecm, e);
        }
        return false;
      }

      const G* g;
      EdgeConstraintMap ecm;
    };

}// namespace internal

/*!
 * \ingroup PkgPolygonMeshProcessing
 *  discovers all the faces in the same connected component as `seed_face` and records them in `out`.
 * `seed_face` will also be added in `out`.
 *  Two faces are recorded in the same connected component if they share an edge that is not marked as constrained.

 *  \tparam PolygonMesh a model of `FaceGraph`
 *  \tparam FaceOutputIterator a model of `OutputIterator` that accepts
        faces of type
        `boost::graph_traits<PolygonMesh>::%face_descriptor`.
 *  \tparam NamedParameters a sequence of \ref namedparameters
 *
 *  \param seed_face a face of `pmesh` from which exploration starts to detect the connected component
           that contains it
 *  \param pmesh the polygon mesh
 *  \param out the output iterator that collects faces from the same connected component as `seed_face`
 *  \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{edge_is_constrained_map}  a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 *  \returns the output iterator.
 *
 */
template <typename PolygonMesh
          , typename FaceOutputIterator
          , typename NamedParameters
          >
FaceOutputIterator
connected_component(typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face
                    , const PolygonMesh& pmesh
                    , FaceOutputIterator out
                    , const NamedParameters& np)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  typedef typename boost::lookup_named_param_def <
    CGAL::edge_is_constrained_t,
    NamedParameters,
    internal::No_constraint<PolygonMesh>//default
  > ::type                                               EdgeConstraintMap;
  EdgeConstraintMap ecmap
    = choose_param(get_param(np, edge_is_constrained), EdgeConstraintMap());

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  std::set<face_descriptor> already_processed;
  std::vector< face_descriptor > stack;
  stack.push_back(seed_face);
  while (!stack.empty())
    {
      seed_face=stack.back();
      stack.pop_back();
      if (!already_processed.insert(seed_face).second) continue;
      *out++=seed_face;
      BOOST_FOREACH(halfedge_descriptor hd,
                    CGAL::halfedges_around_face(halfedge(seed_face, pmesh), pmesh) )
      {
        if(! get(ecmap, edge(hd, pmesh))){
          face_descriptor neighbor = face( opposite(hd, pmesh), pmesh );
          if ( neighbor != boost::graph_traits<PolygonMesh>::null_face() )
            stack.push_back(neighbor);
        }
      }
    }
  return out;
}

template <typename PolygonMesh, typename OutputIterator>
OutputIterator
connected_component(typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face,
                    const PolygonMesh& pmesh,
                    OutputIterator out)
{
  return connected_component(seed_face, pmesh, out,
          CGAL::Polygon_mesh_processing::parameters::all_default());
}

/*!
 * \ingroup PkgPolygonMeshProcessing
 *  computes for each face the index of the corresponding connected component.
 *  Two faces are recorded in the same connected component if they share an edge that is not marked as constrained.
 *
 *  A property map for `CGAL::face_index_t` should be either available as an internal property map 
 *  to `pmesh` or provided as one of the \ref namedparameters.
 *
 *  \tparam PolygonMesh a model of `FaceListGraph`
 *  \tparam FaceComponentMap a model of `WritablePropertyMap` with
        `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
        `boost::face_index` as value type.
 *  \tparam NamedParameters a sequence of \ref namedparameters

 * \param pmesh the polygon mesh
 * \param fcm the property map with indices of components associated to faces in `pmesh`
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 *  \returns the number of connected components.
 */

template <typename PolygonMesh
        , typename FaceComponentMap
        , typename NamedParameters
>
typename boost::property_traits<FaceComponentMap>::value_type
connected_components(const PolygonMesh& pmesh,
                     const FaceComponentMap& fcm,
                     const NamedParameters& np)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  typedef typename boost::lookup_named_param_def <
    CGAL::edge_is_constrained_t,
    NamedParameters,
    internal::No_constraint<PolygonMesh>//default
  > ::type                                               EdgeConstraintMap;
  EdgeConstraintMap ecmap
    = choose_param(get_param(np, edge_is_constrained), EdgeConstraintMap());

  typedef Dual<PolygonMesh>                              Dual;
  typedef boost::filtered_graph<Dual,
    internal::No_border<PolygonMesh,EdgeConstraintMap> > FiniteDual;
  Dual dual(pmesh);

  FiniteDual finite_dual(dual,
    internal::No_border<PolygonMesh, EdgeConstraintMap>(pmesh, ecmap));

  return boost::connected_components(finite_dual,
    fcm,
    boost::vertex_index_map(
      choose_const_pmap(get_param(np, boost::face_index),
                        pmesh,
                        boost::face_index)
    )
  );
}

template <typename PolygonMesh, typename FaceComponentMap>
typename boost::property_traits<FaceComponentMap>::value_type
connected_components(const PolygonMesh& pmesh,
                     const FaceComponentMap& fcm)
{

  return CGAL::Polygon_mesh_processing::connected_components(pmesh, fcm,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}


template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void keep_connected_components(PolygonMesh& pmesh
                              , const ComponentRange& components_to_keep
                              , const FaceComponentMap& fcm
                              , const NamedParameters& np);

/*!
 * \ingroup PkgPolygonMeshProcessing
 *  removes the small connected components and all the isolated vertices.
 *  Keep `nb_components_to_keep` largest connected components. 
 *  Two faces are considered in the same connected component if they share an edge that is not marked as constrained.
 *
 * Property maps for `CGAL::face_index_t` and `CGAL::vertex_index_t`
 * should be either available as internal property maps 
 * to `pmesh` or provided as \ref namedparameters.
 *
 * \tparam PolygonMesh a model of `FaceListGraph`
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param pmesh the polygon mesh
 * \param nb_components_to_keep the number of components to be kept
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 *  \return the number of connected components removed (ignoring isolated vertices).
 */
template <typename PolygonMesh
        , typename NamedParameters>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh
                                            , std::size_t nb_components_to_keep
                                            , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;

  using boost::choose_param;
  using boost::get_param;
  using boost::choose_const_pmap;

  //FaceIndexMap
  typedef typename GetFaceIndexMap<PM, NamedParameters>::type FaceIndexMap;
  FaceIndexMap fim = choose_const_pmap(get_param(np, boost::face_index),
                                       pmesh,
                                       boost::face_index);

  //vector_property_map
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fim);
  std::size_t num = connected_components(pmesh, face_cc, np);

  if((num == 1)|| (nb_components_to_keep > num) )
    return 0;

  std::vector< std::pair<std::size_t, std::size_t> > component_size(num);

  for(std::size_t i=0; i < num; i++)
    component_size[i] = std::make_pair(i,0);

  BOOST_FOREACH(face_descriptor f, faces(pmesh))
    ++component_size[face_cc[f]].second;

  // we sort the range [0, num) by component size
  std::sort(component_size.begin(), component_size.end(), internal::MoreSecond());
  std::vector<std::size_t> cc_to_keep;
  for(std::size_t i=0; i<nb_components_to_keep; ++i)
    cc_to_keep.push_back( component_size[i].first );

  keep_connected_components(pmesh, cc_to_keep, face_cc, np);

  return num - nb_components_to_keep;
}

template <typename PolygonMesh>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh,
                                              std::size_t nb_components_to_keep)
{
  return keep_largest_connected_components(pmesh,
    nb_components_to_keep,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void keep_or_remove_connected_components(PolygonMesh& pmesh
                                        , const ComponentRange& components_to_keep
                                        , const FaceComponentMap& fcm
                                        , bool  keep
                                        , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  using boost::choose_param;
  using boost::get_param;
  using boost::choose_const_pmap;

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor   face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_iterator     face_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_iterator   vertex_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_iterator     edge_iterator;

  //VertexIndexMap
  typedef typename GetVertexIndexMap<PM, NamedParameters>::type VertexIndexMap;
  VertexIndexMap vim = choose_const_pmap(get_param(np, boost::vertex_index),
                                         pmesh,
                                         boost::vertex_index);

  std::set<std::size_t> cc_to_keep;
  BOOST_FOREACH(std::size_t i, components_to_keep)
    cc_to_keep.insert(i);

  boost::vector_property_map<bool, VertexIndexMap> keep_vertex(vim);
  BOOST_FOREACH(vertex_descriptor v, vertices(pmesh)){
    keep_vertex[v] = false;
  }
  BOOST_FOREACH(face_descriptor f, faces(pmesh)){
    if (cc_to_keep.find(get(fcm,f)) != cc_to_keep.end())
      put(fcm, f, keep ? 1 : 0);
    else
      put(fcm, f, keep ? 0 : 1);
  }

  BOOST_FOREACH(face_descriptor f, faces(pmesh)){
    if (get(fcm, f) == 1){
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(f, pmesh), pmesh)){
        vertex_descriptor v = target(h, pmesh);
        keep_vertex[v] = true;
      }
    }
  }

  edge_iterator eb, ee;
  for (boost::tie(eb, ee) = edges(pmesh); eb != ee;)
  {
    edge_descriptor e = *eb;
    ++eb;
    vertex_descriptor v = source(e, pmesh);
    vertex_descriptor w = target(e, pmesh);
    halfedge_descriptor h = halfedge(e, pmesh);
    halfedge_descriptor oh = opposite(h, pmesh);
    if (!keep_vertex[v] && !keep_vertex[w]){
      // don't care about connectivity
      // As vertices are not kept the faces and vertices will be removed later
      remove_edge(e, pmesh);
    }
    else if (keep_vertex[v] && keep_vertex[w]){
      face_descriptor fh = face(h, pmesh), ofh = face(oh, pmesh);
      if (is_border(h, pmesh) && is_border(oh, pmesh)){
#ifdef CGAL_CC_DEBUG
        std::cerr << "null_face on both sides of " << e << " is kept\n";
#endif
      }
      else if ((is_border(oh, pmesh) && get(fcm,fh)) ||
        (is_border(h, pmesh) && get(fcm,ofh)) ||
        (!is_border(oh, pmesh) && !is_border(h, pmesh) && get(fcm,fh) && get(fcm,ofh))){
        // do nothing
      }
      else if (!is_border(h, pmesh) && get(fcm,fh) && !is_border(oh, pmesh) && !get(fcm,ofh)){
        set_face(oh, boost::graph_traits<PolygonMesh>::null_face(), pmesh);
      }
      else if (!is_border(h, pmesh) && !get(fcm,fh) && !is_border(oh, pmesh) && get(fcm,ofh)){
        set_face(h, boost::graph_traits<PolygonMesh>::null_face(), pmesh);
      }
      else {
        // no face kept
        CGAL_assertion((is_border(h, pmesh) || !get(fcm,fh)) && (is_border(oh, pmesh) || !get(fcm,ofh)));
        // vertices pointing to e must change their halfedge
        if (halfedge(v, pmesh) == oh){
          set_halfedge(v, prev(h, pmesh), pmesh);
        }
        if (halfedge(w, pmesh) == h){
          set_halfedge(w, prev(oh, pmesh), pmesh);
        }
        // shortcut the next pointers as e will be removed
        set_next(prev(h, pmesh), next(oh, pmesh), pmesh);
        set_next(prev(oh, pmesh), next(h, pmesh), pmesh);
        remove_edge(e, pmesh);
      }
    }
    else if (keep_vertex[v]){
      if (halfedge(v, pmesh) == oh){
        set_halfedge(v, prev(h, pmesh), pmesh);
      }
      set_next(prev(h, pmesh), next(oh, pmesh), pmesh);
      remove_edge(e, pmesh);
    }
    else {
      CGAL_assertion(keep_vertex[w]);
      if (halfedge(w, pmesh) == h){
        set_halfedge(w, prev(oh, pmesh), pmesh);
      }
      set_next(prev(oh, pmesh), next(h, pmesh), pmesh);
      remove_edge(e, pmesh);
    }
  }

  face_iterator fb, fe;
  // We now can remove all vertices and faces not marked as kept
  for (boost::tie(fb, fe) = faces(pmesh); fb != fe;){
    face_descriptor f = *fb;
    ++fb;
    if (get(fcm,f) != 1){
      remove_face(f, pmesh);
    }
  }
  vertex_iterator b, e;
  for (boost::tie(b, e) = vertices(pmesh); b != e;){
    vertex_descriptor v = *b;
    ++b;
    if (!keep_vertex[v]){
      remove_vertex(v, pmesh);
    }
  }
}

/*!
* \ingroup keep_connected_components_grp
* keeps the connected components designated by theirs ids in `components_to_keep`,
* and removes the other connected components as well as all the isolated vertices.
* The connected component id of a face is given by `fcm`.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* Property maps for `CGAL::vertex_index_t`
* should be either available as internal property map
* to `pmesh` or provided as \ref namedparameters.
*
* \tparam PolygonMesh a model of `FaceListGraph`
* \tparam NamedParameters a sequence of \ref namedparameters
* \tparam ComponentRange a range of ids convertible to `std::size`
* \tparam FaceComponentMap a model of `ReadWritePropertyMap` with
*         `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
*         `boost::face_index` as value type.
*
* \param components_to_keep the range of ids of connected components to keep
* \param pmesh the polygon mesh
* \param fcm the property map with indices of components associated to faces in `pmesh`.
*        After calling this function, the values of `fcm` are undefined.
* \param np optional \ref namedparameters described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void keep_connected_components(PolygonMesh& pmesh
                              , const ComponentRange& components_to_keep
                              , const FaceComponentMap& fcm
                              , const NamedParameters& np)
{
  keep_or_remove_connected_components(pmesh, components_to_keep, fcm, true, np);
}

/*!
* \ingroup remove_connected_components_grp
* Removes in `pmesh` the connected components designated by theirs ids
* in `components_to_remove` as well as all the isolated vertices.
* The connected component id of a face is given by `fcm`.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* Property maps for `CGAL::vertex_index_t`
* should be either available as internal property map
* to `pmesh` or provided as \ref namedparameters.
*
*
* \tparam PolygonMesh a model of `FaceListGraph`
* \tparam NamedParameters a sequence of \ref namedparameters
* \tparam ComponentRange a range of ids convertible to `std::size`
* \tparam FaceComponentMap a model of `ReadWritePropertyMap` with
*         `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
*         `boost::face_index` as value type.
*
* \param components_to_remove the range of ids of connected components to remove
* \param pmesh the polygon mesh
* \param fcm the property map with indices of components associated to faces in `pmesh`.
*        After calling this function, the values of `fcm` are undefined.
* \param np optional \ref namedparameters described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void remove_connected_components(PolygonMesh& pmesh
                                , const ComponentRange& components_to_remove
                                , const FaceComponentMap& fcm
                                , const NamedParameters& np)
{
  if (components_to_remove.empty()) return;
  keep_or_remove_connected_components(pmesh, components_to_remove, fcm, false, np);
}

/*!
* \ingroup remove_connected_components_grp
*  keeps the connected components not designated by the faces in `components_to_remove`,
*  and removes the other connected components and all the isolated vertices.
*  Two faces are considered in the same connected component if they share an edge that is not marked as constrained.
*
* Property maps for `CGAL::face_index_t` and `CGAL::vertex_index_t`
* should be either available as internal property maps
* to `pmesh` or provided as \ref namedparameters.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* \tparam PolygonMesh a model of `FaceListGraph`
* \tparam NamedParameters a sequence of \ref namedparameters
* \tparam FaceRange a range of `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         indicating the connected components to be removed.
*
* \param components_to_remove a face range, including one face or more on each component to be removed
* \param pmesh the polygon mesh
* \param np optional \ref namedparameters described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
*    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename FaceRange
        , typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void remove_connected_components(PolygonMesh& pmesh
                                , const FaceRange& components_to_remove
                                , const CGAL_PMP_NP_CLASS& np)
{
  if (components_to_remove.empty()) return;
  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;
  using boost::choose_param;
  using boost::get_param;
  using boost::choose_const_pmap;

  //FaceIndexMap
  typedef typename GetFaceIndexMap<PM, CGAL_PMP_NP_CLASS>::type FaceIndexMap;
  FaceIndexMap fim = choose_const_pmap(get_param(np, boost::face_index),
                                       pmesh,
                                       boost::face_index);

  //vector_property_map
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fim);

  connected_components(pmesh, face_cc, np);

  std::vector<std::size_t> cc_to_remove;
  BOOST_FOREACH(face_descriptor f, components_to_remove)
    cc_to_remove.push_back( face_cc[f] );

  remove_connected_components(pmesh, cc_to_remove, face_cc, np);
}

/*!
* \ingroup keep_connected_components_grp
*  keeps the connected components designated by the faces in `components_to_keep`,
*  and removes the other connected components and all the isolated vertices.
*  Two faces are considered in the same connected component if they share an edge that is not marked as constrained.
*
* Property maps for `CGAL::face_index_t` and `CGAL::vertex_index_t`
* should be either available as internal property maps
* to `pmesh` or provided as \ref namedparameters.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* \tparam PolygonMesh a model of `FaceListGraph`
* \tparam NamedParameters a sequence of \ref namedparameters
* \tparam FaceRange a range of `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         indicating the connected components to be kept.
*
* \param pmesh the polygon mesh
* \param components_to_keep a face range, including one face or more on each component to be kept
* \param np optional \ref namedparameters described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
*    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename FaceRange
        , typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void keep_connected_components(PolygonMesh& pmesh
                             , const FaceRange& components_to_keep
                             , const CGAL_PMP_NP_CLASS& np)
{
  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;

  using boost::choose_param;
  using boost::get_param;
  using boost::choose_const_pmap;

  //FaceIndexMap
  typedef typename GetFaceIndexMap<PM, CGAL_PMP_NP_CLASS>::type FaceIndexMap;
  FaceIndexMap fim = choose_const_pmap(get_param(np, boost::face_index),
                                       pmesh,
                                       boost::face_index);

  //vector_property_map
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fim);

  connected_components(pmesh, face_cc, np);

  std::vector<std::size_t> cc_to_keep;
  BOOST_FOREACH(face_descriptor f, components_to_keep)
    cc_to_keep.push_back( face_cc[f] );

  keep_connected_components(pmesh, cc_to_keep, face_cc, np);
}

// non-documented overloads so that named parameters can be omitted

template <typename PolygonMesh, typename FaceRange>
void remove_connected_components(PolygonMesh& pmesh
                                , const FaceRange& components_to_remove)
{
  remove_connected_components(pmesh, components_to_remove,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap>
void keep_connected_components(PolygonMesh& pmesh
                              , const ComponentRange& components_to_keep
                              , const FaceComponentMap& fcm)
{
  keep_connected_components(pmesh, components_to_keep, fcm,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap>
void remove_connected_components(PolygonMesh& pmesh
                                , const ComponentRange& components_to_remove
                                , const FaceComponentMap& fcm )
{
    remove_connected_components(pmesh, components_to_remove, fcm,
      CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <typename PolygonMesh, typename FaceRange>
void keep_connected_components(PolygonMesh& pmesh
                             , const FaceRange& components_to_keep)
{
  keep_connected_components(pmesh, components_to_keep,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

} // namespace Polygon_mesh_processing

} // namespace CGAL

#endif //CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H
