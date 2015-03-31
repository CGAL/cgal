// Copyright (c) 2011 GeometryFactory (France).
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

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/tuple.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/internal/corefinement/Polyhedron_constness_types.h>
#include <CGAL/Default.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>

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
    struct LessFirst {
      typedef std::pair<std::size_t,std::size_t> T;
      bool operator()(const T& a, const T& b) const {
        return a.first < b.first;
      }
    };
    struct MoreSecond {
      typedef std::pair<std::size_t,std::size_t> T;
      bool operator()(const T& a, const T& b) const {
        return a.second > b.second;
      }
    };
  } // namespace internal

/*!
 * \ingroup PkgPolygonMeshProcessing
 *  discovers all the faces in the same connected component as `seed_face` and puts them in `out`.
 * `seed_face` will also be added in `out`.
 *  Two faces are put in the same connected component if they share an edge that is not marked as constrained.

 *  \tparam PolygonMesh a model of `FaceGraph`
 *  \tparam FaceOutputIterator a model of `OutputIterator` that accepts
        `boost::graph_traits<PolygonMesh>::%face_descriptor`s.
 *  \tparam EdgeConstraintMap a model of `ReadablePropertyMap` with
        `boost::graph_traits<PolygonMesh>::%edge_descriptor` as key type and
        `bool` as value type. It should be default-constructible.
 *
 *  \param seed_face a face of `pmesh` from which exploration starts to detect the connected component
           that contains it
 *  \param pmesh the polygon mesh
 *  \param out the output iterator that collects faces from the same connected component as `seed_face`
 *  \param ecmap the property map containing information about edges of `pmesh` being constrained or not.
      If this parameter is omitted,
      a default property map where no edge is constrained is provided.
 *  \returns the output iterator.
 *
 * \todo : add named parameters for ecmap?
 */
template <typename PolygonMesh
          , typename FaceOutputIterator
          , typename EdgeConstraintMap
          >
FaceOutputIterator
connected_component(typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face
                    , const PolygonMesh& pmesh
                    , FaceOutputIterator out
                    , EdgeConstraintMap ecmap
#ifdef DOXYGEN_RUNNING
                    = CGAL::Default()
#endif
                    )
{
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

namespace internal {

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
    if(! is_border(e,*g)){
      return ! get(ecm,e);
    }
    return false;
  }

  const G* g;
  EdgeConstraintMap ecm;
};

}// namespace internal

template <typename PolygonMesh, typename OutputIterator>
OutputIterator
connected_component(typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face,
                    const PolygonMesh& pmesh,
                    OutputIterator out)
{
  return connected_component(seed_face, pmesh, out,
              internal::No_constraint<PolygonMesh>());
}

/*!
 * \ingroup PkgPolygonMeshProcessing
 *  computes for each face the index of the connected component to which it belongs.
 *  Two faces are put in the same connected component if they share an edge that is not marked as constrained.
 *  \tparam PolygonMesh a model of `FaceListGraph`
 *  \tparam FaceComponentMap a model of `WritablePropertyMap` with
        `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
        `boost::face_index` as value type.
 *  \tparam EdgeConstraintMap a model of `ReadablePropertyMap` with
        `boost::graph_traits<PolygonMesh>::%edge_descriptor` as key type and
        `bool` as value type. It should be default-constructible.
 * \tparam FaceIndexMap a model of `ReadablePropertyMap` with
       `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
       ` CGAL::face_index_t` as value type.

 * \param pmesh the polygon mesh
 * \param fcm the property map with indices of components associated to faces in `pmesh`
 * \param ecmap the property map containing information about edges of `pmesh` being constrained or not.
          It defaults to a map with all value types being `false`.
 * \param fim the property map containing the index of each face of `pmesh`
 *
 *  \returns the number of connected components.
 */

template <typename PolygonMesh
        , typename FaceComponentMap
        , typename NamedParameters
>
typename boost::property_traits<FaceComponentMap>::value_type
connected_components(const PolygonMesh& pmesh,
                     FaceComponentMap& fcm,
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
                     FaceComponentMap& fcm)
{

  return CGAL::Polygon_mesh_processing::connected_components(pmesh, fcm,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

/*!
 * \ingroup PkgPolygonMeshProcessing
 *  erases the small connected components and the isolated vertices.
 *  Keep `nb_components_to_keep` largest connected components. 
 *  Two faces are put in the same connected component if they share an edge that is not marked as constrained.
 *
 * \tparam PolygonMesh a model of `FaceListGraph`
 * \tparam EdgeConstraintMap a model of `ReadablePropertyMap` with
       `boost::graph_traits<PolygonMesh>::%edge_descriptor` as key type and
       `bool` as value type. It should be default-constructible.
 * \tparam VertexIndexMap a model of `ReadablePropertyMap` with
       `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
       `CGAL::vertex_index_t` as value type
 * \tparam FaceIndexMap a model of `ReadablePropertyMap` with
       `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
       `CGAL::face_index_t` as value type.

 * \param pmesh the polygon mesh
 * \param nb_components_to_keep
 * \param ecmap the property map containing information about edges of `pmesh` being constrained or not.
               It defaults to a map with all value types being `false`.
 * \param vim the property map containing the index of each vertex of `pmesh`
 * \param fim the property map containing the index of each face of `pmesh`

 *  \return the number of connected components erased (ignoring isolated vertices).
 */
template <typename PolygonMesh
        , typename EdgeConstraintMap
        , typename VertexIndexMap
#ifdef DOXYGEN_RUNNING
          = typename boost::property_map<PolygonMesh, CGAL::vertex_index_t>::type
#endif
        , typename FaceIndexMap
#ifdef DOXYGEN_RUNNING
          = typename boost::property_map<PolygonMesh, CGAL::face_index_t>::type
#endif
>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh
                                              , std::size_t nb_components_to_keep
                                              , EdgeConstraintMap ecmap
#ifdef DOXYGEN_RUNNING
                                              = CGAL::Default()
#endif
                                              , VertexIndexMap vim
#ifdef DOXYGEN_RUNNING
                     = get(CGAL::vertex_index_t, pmesh)
#endif
                                              ,FaceIndexMap fim
#ifdef DOXYGEN_RUNNING
                     = get(CGAL::face_index_t, pmesh)
#endif
)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_iterator face_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_iterator edge_iterator;
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fim);
  
  std::size_t num = connected_components(pmesh, face_cc,
    CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(ecmap).
    face_index_map(fim));

  if((num == 1)|| (nb_components_to_keep > num) ){
    return 0;
  }
  boost::vector_property_map<bool, VertexIndexMap> keep_vertex(vim); 
  BOOST_FOREACH(vertex_descriptor v, vertices(pmesh)){
    keep_vertex[v] = false;
  }
  typedef std::pair<std::size_t, std::size_t> Component;
  std::vector<Component> component_size(num);

  for(std::size_t i=0; i < num; i++){
    component_size[i] = std::make_pair(i,0);
  }
  BOOST_FOREACH(face_descriptor f, faces(pmesh)){
    component_size[face_cc[f]].second++;
  }

#ifdef CGAL_CC_DEBUG
  for(int i=0; i < component_size.size(); ++i){
    std::cerr << "component " << i << " has " << component_size[i].second << " faces\n";
  }
#endif

  // we have to sort the range [0, num) by component size
  internal::MoreSecond ls;
  std::sort(component_size.begin(), component_size.end(), ls);
  for(std::size_t i=0; i < num; i++){
    component_size[i].second = (i < nb_components_to_keep)?1:0;
  }
  internal::LessFirst lsinv;
  std::sort(component_size.begin(), component_size.end(), lsinv);

#ifdef CGAL_CC_DEBUG
  for(std::size_t i=0; i < num; i++){
    std::cout << i  << " " << component_size[i].first  << " " << component_size[i].second << std::endl;
  }
#endif

  BOOST_FOREACH(face_descriptor f, faces(pmesh)){
    face_cc[f] = component_size[face_cc[f]].second;
#ifdef CGAL_CC_DEBUG
    if(face_cc[f] == 1){
      std::cerr << "keep " << f << std::endl;
    }
#endif
  }
  // Now face_cc[f] == 1 means that we want to keep the face

  BOOST_FOREACH(face_descriptor f, faces(pmesh)){
    if(face_cc[f] == 1){
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(f,pmesh),pmesh)){
          vertex_descriptor v = target(h,pmesh);
          keep_vertex[v] = true;
#ifdef CGAL_CC_DEBUG
          std::cout << "keep vertex "<< v << std::endl;
#endif
        }
    }
  }

  edge_iterator eb, ee;
  for(boost::tie(eb,ee)= edges(pmesh);eb!= ee;){
    edge_descriptor e = *eb;
    ++eb;
    vertex_descriptor v = source(e,pmesh);
    vertex_descriptor w = target(e,pmesh);
    halfedge_descriptor h = halfedge(e,pmesh);
    halfedge_descriptor oh = opposite(h,pmesh);
    if(! keep_vertex[v] && ! keep_vertex[w]){
      // don't care about connectivity
      // As vertices are not kept the faces and vertices will be removed later
      remove_edge(e,pmesh);
    } else if( keep_vertex[v] &&  keep_vertex[w]){
      face_descriptor fh = face(h,pmesh), ofh = face(oh,pmesh);
      if(is_border(h,pmesh) && is_border(oh,pmesh)){
#ifdef CGAL_CC_DEBUG
        std::cerr << "null_face on both sides of " << e << " is kept\n";
#endif
      } else if( (face_cc[fh] && is_border(oh,pmesh)) ||
                 (face_cc[ofh] && is_border(h,pmesh)) ||
                 (face_cc[fh] && face_cc[ofh]) ){
        // do nothing
      } else if(face_cc[fh] && ! face_cc[ofh]){
        set_face(oh, boost::graph_traits<PolygonMesh>::null_face(), pmesh);
      } else if(! face_cc[fh] && face_cc[ofh]){
        set_face(h, boost::graph_traits<PolygonMesh>::null_face(), pmesh);
      } else {
        // no face kept
        assert( ! face_cc[fh] && ! face_cc[ofh]);
        // vertices pointing to e must change their halfedge
        if(halfedge(v,pmesh) == oh){
          set_halfedge(v,prev(h,pmesh),pmesh);
        }
        if(halfedge(w,pmesh) == h){
          set_halfedge(w,prev(oh,pmesh),pmesh);
        }
        // shortcut the next pointers as e will be removed
        set_next(prev(h,pmesh), next(oh,pmesh),pmesh);
        set_next(prev(oh,pmesh), next(h,pmesh),pmesh);
        remove_edge(e,pmesh);
      }
    } else if( keep_vertex[v] ){
      if(halfedge(v,pmesh) == oh){
        set_halfedge(v,prev(h,pmesh),pmesh);
      }
      set_next(prev(h,pmesh), next(oh,pmesh),pmesh);
        remove_edge(e,pmesh);
    } else {
      assert (keep_vertex[w]);
      if(halfedge(w,pmesh) == h){
        set_halfedge(w,prev(oh,pmesh),pmesh);
      }
      set_next(prev(oh,pmesh), next(h,pmesh),pmesh);
        remove_edge(e,pmesh);
    }
  }
 
  face_iterator fb, fe;
  // We now can remove all vertices and faces not marked as kept
  for(boost::tie(fb,fe)=faces(pmesh); fb!=fe;){
    face_descriptor f = *fb;
    ++fb;
    if(face_cc[f] != 1){
      remove_face(f,pmesh);
    }
  }
  vertex_iterator b,e;
  for(boost::tie(b,e)=vertices(pmesh); b!=e;){
    vertex_descriptor v = *b;
    ++b;
    if(! keep_vertex[v]){
      remove_vertex(v,pmesh);
    }
  }
  return num - nb_components_to_keep;
}

template <typename PolygonMesh
        , typename VertexIndexMap
        , typename FaceIndexMap>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh
                                            , std::size_t nb_components_to_keep
                                            , CGAL::Default
                                            , VertexIndexMap vim
                                            , FaceIndexMap fim)
{
  return keep_largest_connected_components(pmesh, nb_components_to_keep,
                                           internal::No_constraint<PolygonMesh>(),
                                           vim,
                                           fim);
}

template <typename PolygonMesh, typename EdgeConstraintMap, typename VertexIndexMap>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh,
                                              std::size_t nb_components_to_keep,
                                              EdgeConstraintMap ecmap,
                                              VertexIndexMap vim)
{
  return keep_largest_connected_components(pmesh, nb_components_to_keep, ecmap, vim,
                                           get(boost::face_index, pmesh));
}

template <typename PolygonMesh, typename EdgeConstraintMap>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh,
                                              std::size_t nb_components_to_keep,
                                              EdgeConstraintMap ecmap)
{
  return keep_largest_connected_components(pmesh, nb_components_to_keep, ecmap,
                                           get(boost::vertex_index, pmesh),
                                           get(boost::face_index, pmesh));
}

template <typename PolygonMesh>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh,
                                              std::size_t nb_components_to_keep)
{
  return keep_largest_connected_components(pmesh,
                                           nb_components_to_keep,
                                           internal::No_constraint<PolygonMesh>(),
                                           get(boost::vertex_index, pmesh),
                                           get(boost::face_index, pmesh));
}




} // namespace Polygon_mesh_processing

} // namespace CGAL

#endif //CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H
