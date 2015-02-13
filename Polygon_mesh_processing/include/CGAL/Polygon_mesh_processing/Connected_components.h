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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H
#define CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H

#include<set>
#include<vector>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/tuple.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/internal/corefinement/Polyhedron_constness_types.h>

namespace CGAL {
 namespace internal{
   namespace corefinement{

template <class Polyhedron>
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
  template <class T>
  bool operator()(T) const  {return true;}
};

template <class Polyhedron,class HDS=typename Polyhedron::HalfedgeDS>
class Build_polyhedron_subset : public ::CGAL::Modifier_base<HDS> {
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;       
  
  typedef typename HDS::Vertex::Point Point;
  std::list<Vertex_const_handle> points;
  std::list< std::vector<unsigned int> > facets;
    
  template <class Facet_iterator>
  typename Polyhedron::Halfedge_const_handle get_facet_halfedge(Facet_iterator facet_it) const
  {
    return (*facet_it)->halfedge();
  }

  typename Polyhedron::Halfedge_const_handle get_facet_halfedge(typename Polyhedron::Facet_const_handle facet) const
  {
    return facet->halfedge();
  }
  
public:
  template <class Facets_const_iterator>
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



template <class Polyhedron,class Adjacency_criterium,class Face_to_UF_handle_map,class Result>
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

template <class Polyhedron,class Adjacency_criterium,class Output_iterator>
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

template <class Polyhedron, class Adjacency_criterium, class Face_marker>
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

template <class Polyhedron, class Adjacency_criterium, class Face_marker, class OutputIterator>
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

template <class Polyhedron,class Output_iterator>
void extract_connected_components(const Polyhedron& P,Output_iterator out)
{
  extract_connected_components(P,Dummy_true(),out);
}

template <class Polyhedron, class Polyhedron_facet_index_map>
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
template <class Polyhedron, class Adjacency_criterium, class Polyhedron_facet_index_map>
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
  /*!
   * \ingroup PkgPolygonMeshProcessing
   *  Erases the small connected components and the isolated vertices.
   *  Keep `nb_components_to_keep` largest connected components.
   *  \return the number of connected components erased (ignoring isolated vertices).
   * \todo BGLize me
  */
  template <class PolygonMesh>
  std::size_t keep_largest_connected_components(PolygonMesh& pmesh, std::size_t nb_components_to_keep)
  {
    return pmesh.keep_largest_connected_components(nb_components_to_keep);
  }

  /*!
   * \ingroup PkgPolygonMeshProcessing
   *  Discovers all the faces in the same connected component as `seed_face` and puts them in `out`.
   * `seed_face` will also be added in `out`.
   *  Two faces are considered to be in the same connected component if they share an edge.
   *  \tparam PolygonMesh a model of `FaceGraph`
   *  \tparam EdgeConstraintMap a property map with the edge descriptor as key type and `bool` as vaue type
   *  \tparam OutputIterator an output iterator that accepts face descriptors.
   *  \returns the output iterator.
  */
  template <class PolygonMesh, class EdgeConstraintMap, class OutputIterator>
  OutputIterator
  connected_component(
    typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face,
    PolygonMesh& pmesh,
    EdgeConstraintMap ecmap,
    OutputIterator out)
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
        if(! ecmap[hd]){
            face_descriptor neighbor = face( opposite(hd, pmesh), pmesh );
            if ( neighbor != boost::graph_traits<PolygonMesh>::null_face() )
              stack.push_back(neighbor);
          }
      }
    }
    return out;
  }

namespace internal {

template <typename G>
struct No_constraint {
  No_constraint() { }

  No_constraint(G & g)
    : g(&g) { }

  template <typename T>
  bool operator[](const T & ) const {
    return true;
  }

  G* g;
};

  template <typename G, typename EdgeConstraintMap = No_constraint<G> >
struct No_border {
  No_border() 
  {}
  
    No_border(G & g, EdgeConstraintMap ecm = EdgeConstraintMap())
      : g(&g), ecm(ecm)
  { }
  
  bool operator()(typename boost::graph_traits<G>::edge_descriptor e) const {
    if(! is_border(e,*g)){
      return ! ecm[e];
    }
    return false;
  }

  G* g;
  EdgeConstraintMap ecm;
};



}// namespace internal
 

  /*!
   * \ingroup PkgPolygonMeshProcessing
   *  Discovers all the faces in the same connected component as `seed_face` and puts them in `out`.
   * `seed_face` will also be added in `out`.
   *  Two faces are considered to be in the same connected component if they share an edge.
   *  \tparam PolygonMesh a model of `FaceGraph`
   *  \tparam EdgeConstraintMap a property map with the edge descriptor as key type and `bool` as value type
   *  \tparam OutputIterator an output iterator that accepts face descriptors.
   *  \returns the output iterator.
  */
  template <class PolygonMesh, class OutputIterator>
  OutputIterator
  connected_component(
    typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face,
    PolygonMesh& pmesh,
    OutputIterator out)
  {
    return connected_component(seed_face, pmesh, internal::No_constraint<PolygonMesh>(pmesh), out);
  }


  /*!
   * \ingroup PkgPolygonMeshProcessing
   *  computes for each face the index of the connected components to which it belongs.
   *  Two faces are considered to be in the same connected component if they share an edge that is not constrained.
   *  \tparam PolygonMesh a model of `FaceGraph`
   *  \tparam EdgeConstraintMap a property map with the edge descriptor as key type and `bool` as value type.
   *  \tparam FaceIndexMap the property map with the face index as value type, and the index of its connected component
   *  \returns the number of connected components.
  */
  template <class PolygonMesh, class EdgeConstraintMap, class FaceIndexMap>
  typename boost::graph_traits<PolygonMesh>::faces_size_type
  connected_components(PolygonMesh& pmesh,
                       EdgeConstraintMap ecmap,
                       FaceIndexMap& fim)
  {
    typedef Dual<PolygonMesh> Dual;
    typedef boost::filtered_graph<Dual, internal::No_border<PolygonMesh,EdgeConstraintMap> > FiniteDual;

    Dual dual(pmesh);
    FiniteDual finite_dual(dual,internal::No_border<PolygonMesh,EdgeConstraintMap>(pmesh,ecmap));
    return boost::connected_components(finite_dual, fim);
  }

 /*!
   * \ingroup PkgPolygonMeshProcessing
   *  computes for each face the index of the connected components to which it belongs.
   *  Two faces are considered to be in the same connected component if they share an edge.
   *  \tparam PolygonMesh a model of `FaceGraph`
   *  \tparam FaceIndexMap the property map with the face index as value type, and the index of its connected component
   *  \returns the number of connected components.
  */
  template <class PolygonMesh, class FaceIndexMap>
  typename boost::graph_traits<PolygonMesh>::faces_size_type
  connected_components(PolygonMesh& pmesh,
                       FaceIndexMap& fim)
  {
    typedef Dual<PolygonMesh> Dual;
    typedef boost::filtered_graph<Dual, internal::No_border<PolygonMesh> > FiniteDual;

    Dual dual(pmesh);
    FiniteDual finite_dual(dual,internal::No_border<PolygonMesh>(pmesh));
    return boost::connected_components(finite_dual, fim);
  }
} // namespace Polygon_mesh_processing

} // namespace CGAL

#endif //CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H
