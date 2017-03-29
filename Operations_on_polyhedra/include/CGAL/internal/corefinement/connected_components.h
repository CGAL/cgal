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

#include <CGAL/license/Polygon_mesh_processing.h>


#include<set>
#include<vector>

#include <boost/graph/connected_components.hpp>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/assertions.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Union_find.h>
#include <CGAL/internal/corefinement/Polyhedron_constness_types.h>

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


} // namespace CGAL

#endif //CGAL_INTERNAL_POLYHEDRON_SUBSET_EXTRACTION_H
