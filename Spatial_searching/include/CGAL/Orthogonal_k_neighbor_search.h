// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
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
// Author(s)     : Gael Guennebaud (gael.guennebaud@inria.fr), Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#define CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H

#include <CGAL/internal/K_neighbor_search.h>

namespace CGAL {

template <class SearchTraits, 
          class Distance= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
          class Splitter= Sliding_midpoint<SearchTraits> ,
          class Tree= Kd_tree<SearchTraits, Splitter, Tag_true> >
class Orthogonal_k_neighbor_search: public internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> {
  typedef  internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> Base;

  typename SearchTraits::Cartesian_const_iterator_d query_object_it;
  
public:
  typedef typename Base::FT FT;

  Orthogonal_k_neighbor_search(const Tree& tree, const typename Base::Query_item& q,  
                               unsigned int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance(),bool sorted=true)
    : Base(q,k,Eps,Search_nearest,d) 
  {
    if (tree.empty()) return;

    FT distance_to_root;
    if (this->search_nearest) 
      distance_to_root = this->distance_instance.min_distance_to_rectangle(q, tree.bounding_box());
    else 
      distance_to_root = this->distance_instance.max_distance_to_rectangle(q, tree.bounding_box());


    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=tree.traits().construct_cartesian_const_iterator_d_object();
    query_object_it = construct_it(this->query_object);
      
    compute_neighbors_orthogonally(tree.root(), distance_to_root);      
      
    if (sorted) this->queue.sort();
  }
private:

  void compute_neighbors_orthogonally(typename Base::Node_const_handle N, FT rd)
  {
    if (!(N->is_leaf())) 
    {
      this->number_of_internal_nodes_visited++;
      int new_cut_dim=N->cutting_dimension();
      FT old_off, new_rd;
      FT new_off = *(query_object_it + new_cut_dim) - N->cutting_value();
      if ( ((new_off <  FT(0.0)) && (this->search_nearest)) ||
           ((new_off >= FT(0.0)) && (!this->search_nearest)) ) 
      {
        compute_neighbors_orthogonally(N->lower(),rd);
        if (this->search_nearest) {
          old_off= *(query_object_it + new_cut_dim) - N->low_value();
          if (old_off>FT(0.0)) 
            old_off=FT(0.0);
        }
        else {	
          old_off= *(query_object_it + new_cut_dim) - N->high_value();
          if (old_off<FT(0.0)) 
            old_off=FT(0.0);
        }
        new_rd = this->distance_instance.new_distance(rd,old_off,new_off,new_cut_dim);
        if (this->branch(new_rd)) 
          compute_neighbors_orthogonally(N->upper(), new_rd);                               
      }
      else // compute new distance
      {
        compute_neighbors_orthogonally(N->upper(),rd); 
        if (this->search_nearest) {
          old_off= N->high_value() - *(query_object_it + new_cut_dim);
          if (old_off>FT(0.0)) 
            old_off=FT(0.0);
        }
        else  {       
          old_off= N->low_value() - *(query_object_it + new_cut_dim);
          if (old_off<FT(0.0)) 
            old_off=FT(0.0);
        }  
        new_rd = this->distance_instance. new_distance(rd,old_off,new_off,new_cut_dim);
        if (this->branch(new_rd)) 
          compute_neighbors_orthogonally(N->lower(), new_rd);       
      }
    }
    else
    {
      // n is a leaf
      this->number_of_leaf_nodes_visited++;
      if (N->size() > 0)
      {
        for (typename Base::Point_d_iterator it=N->begin(); it != N->end(); it++) 
        {
          this->number_of_items_visited++;
          FT distance_to_query_object=
            this->distance_instance.transformed_distance(this->query_object,**it);
          this->queue.insert(std::make_pair(*it,distance_to_query_object));
        }
      }
    }
  }    
  
}; // class 

} // namespace CGAL

#endif  // CGAL_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
