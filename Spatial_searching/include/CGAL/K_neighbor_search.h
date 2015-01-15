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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef  CGAL_K_NEIGHBOR_SEARCH_H
#define  CGAL_K_NEIGHBOR_SEARCH_H

#include <CGAL/internal/K_neighbor_search.h>
#include <CGAL/internal/Get_dimension_tag.h>




namespace CGAL {
	template <class SearchTraits, class Distance,class Splitter,class Tree>
class K_neighbor_search;

template <class SearchTraits, 
          class Distance= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
          class Splitter= Sliding_midpoint<SearchTraits> ,
          class Tree= Kd_tree<SearchTraits, Splitter, Tag_true> >
class K_neighbor_search: public internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> {
  typedef  internal::K_neighbor_search<SearchTraits,Distance,Splitter,Tree> Base;
  
public:
  typedef typename Base::FT FT;  
  typedef typename internal::Get_dimension_tag<SearchTraits>::Dimension D;

  K_neighbor_search(const Tree& tree, const typename Base::Query_item& q,  
    unsigned int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance(),bool sorted=true)
    : Base(q,k,Eps,Search_nearest,d) 
  {
    if (tree.empty()) return;
    compute_neighbors_general(tree.root(),tree.bounding_box());
    if (sorted) this->queue.sort();    
  };

private:  
  typedef typename Base::Node_const_handle Node_const_handle; 
  using Base::branch;

  void 
  compute_neighbors_general(typename Base::Node_const_handle N, const Kd_tree_rectangle<FT,D>& r) 
  {
    if (!(N->is_leaf())) {
       typename Tree::Internal_node_const_handle node =
        static_cast<typename Tree::Internal_node_const_handle>(N);
      this->number_of_internal_nodes_visited++;
      int new_cut_dim=node->cutting_dimension();
      FT  new_cut_val=node->cutting_value();

      Kd_tree_rectangle<FT,D> r_lower(r);

      // modifies also r_lower to lower half
      Kd_tree_rectangle<FT,D> r_upper(r_lower);
      r_lower.split(r_upper, new_cut_dim, new_cut_val);

      FT distance_to_lower_half;
      FT distance_to_upper_half;

      if (this->search_nearest) { 

        distance_to_lower_half = 
          this->distance_instance. min_distance_to_rectangle(this->query_object, 
                                                       r_lower);
                              
        distance_to_upper_half = 
          this->distance_instance.min_distance_to_rectangle(this->query_object, 
                                                      r_upper);
                      

      } 
      else
        { 

          distance_to_lower_half = 
            this->distance_instance.max_distance_to_rectangle(this->query_object, 
                                                        r_lower);

          distance_to_upper_half = 
            this->distance_instance.max_distance_to_rectangle(this->query_object, 
                                                        r_upper);

        }

      if ( (( this->search_nearest) && 
            (distance_to_lower_half < distance_to_upper_half)) 
           ||
           ((!this->search_nearest) && 
            (distance_to_lower_half >= 
             distance_to_upper_half))  )
        {
          compute_neighbors_general(node->lower(), r_lower);
          if (branch(distance_to_upper_half)) 
            compute_neighbors_general (node->upper(), r_upper);
        }  
      else
        {	compute_neighbors_general(node->upper(), r_upper);
        if (branch(distance_to_lower_half)) 
          compute_neighbors_general (node->lower(), 
                                     r_lower);
        }

    }
    else
      {
        // n is a leaf
         typename Tree::Leaf_node_const_handle node =
          static_cast<typename Tree::Leaf_node_const_handle>(N);
        this->number_of_leaf_nodes_visited++;
        if (node->size() > 0)
          for (typename Tree::iterator it = node->begin(); it != node->end(); it++) {
            this->number_of_items_visited++;
            FT distance_to_query_object =
              this->distance_instance.transformed_distance(this->query_object,*it);
            this->queue.insert(std::make_pair(&(*it),distance_to_query_object));
          }
      }
  }
  
}; // class 

} // namespace CGAL


#endif  // CGAL_K_NEIGHBOR_SEARCH_H
