// Copyright (c) 2002 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef  CGAL_K_NEIGHBOR_SEARCH_H
#define  CGAL_K_NEIGHBOR_SEARCH_H
#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Splitters.h>

namespace CGAL {

  template <class SearchTraits,
            class Distance_=Euclidean_distance<SearchTraits>,
            class Splitter_=Sliding_midpoint<SearchTraits> , 
            class Tree_=Kd_tree<SearchTraits, Splitter_, Tag_false> >

  class K_neighbor_search {

  public:

    typedef Splitter_ Splitter;
    typedef Distance_ Distance;
    typedef Tree_ Tree;
    typedef typename SearchTraits::Point_d Point_d;
    typedef typename SearchTraits::FT FT;
    typedef std::pair<Point_d,FT> Point_with_transformed_distance;
  
    typedef typename Tree::Node_handle Node_handle;
  
    typedef typename Tree::Point_d_iterator Point_d_iterator;
    typedef Kd_tree_rectangle<SearchTraits> Rectangle;
    typedef typename Distance::Query_item Query_item;

  private:

    int number_of_internal_nodes_visited;
    int number_of_leaf_nodes_visited;
    int number_of_items_visited;
  
    bool search_nearest;
  
    FT multiplication_factor;
    Query_item query_object;
    int total_item_number;
    FT distance_to_root;   
  
    typedef std::list<Point_with_transformed_distance> NN_list;

  public:

    typedef typename NN_list::const_iterator iterator;

  private:

    NN_list l;
    int max_k;
    int actual_k;
  
    Distance distance_instance;
  
    bool 
    branch(FT distance) 
    {
      if (actual_k<max_k) return true;
      else 
	if (search_nearest) return 
			      ( distance * multiplication_factor < l.rbegin()->second);
	else return ( distance >
		      l.begin()->second * multiplication_factor);
    
    }
  
    void 
    insert(Point_d* I, FT dist) 
    {
      bool insert;
      if (actual_k<max_k) insert=true;
      else 
	if (search_nearest) insert=
			      ( dist < l.rbegin()->second ); 
	else insert=(dist > l.begin()->second);
      if (insert) {
	actual_k++;	 	
	typename NN_list::iterator it=l.begin();
	for (; (it != l.end()); ++it) 
	  { if (dist < it->second) break;}
	Point_with_transformed_distance NN_Candidate(*I,dist);
	l.insert(it,NN_Candidate);
	if (actual_k > max_k) {
	  actual_k--;
	  if (search_nearest) l.pop_back();
	  else l.pop_front();
	}
      }
    
    }
  
	
  public:


    iterator 
    begin() const
    {
      return l.begin();
    }

    iterator 
    end() const
    {
      return l.end();
    }

    // constructor
    K_neighbor_search(Tree& tree, const Query_item& q, 
		      int k=1, FT Eps=FT(0.0), 
		      bool Search_nearest=true,
		      const Distance& d=Distance())
      : distance_instance(d), max_k(k), actual_k(0), search_nearest(Search_nearest), query_object(q),
	total_item_number(tree.size()), number_of_leaf_nodes_visited(0), number_of_internal_nodes_visited(0),
	number_of_items_visited(0), multiplication_factor(distance_instance.transformed_distance(FT(1.0)+Eps))
    {
      compute_neighbors_general(tree.root(), tree.bounding_box());
    }

    // Print statistics of the general standard search process.
    std::ostream& statistics (std::ostream& s) {
      s << "General search statistics:" << std::endl;
      s << "Number of internal nodes visited:" 
	<< number_of_internal_nodes_visited << std::endl;
      s << "Number of leaf nodes visited:" 
	<< number_of_leaf_nodes_visited << std::endl;
      s << "Number of items visited:" 
	<< number_of_items_visited << std::endl;
      return s;
    }


  private:
   

    void 
    compute_neighbors_general(Node_handle N, const Kd_tree_rectangle<SearchTraits>& r) 
    {
      if (!(N->is_leaf())) {
	number_of_internal_nodes_visited++;
	int new_cut_dim=N->cutting_dimension();
	FT  new_cut_val=N->cutting_value();

	Kd_tree_rectangle<SearchTraits> r_lower(r);

	// modifies also r_lower to lower half
	Kd_tree_rectangle<SearchTraits> r_upper(r_lower);
	r_lower.split(r_upper, new_cut_dim, new_cut_val);

	FT distance_to_lower_half;
	FT distance_to_upper_half;

	if (search_nearest) { 

	  distance_to_lower_half = 
	    distance_instance. min_distance_to_rectangle(query_object, 
							 r_lower);
				
	  distance_to_upper_half = 
	    distance_instance.min_distance_to_rectangle(query_object, 
							r_upper);
			

	} 
	else
	  { 

	    distance_to_lower_half = 
	      distance_instance.max_distance_to_rectangle(query_object, 
							  r_lower);

	    distance_to_upper_half = 
	      distance_instance.max_distance_to_rectangle(query_object, 
							  r_upper);

	  }
 
	if ( (( search_nearest) && 
	      (distance_to_lower_half < distance_to_upper_half)) 
	     ||
	     ((!search_nearest) && 
	      (distance_to_lower_half >= 
	       distance_to_upper_half))  )
	  {
	    compute_neighbors_general(N->lower(), r_lower);
	    if (branch(distance_to_upper_half)) 
	      compute_neighbors_general (N->upper(), r_upper);
	  }  
	else
	  {	compute_neighbors_general(N->upper(), r_upper);
	  if (branch(distance_to_lower_half)) 
	    compute_neighbors_general (N->lower(), 
				       r_lower);
	  }

      }
      else
	{
	  // n is a leaf
	  number_of_leaf_nodes_visited++;
	  if (N->size() > 0)
	    for (Point_d_iterator it = N->begin(); it != N->end(); it++) {
	      number_of_items_visited++;
	      FT distance_to_query_object =
		distance_instance.transformed_distance(query_object,**it);
	      insert(*it,distance_to_query_object);
	    }
	}
    }
        

    
  }; // class 



} // namespace CGAL


#endif  // CGAL_K_NEIGHBOR_SEARCH_H
