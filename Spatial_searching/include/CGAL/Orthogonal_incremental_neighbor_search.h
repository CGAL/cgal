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

#ifndef CGAL_ORTHOGONAL_INCREMENTAL_NEIGHBOR_SEARCH
#define CGAL_ORTHOGONAL_INCREMENTAL_NEIGHBOR_SEARCH

#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>

namespace CGAL {

  template <class SearchTraits, 
            class Distance_= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
            class Splitter_ = Sliding_midpoint<SearchTraits>,
            class Tree_= Kd_tree<SearchTraits, Splitter_, Tag_true> >
  class Orthogonal_incremental_neighbor_search {

  public:
    typedef Splitter_ Splitter;
    typedef Tree_  Tree;
    typedef Distance_ Distance;
    typedef typename SearchTraits::Point_d Point_d;
    typedef typename Distance::Query_item Query_item;
    typedef typename SearchTraits::FT FT;
    typedef typename Tree::Point_d_iterator Point_d_iterator;
    typedef typename Tree::Node_const_handle Node_const_handle;

    typedef std::pair<Point_d,FT> Point_with_transformed_distance;
    typedef std::pair<Node_const_handle,FT> Node_with_distance;
    typedef std::vector<Node_with_distance*> Node_with_distance_vector;
    typedef std::vector<Point_with_transformed_distance*> Point_with_transformed_distance_vector;

    template<class T>
    struct Object_wrapper
    {   
      T object;
      Object_wrapper(const T& t):object(t){}
      const T& operator* () const { return object; }
      const T* operator-> () const { return &object; }
    };

    class Iterator_implementation {
      SearchTraits traits;
    public:

      int number_of_neighbours_computed;
      int number_of_internal_nodes_visited;
      int number_of_leaf_nodes_visited;
      int number_of_items_visited;

    private:

      Distance Orthogonal_distance_instance;
    
      FT multiplication_factor;

      Query_item query_point;

      FT distance_to_root;

      bool search_nearest_neighbour;

      FT rd;

      class Priority_higher {
      public:

        bool search_nearest;

        Priority_higher(bool search_the_nearest_neighbour) 
	  : search_nearest(search_the_nearest_neighbour)
	{} 

        //highest priority is smallest distance
        bool 
	operator() (Node_with_distance* n1, Node_with_distance* n2) const
	{
	  return (search_nearest) ? (n1->second > n2->second) : (n2->second > n1->second);
        }
      };

      class Distance_smaller {

      public:

        bool search_nearest;

        Distance_smaller(bool search_the_nearest_neighbour) 
	  : search_nearest(search_the_nearest_neighbour)
	{}

        //highest priority is smallest distance
        bool operator() (Point_with_transformed_distance* p1, Point_with_transformed_distance* p2) const
	{
	  return (search_nearest) ? (p1->second > p2->second) : (p2->second > p1->second);
        }
      };


      std::priority_queue<Node_with_distance*, Node_with_distance_vector,
	                  Priority_higher> PriorityQueue;

    public:
      std::priority_queue<Point_with_transformed_distance*, Point_with_transformed_distance_vector,
	                  Distance_smaller> Item_PriorityQueue;


    public:

      int reference_count;

    

      // constructor
      Iterator_implementation(const Tree& tree,const Query_item& q, const Distance& tr,
			      FT Eps=FT(0.0), bool search_nearest=true)
	: traits(tree.traits()),number_of_neighbours_computed(0), number_of_internal_nodes_visited(0), 
	number_of_leaf_nodes_visited(0), number_of_items_visited(0),
	Orthogonal_distance_instance(tr), multiplication_factor(Orthogonal_distance_instance.transformed_distance(FT(1.0)+Eps)), 
	query_point(q), search_nearest_neighbour(search_nearest), 
	PriorityQueue(Priority_higher(search_nearest)), Item_PriorityQueue(Distance_smaller(search_nearest)),
	reference_count(1)
	  
	  
      {
        if (tree.empty()) return;
        
        // if (search_nearest) 
	distance_to_root=
	  Orthogonal_distance_instance.min_distance_to_rectangle(q,
								  tree.bounding_box());
        // else distance_to_root=
   	// Orthogonal_Distance_instance->max_distance_to_rectangle(q,
	//					tree.bounding_box());

        Node_with_distance *The_Root = new Node_with_distance(tree.root(),
							      distance_to_root);
        PriorityQueue.push(The_Root);

        // rd is the distance of the top of the priority queue to q
        rd=The_Root->second;
        Compute_the_next_nearest_neighbour();
      }

      // * operator
      const Point_with_transformed_distance& 
      operator* () const 
      {
	return *(Item_PriorityQueue.top());
      }

      // prefix operator
      Iterator_implementation& 
      operator++() 
      {
        Delete_the_current_item_top();
        Compute_the_next_nearest_neighbour();
        return *this;
      }

      // postfix operator
      Object_wrapper<Point_with_transformed_distance>
      operator++(int) 
      {
        Object_wrapper<Point_with_transformed_distance> result( *(Item_PriorityQueue.top()) );
        ++*this;
        return result;
      }

      // Print statistics of the general priority search process.
      std::ostream& 
      statistics (std::ostream& s) const {
    	s << "Orthogonal priority search statistics:" 
	  << std::endl;
    	s << "Number of internal nodes visited:" 
	  << number_of_internal_nodes_visited << std::endl;
    	s << "Number of leaf nodes visited:" 
	  << number_of_leaf_nodes_visited << std::endl;
    	s << "Number of items visited:" 
	  << number_of_items_visited << std::endl;
        s << "Number of neighbours computed:" 
	  << number_of_neighbours_computed << std::endl;
        return s;
      }


      //destructor
      ~Iterator_implementation() 
      {
	while (!PriorityQueue.empty()) {
	  Node_with_distance* The_top=PriorityQueue.top();
	  PriorityQueue.pop();
	  delete The_top;
	}
	while (!Item_PriorityQueue.empty()) {
	  Point_with_transformed_distance* The_top=Item_PriorityQueue.top();
	  Item_PriorityQueue.pop();
	  delete The_top;
        }
      }

    private:

      void 
      Delete_the_current_item_top() 
      {
        Point_with_transformed_distance* The_item_top=Item_PriorityQueue.top();
        Item_PriorityQueue.pop();
        delete The_item_top;
      }

      void 
      Compute_the_next_nearest_neighbour() 
      {
        // compute the next item
        bool next_neighbour_found=false;
        if (!(Item_PriorityQueue.empty())) {
	  if (search_nearest_neighbour)
	    next_neighbour_found=
	      (multiplication_factor*rd > Item_PriorityQueue.top()->second);
	  else
	    next_neighbour_found=
	      (rd < multiplication_factor*Item_PriorityQueue.top()->second);
        }
	typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
	typename SearchTraits::Cartesian_const_iterator_d query_point_it = construct_it(query_point);
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue.empty())) {
	  Node_with_distance* The_node_top=PriorityQueue.top();
	  Node_const_handle N= The_node_top->first;
	  PriorityQueue.pop();
	  delete The_node_top;
	  FT copy_rd=rd;
	  while (!(N->is_leaf())) { // compute new distance
	    number_of_internal_nodes_visited++;
	    int new_cut_dim=N->cutting_dimension();
	    FT old_off, new_rd;
	    FT new_off =
	      *(query_point_it + new_cut_dim) -
	      N->cutting_value();
	    if (new_off < FT(0.0)) {
	      old_off=
		*(query_point_it + new_cut_dim)-N->low_value();
	      if (old_off>FT(0.0)) old_off=FT(0.0);
	      new_rd=
		Orthogonal_distance_instance.new_distance(copy_rd,old_off,new_off,new_cut_dim);
	      CGAL_assertion(new_rd >= copy_rd);
	      if (search_nearest_neighbour) {
		Node_with_distance *Upper_Child =
		  new Node_with_distance(N->upper(), new_rd);
		PriorityQueue.push(Upper_Child);
		N=N->lower();
	      }
	      else {
		Node_with_distance *Lower_Child =
		  new Node_with_distance(N->lower(), copy_rd);
		PriorityQueue.push(Lower_Child);
		N=N->upper();
		copy_rd=new_rd;
	      }

	    }
	    else { // compute new distance
	      old_off= N->high_value() -
		*(query_point_it+new_cut_dim);
	      if (old_off>FT(0.0)) old_off=FT(0.0);
	      new_rd=Orthogonal_distance_instance.new_distance(copy_rd,old_off,new_off,new_cut_dim);  
	      CGAL_assertion(new_rd >= copy_rd);
	      if (search_nearest_neighbour) {
		Node_with_distance *Lower_Child =
		  new Node_with_distance(N->lower(), new_rd);
		PriorityQueue.push(Lower_Child);
		N=N->upper();
	      }
	      else {
		Node_with_distance *Upper_Child =
		  new Node_with_distance(N->upper(), copy_rd);
		PriorityQueue.push(Upper_Child);
		N=N->lower();
		copy_rd=new_rd;
	      }
	    }
	  }
	  // n is a leaf
	  number_of_leaf_nodes_visited++;
	  if (N->size() > 0) {
	    for (Point_d_iterator it=N->begin(); it != N->end(); it++) {
	      number_of_items_visited++;
	      FT distance_to_query_point=
		Orthogonal_distance_instance.transformed_distance(query_point,**it);
	      Point_with_transformed_distance *NN_Candidate=
		new Point_with_transformed_distance(**it,distance_to_query_point);
	      Item_PriorityQueue.push(NN_Candidate);
	    }
	    // old top of PriorityQueue has been processed,
	    // hence update rd
                
	    if (!(PriorityQueue.empty()))  {
	      rd = PriorityQueue.top()->second;
	      if (search_nearest_neighbour)
		next_neighbour_found =
                  (multiplication_factor*rd > 
		   Item_PriorityQueue.top()->second);
	      else
		next_neighbour_found =
                  (multiplication_factor*rd < 
		   Item_PriorityQueue.top()->second);
	    }
	    else // priority queue empty => last neighbour found
	      {
		next_neighbour_found=true;
	      }

	    number_of_neighbours_computed++;
	  }
        }   // next_neighbour_found or priority queue is empty
        // in the latter case also the item priority quee is empty
      }
    }; // class Iterator_implementaion
  





    class iterator;

    typedef std::vector<FT> Distance_vector;

  public:

    // constructor
    Orthogonal_incremental_neighbor_search(const Tree& tree,  
					   const Query_item& q, FT Eps = FT(0.0), 
					   bool search_nearest=true, const Distance& tr=Distance())
      : m_tree(tree),m_query(q),m_dist(tr),m_Eps(Eps),m_search_nearest(search_nearest)
    {}

    iterator 
    begin() 
    {
      return iterator(m_tree,m_query,m_dist,m_Eps,m_search_nearest);
    }

    iterator 
    end() 
    {
      return iterator();
    }

    std::ostream& 
    statistics(std::ostream& s) 
    {
      begin()->statistics(s);
      return s;
    }




    class iterator {

    public:

      typedef std::input_iterator_tag iterator_category;
      typedef Point_with_transformed_distance       value_type;
      typedef Point_with_transformed_distance*      pointer;
      typedef const Point_with_transformed_distance&      reference;
      typedef std::size_t               size_type;
      typedef std::ptrdiff_t            difference_type;
      typedef int distance_type;

      //class Iterator_implementation;
      Iterator_implementation *Ptr_implementation;


    public:

      // default constructor
      iterator() 
      : Ptr_implementation(0)
      {}

      int 
      the_number_of_items_visited() 
      {
        return Ptr_implementation->number_of_items_visited;
      }

      // constructor
      iterator(const Tree& tree,const Query_item& q, const Distance& tr=Distance(), FT eps=FT(0.0), 
	       bool search_nearest=true)
	: Ptr_implementation(new Iterator_implementation(tree, q, tr, eps, search_nearest))
	{}

      // copy constructor
      iterator(const iterator& Iter) 
      {
        Ptr_implementation = Iter.Ptr_implementation;
        if (Ptr_implementation != 0) Ptr_implementation->reference_count++;
      }

      iterator& operator=(const iterator& Iter)
      {
        if (Ptr_implementation != Iter.Ptr_implementation){
          if (Ptr_implementation != 0 && --(Ptr_implementation->reference_count)==0) {
              delete Ptr_implementation;
          }
          Ptr_implementation = Iter.Ptr_implementation;
          if (Ptr_implementation != 0) Ptr_implementation->reference_count++;
        }
        return *this;
      }      
      
      
      const Point_with_transformed_distance& 
      operator* () const 
      {
	return *(*Ptr_implementation);
      }
      
      // -> operator
      const Point_with_transformed_distance*
      operator-> () const 
      {
	return &*(*Ptr_implementation);
      }

      // prefix operator
      iterator& 
      operator++() 
      {
        ++(*Ptr_implementation);
        return *this;
      }

      // postfix operator
      Object_wrapper<Point_with_transformed_distance>
      operator++(int) 
      {
	return (*Ptr_implementation)++;
      }


      bool 
      operator==(const iterator& It) const 
      {
        if (
	    ((Ptr_implementation == 0) || 
	     Ptr_implementation->Item_PriorityQueue.empty()) &&
	    ((It.Ptr_implementation == 0) ||  
	     It.Ptr_implementation->Item_PriorityQueue.empty())
	    )
	  return true;
        // else
        return (Ptr_implementation == It.Ptr_implementation);
      }

      bool 
      operator!=(const iterator& It) const 
      {
        return !(*this == It);
      }

      std::ostream& 
      statistics (std::ostream& s) 
      {
    	Ptr_implementation->statistics(s);
        return s;
      }

      ~iterator() 
      {
        if (Ptr_implementation != 0) {
	  Ptr_implementation->reference_count--;
	  if (Ptr_implementation->reference_count==0) {
	    delete Ptr_implementation;
	    Ptr_implementation = 0;
	  }
        }
      }


    }; // class iterator

    //data members
    const Tree& m_tree;
    Query_item m_query;
    Distance m_dist;
    FT m_Eps; 
    bool m_search_nearest;
  }; // class 

  template <class Traits, class Query_item, class Distance>
  void swap (typename Orthogonal_incremental_neighbor_search<Traits, 
	     Query_item, Distance>::iterator& x,
	     typename Orthogonal_incremental_neighbor_search<Traits, 
	     Query_item, Distance>::iterator& y) 
  {
    typename Orthogonal_incremental_neighbor_search<Traits, 
      Query_item, Distance>::iterator::Iterator_implementation
      *tmp = x.Ptr_implementation;
    x.Ptr_implementation  = y.Ptr_implementation;
    y.Ptr_implementation = tmp;
  }

} // namespace CGAL

#endif // CGAL_ORTHOGONAL_INCREMENTAL_NEIGHBOR_SEARCH_H
