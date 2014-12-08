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

#ifndef CGAL_INCREMENTAL_NEIGHBOR_SEARCH_H
#define CGAL_INCREMENTAL_NEIGHBOR_SEARCH_H

#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Euclidean_distance.h>

namespace CGAL {

  template <class SearchTraits, 
            class Distance_=typename internal::Spatial_searching_default_distance<SearchTraits>::type,
            class Splitter_ = Sliding_midpoint<SearchTraits>,
            class Tree_=Kd_tree<SearchTraits, Splitter_, Tag_false> >
  class Incremental_neighbor_search { 

  public:

    typedef Distance_ Distance;
    typedef Tree_     Tree;
    typedef typename SearchTraits::Point_d Point_d;
    typedef typename SearchTraits::FT FT;
    typedef typename SearchTraits::Dimension Dimension;
    typedef typename Tree::Point_d_iterator Point_d_iterator;
    typedef typename Tree::Node_const_handle Node_const_handle;
    typedef typename Tree::Splitter Splitter;
    typedef Kd_tree_rectangle<FT,Dimension> Node_box;
    typedef typename Distance::Query_item Query_item;

    class Cell {

    private:

      Node_box* the_box;
      Node_const_handle the_node;

    public:

      // constructor
      Cell (Node_box* Nb, Node_const_handle N)
      :the_box(Nb), the_node(N)
      {}

      Node_box* 
      box() 
      {
	return the_box;
      }

      Node_const_handle    
      node() 
      {
	return the_node;
      }
    };
    
    

    typedef std::pair<Point_d,FT> Point_with_transformed_distance;
    typedef std::pair<Cell*,FT> Cell_with_distance;

    class iterator;

    typedef std::vector<Cell_with_distance*> Cell_with_distance_vector;
    typedef std::vector<Point_with_transformed_distance*> Point_with_distance_vector;
    typedef std::vector<FT> Distance_vector;

    //data members
    const Tree& m_tree;
    Query_item m_query;
    Distance m_dist;
    FT m_Eps; 
    bool m_search_nearest;    

  public:

    // constructor
    Incremental_neighbor_search(const Tree& tree, const Query_item& q,
				FT Eps=FT(0.0), bool search_nearest=true, 
				const Distance& tr=Distance()): 
          m_tree(tree),m_query(q),m_dist(tr),m_Eps(Eps),m_search_nearest(search_nearest)
    {}

    iterator 
    begin() const
    {
      return iterator(m_tree,m_query,m_dist,m_Eps,m_search_nearest);
    }

    iterator 
    end() const
    {
      return iterator();
    }

    std::ostream&  
    statistics(std::ostream& s) 
    {
      begin()->statistics(s);
      return s;
    }


    template<class T>
    struct Object_wrapper
    {
      T object;
      Object_wrapper(const T& t):object(t){}
      const T& operator* () const { return object; }
      const T* operator-> () const { return &object; }
    };    


    class iterator {

    public:

      typedef std::input_iterator_tag iterator_category;
      typedef Point_with_transformed_distance       value_type;
      typedef Point_with_transformed_distance*      pointer;
      typedef const Point_with_transformed_distance&      reference;
      typedef std::size_t               size_type;
      typedef std::ptrdiff_t            difference_type;
      typedef int distance_type;
    
    
      class Iterator_implementation;
      Iterator_implementation *ptr;


      // default constructor
      iterator() 
	: ptr(0)
      {}

      int 
      the_number_of_items_visited() 
      {
        return ptr->number_of_items_visited;
      }

      // constructor
      iterator(const Tree& tree, const Query_item& q, const Distance& tr, FT eps, 
	       bool search_nearest)
        : ptr(new Iterator_implementation(tree, q, tr, eps, search_nearest))
      {}
      
      // copy constructor
      iterator(const iterator& Iter)
	: ptr(Iter.ptr)
      {
	if (ptr != 0) ptr->reference_count++;
      }
      
      iterator& operator=(const iterator& Iter)
      {
        if (ptr!=Iter.ptr){
          if (ptr != 0 && --(ptr->reference_count)==0) {
              delete ptr;
          }
          ptr = Iter.ptr;
          if (ptr != 0) ptr->reference_count++;
        }
        return *this;
      }      
      
      const Point_with_transformed_distance& 
      operator* () const
      {
	return *(*ptr);
      }

      // -> operator
      const Point_with_transformed_distance*
      operator-> () const 
      {
	return &*(*ptr);
      }      
      
      // prefix operator
      iterator& operator++() 
      {
        ++(*ptr);
        return *this;
      }

      // postfix operator
      Object_wrapper<Point_with_transformed_distance>
      operator++(int) 
      {
	return (*ptr)++;
      }

      bool 
      operator==(const iterator& It) const 
      {
        if ( ((ptr == 0) || 
	      ptr->Item_PriorityQueue.empty()) &&
	     ((It.ptr == 0) ||  
	      It.ptr->Item_PriorityQueue.empty())
	     )
	  return true;
        // else
        return (ptr == It.ptr);
      }

      bool 
      operator!=(const iterator& It) const
      {
        return !(*this == It);
      }

      std::ostream& 
      statistics (std::ostream& s) 
      {
    	ptr->statistics(s);
        return s;
      }

      ~iterator() 
      {
        if (ptr != 0) {
	  ptr->reference_count--;
	  if (ptr->reference_count==0) {
	    delete ptr;
	    ptr = 0;
	  }
        }
      }
      

      class Iterator_implementation {

      private:

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
	  bool operator() (Cell_with_distance* n1, Cell_with_distance* n2) const 
	  {
	    return (search_nearest)? (n1->second > n2->second) : (n2->second > n1->second);
	  }
	};


	class Distance_smaller {

	public:

	  bool search_nearest;

	  Distance_smaller(bool search_the_nearest_neighbour)
	    :search_nearest(search_the_nearest_neighbour)
	  {}

	  //highest priority is smallest distance
	  bool 
	  operator() (Point_with_transformed_distance* p1, Point_with_transformed_distance* p2) const 
	  {
	    return (search_nearest) ? (p1->second > p2->second) : (p2->second > p1->second);
	  }
	};


	std::priority_queue<Cell_with_distance*, Cell_with_distance_vector,
	                    Priority_higher> PriorityQueue;
      public:
	std::priority_queue<Point_with_transformed_distance*, Point_with_distance_vector,
	                    Distance_smaller> Item_PriorityQueue;
    

	Distance distance;

      public:

	int reference_count;

	int number_of_internal_nodes_visited;
	int number_of_leaf_nodes_visited;
	int number_of_items_visited;
	int number_of_neighbours_computed;

	// constructor
	Iterator_implementation(const Tree& tree, const Query_item& q,const Distance& tr,
				FT Eps, bool search_nearest)
	  : query_point(q), search_nearest_neighbour(search_nearest), 
	  PriorityQueue(Priority_higher(search_nearest)),
	  Item_PriorityQueue(Distance_smaller(search_nearest)),
	  distance(tr), reference_count(1), number_of_internal_nodes_visited(0), 
	  number_of_leaf_nodes_visited(0), number_of_items_visited(0),
	  number_of_neighbours_computed(0)
	{
          if (tree.empty()) return;
            
	  multiplication_factor= distance.transformed_distance(FT(1)+Eps);

	  Node_box *bounding_box = new Node_box((tree.bounding_box()));

	  if (search_nearest) distance_to_root=
				distance.min_distance_to_rectangle(q,*bounding_box);
	  else distance_to_root=
		 distance.max_distance_to_rectangle(q,*bounding_box);

        

	  Cell *Root_Cell = new Cell(bounding_box,tree.root());
	  Cell_with_distance  *The_Root = 
	    new Cell_with_distance(Root_Cell,distance_to_root);

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
	statistics (std::ostream& s) const
	{
	  s << "General priority search statistics:" << std::endl;
	  s << "Number of internal nodes visited:" << 
	    number_of_internal_nodes_visited << std::endl;
	  s << "Number of leaf nodes visited:" << 
	    number_of_leaf_nodes_visited << std::endl;
	  s << "Number of points visited:" << 
	    number_of_items_visited << std::endl;
	  s << "Number of neighbours computed:" << 
	    number_of_neighbours_computed << std::endl;
	  return s;
	}

	//destructor
	~Iterator_implementation() 
	{
	  while (! PriorityQueue.empty()) {
	    Cell_with_distance* The_top=PriorityQueue.top();
	    PriorityQueue.pop();
	    delete The_top->first->box();
	    delete The_top->first;
	    delete The_top;
	  }
	  while (! Item_PriorityQueue.empty()) {
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
	      next_neighbour_found =
		(multiplication_factor*rd > Item_PriorityQueue.top()->second);
	    else
	      next_neighbour_found =
		(rd < multiplication_factor*Item_PriorityQueue.top()->second);
	  }
	  while ((!next_neighbour_found) && (!PriorityQueue.empty())) {
                
	    Cell_with_distance* The_node_top = PriorityQueue.top();
	    Node_const_handle N = The_node_top->first->node();
	    Node_box* B = The_node_top->first->box();
	    PriorityQueue.pop();
	    delete The_node_top->first;
	    delete The_node_top;

	    while (!(N->is_leaf())) {
              typename Tree::Internal_node_const_handle node =
                static_cast<typename Tree::Internal_node_const_handle>(N);
	      number_of_internal_nodes_visited++;
	      int new_cut_dim = node->cutting_dimension();
	      FT  new_cut_val = node->cutting_value();
                        
	      Node_box* lower_box = new Node_box(*B);
	      Node_box* upper_box = new Node_box(*B); 
		lower_box->split(*upper_box,new_cut_dim, new_cut_val);
	      delete B;
	      if (search_nearest_neighbour) {
		FT distance_to_box_lower =
		  distance.min_distance_to_rectangle(query_point, *lower_box);
		FT distance_to_box_upper =
		  distance.min_distance_to_rectangle(query_point, *upper_box);
		if (distance_to_box_lower <= distance_to_box_upper) {

		  Cell* C_upper = new Cell(upper_box, node->upper());
		  Cell_with_distance *Upper_Child =
		    new Cell_with_distance(C_upper,distance_to_box_upper);
		  PriorityQueue.push(Upper_Child);
		  N=node->lower();
		  B=lower_box;
		} else {
		  Cell* C_lower = new Cell(lower_box, node->lower());
		  Cell_with_distance *Lower_Child =
		    new Cell_with_distance(C_lower,distance_to_box_lower);
		  PriorityQueue.push(Lower_Child);
		  N=node->upper();
		  B=upper_box;
		}
	      }
	      else { // search furthest
		FT distance_to_box_lower =
		  distance.max_distance_to_rectangle(query_point, *lower_box);
		FT distance_to_box_upper =
		  distance.max_distance_to_rectangle(query_point, *upper_box);
		if (distance_to_box_lower >= distance_to_box_upper) {
		  Cell* C_upper = new Cell(upper_box, node->upper());
		  Cell_with_distance *Upper_Child =
		    new Cell_with_distance(C_upper,distance_to_box_upper);
		  PriorityQueue.push(Upper_Child);
		  N=node->lower();
		  B=lower_box;
		}
		else {
		  Cell* C_lower = new Cell(lower_box, node->lower());
		  Cell_with_distance *Lower_Child =
		    new Cell_with_distance(C_lower,distance_to_box_lower);
		  PriorityQueue.push(Lower_Child);
		  N=node->upper();
		  B=upper_box;
		}
	      }
	    }
	    delete B;
             typename Tree::Leaf_node_const_handle node =
              static_cast<typename Tree::Leaf_node_const_handle>(N);
	    number_of_leaf_nodes_visited++;
	    if (node->size() > 0) {
	      for (typename Tree::iterator it = node->begin(); it != node->end(); it++) {
		number_of_items_visited++;
		FT distance_to_query_point=
		  distance.transformed_distance(query_point,*it);
		Point_with_transformed_distance *NN_Candidate=
		  new Point_with_transformed_distance(*it,distance_to_query_point);
		Item_PriorityQueue.push(NN_Candidate);
	      }
	      // old top of PriorityQueue has been processed,
	      // hence update rd
	      if (!PriorityQueue.empty()) {
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
		};
	      number_of_neighbours_computed++;
	    }
	  }   // next_neighbour_found or priority queue is empty
	  // in the latter case also the item priority queue is empty
        
	}
      }; // class Iterator_implementation
    }; // class iterator
  }; // class 

} // namespace CGAL

#endif  // CGAL_INCREMENTAL_NEIGHBOR_SEARCH_H
