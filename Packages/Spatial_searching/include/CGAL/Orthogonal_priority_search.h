// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-99 $
// release_date  : $CGAL_Date: 2003/05/23 $
//
// file          : include/CGAL/Orthogonal_priority_search.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef  ORTHOGONAL_PRIORITY_SEARCH
#define  ORTHOGONAL_PRIORITY_SEARCH
#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Euclidean_distance.h>

namespace CGAL {

template <class TreeTraits, 
	  class Distance=Euclidean_distance<typename TreeTraits::Point>, 
	  class Tree=Kd_tree<TreeTraits> >
class Orthogonal_priority_search {

public:

typedef typename TreeTraits::Point Point;
typedef typename TreeTraits::Point Query_item;
typedef typename TreeTraits::NT NT;
typedef typename Tree::Point_iterator Point_iterator;
typedef typename Tree::Node_handle Node_handle;

typedef std::pair<Point*,NT> Point_with_distance;
typedef std::pair<Node_handle,NT> Node_with_distance;


class iterator;

    typedef std::vector<Node_with_distance*> Node_with_distance_vector;

    typedef std::vector<Point_with_distance*> Point_with_distance_vector;

    typedef std::vector<NT> Distance_vector;

    iterator *start;
    iterator *past_the_end;

    public:

    // constructor
    Orthogonal_priority_search(Tree& tree,  
	Query_item& q, const Distance& tr=Distance(), NT Eps = NT(0.0), 
        bool search_nearest=true) 
    {
	start = new iterator(tree,q,tr,Eps,search_nearest);
        past_the_end = new iterator();
        
    };

    // destructor
    ~Orthogonal_priority_search() {
		delete start;
                delete past_the_end;
    };

    iterator begin() {
		return *start;
    }

    iterator end() {
		return *past_the_end;
    }

    void statistics() {
	start->statistics();
    }

    class iterator {

    public:

    typedef std::input_iterator_tag iterator_category;
    typedef Point_with_distance value_type;
    typedef int distance_type;

    class Iterator_implementation;
    Iterator_implementation *Ptr_implementation;


    public:

    // default constructor
    iterator() {Ptr_implementation=0;}

    int the_number_of_items_visited() {
        return Ptr_implementation->number_of_items_visited;
    }

    // constructor
    iterator(Tree& tree, Query_item& q, const Distance& tr=Distance(), NT eps=NT(0.0), 
    bool search_nearest=true){
        Ptr_implementation =
        new Iterator_implementation(tree, q, tr, eps, search_nearest);
    }

    // copy constructor
    iterator(const iterator& Iter) {
        Ptr_implementation = Iter.Ptr_implementation;
        if (Ptr_implementation != 0) Ptr_implementation->reference_count++;
    }

    Point_with_distance& operator* () {
                return *(*Ptr_implementation);
    }

    // prefix operator
    iterator& operator++() {
        ++(*Ptr_implementation);
        return *this;
    }

    // postfix operator
    std::auto_ptr<Point_with_distance> operator++(int) {
        std::auto_ptr<Point_with_distance> result = (*Ptr_implementation)++;
        return result;
    }


    bool operator==(const iterator& It) const {

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

    bool operator!=(const iterator& It) const {
        return !(*this == It);
    }

    void statistics () {
    	Ptr_implementation->statistics();
    }

    ~iterator() {
        if (Ptr_implementation != 0) {
                Ptr_implementation->reference_count--;
                if (Ptr_implementation->reference_count==0) {
                        delete Ptr_implementation;
                        Ptr_implementation = 0;
                }
        }
    }


    class Iterator_implementation {

    public:

    int number_of_neighbours_computed;
    int number_of_internal_nodes_visited;
    int number_of_leaf_nodes_visited;
    int number_of_items_visited;

    private:

    
    NT multiplication_factor;

    Point* query_point;

    int total_item_number;

    NT distance_to_root;

    bool search_nearest_neighbour;

    NT rd;

    class Priority_higher
    {
    public:

        bool search_nearest;

        Priority_higher(bool search_the_nearest_neighbour) {
                search_nearest = search_the_nearest_neighbour;
        } 

        //highest priority is smallest distance
        bool operator() (Node_with_distance* n1, Node_with_distance* n2) const
	{
                if (search_nearest) { return (n1->second > n2->second);}
                else {return (n2->second > n1->second);}
        }
    };

    class Distance_smaller
    {

    public:

        bool search_nearest;

        Distance_smaller(bool search_the_nearest_neighbour) {
                // Search_traits s;
                search_nearest = search_the_nearest_neighbour;
        } 

        //highest priority is smallest distance
        bool operator() (Point_with_distance* p1, Point_with_distance* p2) const
	{
		if (search_nearest) {return (p1->second > p2->second);}
                else {return (p2->second > p1->second);}
        }
    };

    std::priority_queue<Node_with_distance*, Node_with_distance_vector,
    Priority_higher>* PriorityQueue;

    std::priority_queue<Point_with_distance*, Point_with_distance_vector,
    Distance_smaller>* Item_PriorityQueue;

    Distance* Orthogonal_distance_instance;

    public:

    int reference_count;

    

    // constructor
    Iterator_implementation(Tree& tree, Query_item& q, const Distance& tr,
        NT Eps=NT(0.0), bool search_nearest=true)
    {
        PriorityQueue= new std::priority_queue<Node_with_distance*, 
	Node_with_distance_vector,
    	Priority_higher> 
        (Priority_higher(search_nearest));

        Item_PriorityQueue = new std::priority_queue<Point_with_distance*, 
	Point_with_distance_vector,
    	Distance_smaller> 
       (Distance_smaller(search_nearest));
       
        search_nearest_neighbour=search_nearest;
	reference_count=1;
        Orthogonal_distance_instance= new Distance(tr);
        multiplication_factor=
	Orthogonal_distance_instance->transformed_distance(NT(1.0)+Eps);

        // if (search_nearest) 
	distance_to_root=
	Orthogonal_distance_instance->min_distance_to_queryitem(q,
						*(tree.bounding_box()));
        // else distance_to_root=
   	// Orthogonal_Distance_instance->max_distance_to_queryitem(q,
	//					*(tree.bounding_box()));

        
        query_point = &q;

        total_item_number=tree.item_number();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
        number_of_neighbours_computed=0;



        Node_with_distance *The_Root = new Node_with_distance(tree.root(),
						distance_to_root);
        PriorityQueue->push(The_Root);

        // rd is the distance of the top of the priority queue to q
        rd=The_Root->second;
        Compute_the_next_nearest_neighbour();
    }

    // * operator
    Point_with_distance& operator* () {
			return *(Item_PriorityQueue->top());
    }

    // prefix operator
    Iterator_implementation& operator++() {
        // std::cout << "called ++" << std::endl;
        Delete_the_current_item_top();
        Compute_the_next_nearest_neighbour();
        return *this;
    }

    // postfix operator
    std::auto_ptr<Point_with_distance> operator++(int) {
        Point_with_distance Value = *(Item_PriorityQueue->top());
        std::auto_ptr<Point_with_distance>
        result(new Point_with_distance(Value));
        ++*this;
        return result;
    }

    // Print statistics of the general priority search process.
    void statistics () {
    	std::cout << "Orthogonal priority search statistics:" 
	<< std::endl;
    	std::cout << "Number of internal nodes visited:" 
	<< number_of_internal_nodes_visited << std::endl;
    	std::cout << "Number of leaf nodes visited:" 
	<< number_of_leaf_nodes_visited << std::endl;
    	std::cout << "Number of items visited:" 
	<< number_of_items_visited << std::endl;
        std::cout << "Number of neighbours computed:" 
	<< number_of_neighbours_computed << std::endl;
    }


    //destructor
    ~Iterator_implementation() {

        
	while (PriorityQueue->size()>0) {
                Node_with_distance* The_top=PriorityQueue->top();
                PriorityQueue->pop();
                delete The_top;
	};
	while (Item_PriorityQueue->size()>0) {
                Point_with_distance* The_top=Item_PriorityQueue->top();
                Item_PriorityQueue->pop();
                delete The_top;
        };
        delete PriorityQueue;
        delete Item_PriorityQueue;
	delete Orthogonal_distance_instance;
    }

    private:

    void Delete_the_current_item_top() {
        Point_with_distance* The_item_top=Item_PriorityQueue->top();
        Item_PriorityQueue->pop();
        delete The_item_top;
    }

    void Compute_the_next_nearest_neighbour() {

        // compute the next item
        bool next_neighbour_found=false;
        if (!(Item_PriorityQueue->empty())) {
        if (search_nearest_neighbour)
        	next_neighbour_found=
		(multiplication_factor*rd > Item_PriorityQueue->top()->second);
        else
		next_neighbour_found=
		(rd < multiplication_factor*Item_PriorityQueue->top()->second);
        }
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue->empty())) {
                Node_with_distance* The_node_top=PriorityQueue->top();
                Node_handle N= The_node_top->first;
                PriorityQueue->pop();
                delete The_node_top;
		NT copy_rd=rd;
                while (!(N->is_leaf())) { // compute new distance
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->cutting_dimension();
                        NT old_off, new_rd;
                        NT new_off =
                        (*query_point)[new_cut_dim] -
                        N->cutting_value();
                        if (new_off < NT(0.0)) {
				old_off=
                                (*query_point)[new_cut_dim]-N->low_value();
                                if (old_off>NT(0.0)) old_off=NT(0.0);
                                new_rd=
                                Orthogonal_distance_instance->
                                new_distance(copy_rd,old_off,new_off,new_cut_dim);
				assert(new_rd >= copy_rd);
				if (search_nearest_neighbour) {
                                	Node_with_distance *Upper_Child =
                                	new Node_with_distance(N->upper(), new_rd);
					PriorityQueue->push(Upper_Child);
                                	N=N->lower();
				}
				else {
                                	Node_with_distance *Lower_Child =
                                	new Node_with_distance(N->lower(), copy_rd);
					PriorityQueue->push(Lower_Child);
                                	N=N->upper();
					copy_rd=new_rd;
				}

                        }
                        else { // compute new distance
				old_off= N->high_value() -
                                (*query_point)[new_cut_dim];
                                if (old_off>NT(0.0)) old_off=NT(0.0);
                                new_rd=Orthogonal_distance_instance->
                                new_distance(copy_rd,old_off,new_off,new_cut_dim);  
				assert(new_rd >= copy_rd);
				if (search_nearest_neighbour) {
                                	Node_with_distance *Lower_Child =
                                	new Node_with_distance(N->lower(), new_rd);
                                	PriorityQueue->push(Lower_Child);
                                	N=N->upper();
				}
				else {
                                	Node_with_distance *Upper_Child =
                                	new Node_with_distance(N->upper(), copy_rd);
                                	PriorityQueue->push(Upper_Child);
                                	N=N->lower();
					copy_rd=new_rd;
				}
                        }
                }
                // n is a leaf
                number_of_leaf_nodes_visited++;
                if (N->size() > 0) {
                  for (Point_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
                        NT distance_to_query_point=
                        Orthogonal_distance_instance->
                        distance(*query_point,**it);
                        Point_with_distance *NN_Candidate=
                        new Point_with_distance(*it,distance_to_query_point);
                        Item_PriorityQueue->push(NN_Candidate);
                  };
                  // old top of PriorityQueue has been processed,
                  // hence update rd
                
                  if (!(PriorityQueue->empty()))  {
                        rd = PriorityQueue->top()->second;
                        if (search_nearest_neighbour)
				next_neighbour_found =
                  (multiplication_factor*rd > 
		   Item_PriorityQueue->top()->second);
                        else
				next_neighbour_found =
                  (multiplication_factor*rd < 
		   Item_PriorityQueue->top()->second);
                  }
                  else // priority queue empty => last neighbour found
                  {
                        next_neighbour_found=true;
                  };

                  number_of_neighbours_computed++;
           }
        }   // next_neighbour_found or priority queue is empty
        // in the latter case also the item priority quee is empty
    }
}; // class Iterator_implementaion
}; // class iterator
}; // class 

template <class Traits, class Query_item, class Distance>
void swap (typename Orthogonal_priority_search<Traits, 
				Query_item, Distance>::iterator& x,
        typename Orthogonal_priority_search<Traits, 
				Query_item, Distance>::iterator& y) {
        typename Orthogonal_priority_search<Traits, 
		Query_item, Distance>::iterator::Iterator_implementation
        *tmp = x.Ptr_implementation;
        x.Ptr_implementation  = y.Ptr_implementation;
        y.Ptr_implementation = tmp;
}



} // namespace CGAL


#endif  // ORTHOGONAL_PRIORITY_SEARCH_H
