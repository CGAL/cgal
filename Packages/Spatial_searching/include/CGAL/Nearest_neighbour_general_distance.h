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
// release       :
// release_date  :
//
// file          : include/CGAL/Nearest_neighbour_general_distance.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef  NEAREST_NEIGHBOUR_GENERAL_DISTANCE_H
#define  NEAREST_NEIGHBOUR_GENERAL_DISTANCE_H

#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Internal_node.h>
#include <CGAL/Box.h>

// copy to example.cpp
// #include <CGAL/Kd_tree_traits_point.h>
// #include <CGAL/Weighted_Minkowski_distance.h>

using std::list; // to avoid compiler crash on MSVC++

namespace CGAL {

template <class Tree_traits, class Search_traits, class Distance>
class Nearest_neighbour_general_distance {

public:

typedef typename Tree_traits::Item Item;
typedef typename Tree_traits::NT NT;
typedef Item** Item_iterator;
typedef Base_node<Tree_traits> Node;
typedef Binary_search_tree<Tree_traits> Tree;
typedef Box<NT> Node_box;    


class Cell 
    {
    private:

    	Node_box* the_box;
    	Node*     the_node;

    public:

        // constructor
        Cell (Node_box* Nb, Node* N)
	{
		the_box = Nb;
		the_node = N;
	}

        Node_box* box() {return the_box;};
        Node*    node() {return the_node;};
        

	~Cell() {}

    };

    
typedef Tree_traits::Item_with_distance Item_with_distance;
typedef std::pair<Cell*,NT> Cell_with_distance;

class iterator;

class Priority_higher
    {
    public:

        bool search_nearest;

        Priority_higher() {
                Search_traits s;
                search_nearest = s.Search_nearest();
        }
        //highest priority is smallest distance
        bool operator() (Cell_with_distance* n1, Cell_with_distance* n2) const {
                if (search_nearest) { return (n1->second > n2->second);}
                else {return (n2->second > n1->second);}
        }
    };

class Distance_smaller
    {

    public:

        bool search_nearest;

        Distance_smaller() {
                Search_traits s;
                search_nearest = s.Search_nearest();
        }

        //highest priority is smallest distance
        bool operator() (Item_with_distance* p1, Item_with_distance* p2) const {
		if (search_nearest) {return (p1->second > p2->second);}
                else {return (p2->second > p1->second);}
        }
    };

    typedef std::vector<Cell_with_distance*> Cell_with_distance_vector;

    typedef std::vector<Item_with_distance*> Item_with_distance_vector;

    typedef std::vector<NT> Distance_vector;

    iterator *start;
    iterator *past_the_end;

    public:

    // constructor
    Nearest_neighbour_general_distance(Tree& tree, Item& q, Distance& tr,
    NT Eps=0.0)
    {
        start = new iterator(tree,q,tr,Eps);
        past_the_end = new iterator();
    };

    // destructor
    ~Nearest_neighbour_general_distance() {
		delete start;
        delete past_the_end;
    };

    iterator begin() {
		return *start;
    }

    iterator end() {
		return *past_the_end;
    }

    class iterator {

    public:

    typedef std::input_iterator_tag iterator_category;
    typedef Item_with_distance value_type;
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
    iterator(Tree& tree, Item& q, Distance& tr, NT eps=0.0){
        Ptr_implementation =
        new Iterator_implementation(tree, q, tr, eps);
    }

    // copy constructor
    iterator(iterator& Iter) {
        Ptr_implementation = Iter.Ptr_implementation;
        if (Ptr_implementation != 0) Ptr_implementation->reference_count++;
    }

    Item_with_distance& operator* () {
                return *(*Ptr_implementation);
    }

    // prefix operator
    iterator& operator++() {
        ++(*Ptr_implementation);
        return *this;
    }

    // postfix operator
    std::auto_ptr<Item_with_distance> operator++(int) {
        std::auto_ptr<Item_with_distance> result = (*Ptr_implementation)++;
        return result;
    }

    bool operator==(const iterator& It) const {

        if (
                ((Ptr_implementation == 0) || Ptr_implementation->Item_PriorityQueue.empty()) &&
                ((It.Ptr_implementation == 0) ||  It.Ptr_implementation->Item_PriorityQueue.empty())
        )
        return true;
        // else
        return (Ptr_implementation == It.Ptr_implementation);
    }

    bool operator!=(const iterator& It) const {
        return !(*this == It);
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

    Item* query_point;

    int total_item_number;

    NT distance_to_root;

    bool search_nearest;

    NT rd;

    std::priority_queue<Cell_with_distance*, Cell_with_distance_vector,
    Priority_higher> PriorityQueue;

    Distance* Distance_instance;

    public:

    int reference_count;

    std::priority_queue<Item_with_distance*, Item_with_distance_vector,
    Distance_smaller> Item_PriorityQueue;

    // constructor
    Iterator_implementation(Tree& tree, Item& q, Distance& tr,
        NT Eps=0.0)
    {

        Search_traits s;
        search_nearest = s.Search_nearest();

        std::cout << "search_nearest=" << search_nearest << std::endl;

	reference_count=1;
        Distance_instance=&tr;
        multiplication_factor=Distance_instance->transformed_distance(1.0+Eps);

        Node_box *bounding_box = new Node_box(*(tree.bounding_box()));

        std::cout << "Bounding box tree is " << *(tree.bounding_box()) << std::endl;

        if (search_nearest) distance_to_root=
        Distance_instance->lower_bound_distance_to_box(q,*bounding_box);
        else distance_to_root=
   	Distance_instance->upper_bound_distance_to_box(q,*bounding_box);

        std::cout << "distance_to_root=" << distance_to_root << std::endl;

        query_point = &q;

        total_item_number=tree.item_number();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
        number_of_neighbours_computed=0;


        // code form Nearest_Neighbour_PQ
        // Node_with_distance *The_Root = new Node_with_distance(tree.root(),distance_to_root);

        Cell *Root_Cell = new Cell(bounding_box,tree.root());
        Cell_with_distance  *The_Root = new Cell_with_distance(Root_Cell,distance_to_root);

        PriorityQueue.push(The_Root);

        // rd is the distance of the top of the priority queue to q
        rd=The_Root->second;
        Compute_the_next_nearest_neighbour();
    }

    // * operator
    Item_with_distance& operator* () {
			return *(Item_PriorityQueue.top());
    }

    // prefix operator
    Iterator_implementation& operator++() {
        // std::cout << "called ++" << std::endl;
        Delete_the_current_item_top();
        Compute_the_next_nearest_neighbour();
        return *this;
    }

    // postfix operator
    std::auto_ptr<Item_with_distance> operator++(int) {
        Item_with_distance Value = *(Item_PriorityQueue.top());
        std::auto_ptr<Item_with_distance>
        result(new Item_with_distance(Value));
        ++*this;
        return result;
    }


    //destructor
    ~Iterator_implementation() {

        // std::cout << "called destructor" << std::endl;
	while (PriorityQueue.size()>0) {
                Cell_with_distance* The_top=PriorityQueue.top();
                PriorityQueue.pop();
                delete The_top->first->box();
		delete The_top->first;
                delete The_top;
	};
	while (Item_PriorityQueue.size()>0) {
                Item_with_distance* The_top=Item_PriorityQueue.top();
                Item_PriorityQueue.pop();
                delete The_top;
        };
    }

    private:

    void Delete_the_current_item_top() {
        Item_with_distance* The_item_top=Item_PriorityQueue.top();
        Item_PriorityQueue.pop();
        delete The_item_top;
    }

    void Compute_the_next_nearest_neighbour() {

        // compute the next item
        bool next_neighbour_found=false;
        if (!(Item_PriorityQueue.empty())) {
        if (search_nearest)
        	next_neighbour_found=(multiplication_factor*rd > Item_PriorityQueue.top()->second);
        else
		next_neighbour_found=(multiplication_factor*rd < Item_PriorityQueue.top()->second);
        }
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue.empty())) {
                
                Cell_with_distance* The_node_top=PriorityQueue.top();
                Node* N= The_node_top->first->node();
                Node_box* B= The_node_top->first->box();
                PriorityQueue.pop();
                delete The_node_top->first;
                delete The_node_top;

                while (!(N->is_leaf())) { // compute new distances
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->separator()->cutting_dimension();
                        NT  new_cut_val=N->separator()->cutting_value();
                        // if (!(search_nearest)) std::cout << "Box is " << *B << std::endl;
			Node_box* lower_box = new Node_box(*B);
                        Node_box* upper_box = lower_box->split(new_cut_dim, new_cut_val);
			delete B;
			if (search_nearest) {
                        	NT distance_to_box_lower =
				Distance_instance->lower_bound_distance_to_box(*query_point, *lower_box);
                                NT distance_to_box_upper =
				Distance_instance->lower_bound_distance_to_box(*query_point, *upper_box);
				if (distance_to_box_lower <= distance_to_box_upper) {
					Cell* C_upper = new Cell(upper_box, N->upper());
					Cell_with_distance *Upper_Child =
                                	new Cell_with_distance(C_upper,distance_to_box_upper);
					PriorityQueue.push(Upper_Child);
                                	N=N->lower();
				        B=lower_box;
				}
				else {
					Cell* C_lower = new Cell(lower_box, N->lower());
					Cell_with_distance *Lower_Child =
                                	new Cell_with_distance(C_lower,distance_to_box_lower);
					PriorityQueue.push(Lower_Child);
                                	N=N->upper();
				        B=upper_box;
				}
                        }
			else { // search furthest
                        	NT distance_to_box_lower =
				Distance_instance->upper_bound_distance_to_box(*query_point, *lower_box);
                                NT distance_to_box_upper =
				Distance_instance->upper_bound_distance_to_box(*query_point, *upper_box);
				if (distance_to_box_lower >= distance_to_box_upper) {
					Cell* C_upper = new Cell(upper_box, N->upper());
					Cell_with_distance *Upper_Child =
                                	new Cell_with_distance(C_upper,distance_to_box_upper);
					PriorityQueue.push(Upper_Child);
                                	N=N->lower();
				        B=lower_box;
				}
				else {
					Cell* C_lower = new Cell(lower_box, N->lower());
					Cell_with_distance *Lower_Child =
                                	new Cell_with_distance(C_lower,distance_to_box_lower);
					PriorityQueue.push(Lower_Child);
                                	N=N->upper();
				        B=upper_box;
				}
                        }
                }
                // n is a leaf
                delete B;
                number_of_leaf_nodes_visited++;
                for (Item_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
                        NT distance_to_query_point=
                        Distance_instance->
                        distance(*query_point,**it);
                        Item_with_distance *NN_Candidate=
                        new Item_with_distance(*it,distance_to_query_point);
                        Item_PriorityQueue.push(NN_Candidate);
                }
                // old top of PriorityQueue has been processed,
                // hence update rd
                if (!PriorityQueue.empty()) {
                        rd = PriorityQueue.top()->second;
                        if (search_nearest)
				next_neighbour_found =
                (multiplication_factor*rd > Item_PriorityQueue.top()->second);
                        else
				next_neighbour_found =
                (multiplication_factor*rd < Item_PriorityQueue.top()->second);
                }
                else // priority queue empty => last neighbour found
                {
                        next_neighbour_found=true;
                };
                number_of_neighbours_computed++;
        }   // next_neighbour_found or priority queue is empty
        // in the latter case also the item priority quee is empty
    }
}; // class Iterator_implementaion
}; // class iterator
}; // class Nearest neighbour_L2

template <class Tree_traits, class Search_traits, class Distance>
void swap (typename Nearest_neighbour_general_distance<Tree_traits, Search_traits, Distance>::iterator& x,
        typename Nearest_neighbour_general_distance<Tree_traits, Search_traits, Distance>::iterator& y) {
        typename Nearest_neighbour_general_distance<Tree_traits, Search_traits, Distance>::iterator::Iterator_implementation
        *tmp = x.Ptr_implementation;
        x.Ptr_implementation  = y.Ptr_implementation;
        y.Ptr_implementation = tmp;
};

} // namespace CGAL


#endif  // NEAREST_NEIGHBOUR_GENERAL_DISTANCE_H
