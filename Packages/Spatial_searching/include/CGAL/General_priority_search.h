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
// file          : include/CGAL/General_priority_search.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef  GENERAL_PRIORITY_SEARCH_H
#define  GENERAL_PRIORITY_SEARCH_H
#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree_rectangle.h>
namespace CGAL {

template <class Tree_traits, class Query_item, class Distance>
class General_priority_search {

public:

typedef typename Tree_traits::Item Item;
typedef typename Tree_traits::NT NT;
typedef Item** Item_iterator;
typedef Kd_tree_node<Tree_traits> Node;
typedef Kd_tree<Tree_traits> Tree;
typedef Kd_tree_rectangle<NT> Node_box;   

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

    
typedef typename Tree_traits::Item_with_distance Item_with_distance;
typedef std::pair<Cell*,NT> Cell_with_distance;

// this forward declaration may problems for g++ 
class iterator;


    typedef std::vector<Cell_with_distance*> Cell_with_distance_vector;

    typedef std::vector<Item_with_distance*> Item_with_distance_vector;

    typedef std::vector<NT> Distance_vector;

    iterator *start;
    iterator *past_the_end;

    public:

    // constructor
    General_priority_search(Tree& tree, Query_item& q, Distance& tr,
    NT Eps=0.0, bool search_nearest=true)
    {
        start = new iterator(tree,q,tr,Eps,search_nearest);
        past_the_end = new iterator();
    };

    // destructor
    ~General_priority_search() {
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
    iterator(Tree& tree, Query_item& q, Distance& tr, NT eps=0.0, bool search_nearest=true) {
        Ptr_implementation =
        new Iterator_implementation(tree, q, tr, eps, search_nearest);
    }

    // copy constructor
    iterator(const iterator& Iter) {
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
                ((Ptr_implementation == 0) || 
		Ptr_implementation->Item_PriorityQueue->empty()) &&
                ((It.Ptr_implementation == 0) ||  
		It.Ptr_implementation->Item_PriorityQueue->empty())
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

    



    private:

    NT multiplication_factor;

    Query_item* query_point;

    int total_item_number;

    NT distance_to_root;

    bool search_nearest;

    NT rd;

    class Priority_higher
    {
    public:

        bool search_nearest;

        Priority_higher(bool search_the_nearest_neighbour) {
                search_nearest = search_the_nearest_neighbour;
        }

        //highest priority is smallest distance
        bool operator() (Cell_with_distance* n1, Cell_with_distance* n2) const 
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
        bool operator() (Item_with_distance* p1, Item_with_distance* p2) const 
        {
		if (search_nearest) {return (p1->second > p2->second);}
                else {return (p2->second > p1->second);}
        }
    };

    std::priority_queue<Cell_with_distance*, Cell_with_distance_vector,
    Priority_higher>* PriorityQueue;

    std::priority_queue<Item_with_distance*, Item_with_distance_vector,
    Distance_smaller>* Item_PriorityQueue;
    

    Distance* Distance_instance;

    public:

    int reference_count;

    int number_of_internal_nodes_visited;
    int number_of_leaf_nodes_visited;
    int number_of_items_visited;
    int number_of_neighbours_computed;

    // constructor
    Iterator_implementation(Tree& tree, Query_item& q, Distance& tr,
        NT Eps=0.0, bool search_nearest=true)
    {
        
        PriorityQueue = new std::priority_queue<Cell_with_distance*, 
	Cell_with_distance_vector, Priority_higher> 
        (Priority_higher(search_nearest));

        Item_PriorityQueue= new std::priority_queue<Item_with_distance*, 
	Item_with_distance_vector,
    	Distance_smaller>
	(Distance_smaller(search_nearest));

	reference_count=1;
        Distance_instance=&tr;
        multiplication_factor=
	Distance_instance->transformed_distance(1.0+Eps);

        Node_box *bounding_box = new Node_box(*(tree.bounding_box()));


        if (search_nearest) distance_to_root=
        Distance_instance->min_distance_to_queryitem(q,*bounding_box);
        else distance_to_root=
   	Distance_instance->max_distance_to_queryitem(q,*bounding_box);

        

        query_point = &q;

        total_item_number=tree.item_number();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
        number_of_neighbours_computed=0;


        Cell *Root_Cell = new Cell(bounding_box,tree.root());
        Cell_with_distance  *The_Root = 
	new Cell_with_distance(Root_Cell,distance_to_root);

        PriorityQueue->push(The_Root);

        // rd is the distance of the top of the priority queue to q
        rd=The_Root->second;
        Compute_the_next_nearest_neighbour();
    }

    // * operator
    Item_with_distance& operator* () {
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
    std::auto_ptr<Item_with_distance> operator++(int) {
        Item_with_distance Value = *(Item_PriorityQueue->top());
        std::auto_ptr<Item_with_distance>
        result(new Item_with_distance(Value));
        ++*this;
        return result;
    }


    // Print statistics of the general priority search process.
    void statistics () {
    	std::cout << "General priority search statistics:" << std::endl;
    	std::cout << "Number of internal nodes visited:" << number_of_internal_nodes_visited << std::endl;
    	std::cout << "Number of leaf nodes visited:" << number_of_leaf_nodes_visited << std::endl;
    	std::cout << "Number of items visited:" << number_of_items_visited << std::endl;
        std::cout << "Number of neighbours computed:" << number_of_neighbours_computed << std::endl;
    }

    //destructor
    ~Iterator_implementation() {

        // std::cout << "called destructor" << std::endl;
	while (PriorityQueue->size()>0) {
                Cell_with_distance* The_top=PriorityQueue->top();
                PriorityQueue->pop();
                delete The_top->first->box();
		delete The_top->first;
                delete The_top;
	};
	while (Item_PriorityQueue->size()>0) {
                Item_with_distance* The_top=Item_PriorityQueue->top();
                Item_PriorityQueue->pop();
                delete The_top;
        };
    }

    private:

    void Delete_the_current_item_top() {
        Item_with_distance* The_item_top=Item_PriorityQueue->top();
        Item_PriorityQueue->pop();
        delete The_item_top;
    }

    void Compute_the_next_nearest_neighbour() {

        // compute the next item
        bool next_neighbour_found=false;
        if (!(Item_PriorityQueue->empty())) {
        if (search_nearest)
        	next_neighbour_found=
		(multiplication_factor*rd > Item_PriorityQueue->top()->second);
        else
		next_neighbour_found=
		(multiplication_factor*rd < Item_PriorityQueue->top()->second);
        }
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue->empty())) {
                
                Cell_with_distance* The_node_top=PriorityQueue->top();
                Node* N= The_node_top->first->node();
                Node_box* B= The_node_top->first->box();
                PriorityQueue->pop();
                delete The_node_top->first;
                delete The_node_top;

                while (!(N->is_leaf())) { // compute new distances
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->separator()->cutting_dimension();
                        NT  new_cut_val=N->separator()->cutting_value();
                        
			Node_box* lower_box = new Node_box(*B);
                        Node_box* upper_box = 
			lower_box->split(new_cut_dim, new_cut_val);
			delete B;
			if (search_nearest) {
NT distance_to_box_lower =
Distance_instance->min_distance_to_queryitem(*query_point, *lower_box);
NT distance_to_box_upper =
Distance_instance->max_distance_to_queryitem(*query_point, *upper_box);
if (distance_to_box_lower <= distance_to_box_upper) {
	Cell* C_upper = new Cell(upper_box, N->upper());
	Cell_with_distance *Upper_Child =
	new Cell_with_distance(C_upper,distance_to_box_upper);
	PriorityQueue->push(Upper_Child);
	N=N->lower();
	B=lower_box;
}
else {
	Cell* C_lower = new Cell(lower_box, N->lower());
	Cell_with_distance *Lower_Child =
	new Cell_with_distance(C_lower,distance_to_box_lower);
	PriorityQueue->push(Lower_Child);
	N=N->upper();
	B=upper_box;
}
                        }
			else { // search furthest
NT distance_to_box_lower =
Distance_instance->max_distance_to_queryitem(*query_point, *lower_box);
NT distance_to_box_upper =
Distance_instance->max_distance_to_queryitem(*query_point, *upper_box);
if (distance_to_box_lower >= distance_to_box_upper) {
	Cell* C_upper = new Cell(upper_box, N->upper());
	Cell_with_distance *Upper_Child =
	new Cell_with_distance(C_upper,distance_to_box_upper);
	PriorityQueue->push(Upper_Child);
	N=N->lower();
	B=lower_box;
}
else {
	Cell* C_lower = new Cell(lower_box, N->lower());
	Cell_with_distance *Lower_Child =
        new Cell_with_distance(C_lower,distance_to_box_lower);
	PriorityQueue->push(Lower_Child);
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
                        Item_PriorityQueue->push(NN_Candidate);
                }
                // old top of PriorityQueue has been processed,
                // hence update rd
                if (!PriorityQueue->empty()) {
                        rd = PriorityQueue->top()->second;
                        if (search_nearest)
				next_neighbour_found =
                (multiplication_factor*rd > Item_PriorityQueue->top()->second);
                        else
				next_neighbour_found =
                (multiplication_factor*rd < Item_PriorityQueue->top()->second);
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
}; // class 

template <class Traits, class Query_item, class Distance>
void swap (typename General_priority_search<Traits, 
				Query_item, Distance>::iterator& x,
        typename General_priority_search<Traits, 
				Query_item, Distance>::iterator& y) {
        typename General_priority_search<Traits, 
		Query_item, Distance>::iterator::Iterator_implementation
        *tmp = x.Ptr_implementation;
        x.Ptr_implementation  = y.Ptr_implementation;
        y.Ptr_implementation = tmp;
}

} // namespace CGAL


#endif  // GENERAL_PRIORITY_SEARCH_H
