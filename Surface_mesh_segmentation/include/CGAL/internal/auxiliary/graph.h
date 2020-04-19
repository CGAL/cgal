// Copyright (c) 2001  Yuri Boykov
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// Re-licensed for CGAL distribution to:
// SPDX-License-Identifier: GPL-3.0-or-later
// Original license is:
// SPDX-License-Identifier: GPL-2.0-or-later
/*
###################################################################
#                                                                 #
#    MAXFLOW - software for computing mincut/maxflow in a graph   #
#                        Version 2.21                             #
#    http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html     #
#                                                                 #
#    Yuri Boykov (yuri@csd.uwo.ca)                                #
#    Vladimir Kolmogorov (v.kolmogorov@cs.ucl.ac.uk)              #
#    2001                                                         #
#                                                                 #
###################################################################

1. Introduction.

This software library implements the maxflow algorithm
described in

        An Experimental Comparison of Min-Cut/Max-Flow Algorithms
        for Energy Minimization in Vision.
        Yuri Boykov and Vladimir Kolmogorov.
        In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI),
        September 2004

This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
at Siemens Corporate Research. To make it available for public use,
it was later reimplemented by Vladimir Kolmogorov based on open publications.

If you use this software for research purposes, you should cite
the aforementioned paper in any resulting publication.

Tested under windows, Visual C++ 6.0 compiler and unix (SunOS 5.8
and RedHat Linux 7.0, GNU c++ compiler).

##################################################################

2. Licence.

Copyright UCL Business PLC

This program is available under dual licence:

1) Under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
Note that any program that incorporates the code under this licence must, under the terms of the GNU GPL, be released under a licence compatible with the GPL. GNU GPL does not permit incorporating this program into proprietary programs. If you wish to do this, please see the alternative licence available below.
GNU General Public License can be found at http://www.gnu.org/licenses/old-licenses/gpl-2.0.html

2) Proprietary Licence from UCL Business PLC.
To enable programers to include the MaxFlow software in a proprietary system (which is not allowed by the GNU GPL), this licence gives you the right to incorporate the software in your program and distribute under any licence of your choosing. The full terms of the licence and applicable fee, are available from the Licensors at: http://www.uclb-elicensing.com/optimisation_software/maxflow_computervision.html

##################################################################

3. Graph representation.

There are two versions of the algorithm using different
graph representations (adjacency list and forward star).
The former one uses more than twice as much memory as the
latter one but is 10-20% faster.

Memory allocation (assuming that all capacities are 'short' - 2 bytes):

                 |   Nodes    |   Arcs
------------------------------------------
Adjacency list   | *24 bytes  | *14 bytes
Forward star     | *28 bytes  |  6 bytes

(* means that often it should be rounded up to be a multiple of 4
- some compilers (e.g. Visual C++) seem to round up elements
of arrays unless the are structures containing only char[].)

Note that arcs are always added in pairs - in forward and reverse directions.
Arcs between nodes and terminals (the source and the sink) are
not stored as arcs, but rather as a part of nodes.

The assumption for the forward star representation is that
the maximum number of arcs per node (except the source
and the sink) is much less than MF_ARC_BLOCK_SIZE (1024 by default).

Both versions have the same interface.

##################################################################

4. Example usage.

This section shows how to use the library to compute
a minimum cut on the following graph:

                        SOURCE
                       /       \
                     1/         \2
                     /      3    \
                   node0 -----> node1
                     |   <-----   |
                     |      4     |
                     \            /
                     5\          /6
                       \        /
                          SINK

///////////////////////////////////////////////////

#include <stdio.h>
#include "graph.h"

void main()
{
        Graph::node_id nodes[2];
        Graph *g = new Graph();

        nodes[0] = g -> add_node();
        nodes[1] = g -> add_node();
        g -> set_tweights(nodes[0], 1, 5);
        g -> set_tweights(nodes[1], 2, 6);
        g -> add_edge(nodes[0], nodes[1], 3, 4);

        Graph::flowtype flow = g -> maxflow();

        printf("Flow = %d\n", flow);
        printf("Minimum cut:\n");
        if (g->what_segment(nodes[0]) == Graph::SOURCE)
                printf("node0 is in the SOURCE set\n");
        else
                printf("node0 is in the SINK set\n");
        if (g->what_segment(nodes[1]) == Graph::SOURCE)
                printf("node1 is in the SOURCE set\n");
        else
                printf("node1 is in the SINK set\n");

        delete g;
}

///////////////////////////////////////////////////
*/

/* block.h */
/*
        Template classes Block and DBlock
        Implement adding and deleting items of the same type in blocks.

        If there there are many items then using Block or DBlock
        is more efficient than using 'new' and 'delete' both in terms
        of memory and time since
        (1) On some systems there is some minimum amount of memory
            that 'new' can allocate (e.g., 64), so if items are
            small that a lot of memory is wasted.
        (2) 'new' and 'delete' are designed for items of varying size.
            If all items has the same size, then an algorithm for
            adding and deleting can be made more efficient.
        (3) All Block and DBlock functions are inline, so there are
            no extra function calls.

        Differences between Block and DBlock:
        (1) DBlock allows both adding and deleting items,
            whereas Block allows only adding items.
        (2) Block has an additional operation of scanning
            items added so far (in the order in which they were added).
        (3) Block allows to allocate several consecutive
            items at a time, whereas DBlock can add only a single item.

        Note that no constructors or destructors are called for items.

        Example usage for items of type 'MyType':

        ///////////////////////////////////////////////////
        #include "block.h"
        #define BLOCK_SIZE 1024

#include <CGAL/license/Surface_mesh_segmentation.h>

        typedef struct { int a, b; } MyType;
        MyType *ptr, *array[10000];

        ...

        Block<MyType> *block = new Block<MyType>(BLOCK_SIZE);

        // adding items
        for (int i=0; i<sizeof(array); i++)
        {
                ptr = block -> New();
                ptr -> a = ptr -> b = rand();
        }

        // reading items
        for (ptr=block->ScanFirst(); ptr; ptr=block->ScanNext())
        {
                printf("%d %d\n", ptr->a, ptr->b);
        }

        delete block;

        ...

        DBlock<MyType> *dblock = new DBlock<MyType>(BLOCK_SIZE);

        // adding items
        for (int i=0; i<sizeof(array); i++)
        {
                array[i] = dblock -> New();
        }

        // deleting items
        for (int i=0; i<sizeof(array); i+=2)
        {
                dblock -> Delete(array[i]);
        }

        // adding items
        for (int i=0; i<sizeof(array); i++)
        {
                array[i] = dblock -> New();
        }

        delete dblock;

        ///////////////////////////////////////////////////

        Note that DBlock deletes items by marking them as
        empty (i.e., by adding them to the list of free items),
        so that this memory could be used for subsequently
        added items. Thus, at each moment the memory allocated
        is determined by the maximum number of items allocated
        simultaneously at earlier moments. All memory is
        deallocated only when the destructor is called.
*/

#ifndef __MAXFLOW_BLOCK_H__
#define __MAXFLOW_BLOCK_H__

#  if defined(CGAL_LICENSE_WARNING)
     CGAL_pragma_warning("\nYou use the MaxFlow package of Vladimir Kolmogorov under the terms of the GPLv2+.")
#  endif // CGAL_LICENSE_WARNING

#  ifdef CGAL_LICENSE_ERROR
#    error "You use the the MaxFlow package of Vladimir Kolmogorov under the terms of the GPLv2+.\
You get this error, as you defined CGAL_LICENSE_ERROR."
#  endif // CGAL_LICENSE_ERROR


#include <stdlib.h>

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

template <class Type> class Block
{
public:
  /* Constructor. Arguments are the block size and
     (optionally) the pointer to the function which
     will be called if allocation failed; the message
     passed to this function is "Not enough memory!" */
  Block(int size, void (*err_function)(const char *) = nullptr) {
    first = last = nullptr;
    block_size = size;
    error_function = err_function;
  }

  /* Destructor. Deallocates all items added so far */
  ~Block() {
    while (first) {
      block *next = first -> next;
      delete[] ((char*)first);
      first = next;
    }
  }

  /* Allocates 'num' consecutive items; returns pointer
     to the first item. 'num' cannot be greater than the
     block size since items must fit in one block */
  Type *New(int num = 1) {
    Type *t;

    if (!last || last->current + num > last->last) {
      if (last && last->next) last = last -> next;
      else {
        block *next = (block *) new char [sizeof(block) + (block_size-1)*sizeof(Type)];
        if (!next) {
          if (error_function) (*error_function)("Not enough memory!");
          exit(1);
        }
        if (last) last -> next = next;
        else first = next;
        last = next;
        last -> current = & ( last -> data[0] );
        last -> last = last -> current + block_size;
        last -> next = nullptr;
      }
    }

    t = last -> current;
    last -> current += num;
    return t;
  }

  /* Returns the first item (or nullptr, if no items were added) */
  Type *ScanFirst() {
    for (scan_current_block=first; scan_current_block;
         scan_current_block = scan_current_block->next) {
      scan_current_data = & ( scan_current_block -> data[0] );
      if (scan_current_data < scan_current_block -> current) return scan_current_data
            ++;
    }
    return nullptr;
  }

  /* Returns the next item (or nullptr, if all items have been read)
     Can be called only if previous ScanFirst() or ScanNext()
     call returned not nullptr. */
  Type *ScanNext() {
    while (scan_current_data >= scan_current_block -> current) {
      scan_current_block = scan_current_block -> next;
      if (!scan_current_block) return nullptr;
      scan_current_data = & ( scan_current_block -> data[0] );
    }
    return scan_current_data ++;
  }

  /* Marks all elements as empty */
  void Reset() {
    block *b;
    if (!first) return;
    for (b=first; ; b=b->next) {
      b -> current = & ( b -> data[0] );
      if (b == last) break;
    }
    last = first;
  }

  /***********************************************************************/

private:

  typedef struct block_st {
    Type                                        *current, *last;
    struct block_st                        *next;
    Type                                        data[1];
  } block;

  int                block_size;
  block        *first;
  block        *last;

  block        *scan_current_block;
  Type        *scan_current_data;

  void        (*error_function)(const char *);
};

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

template <class Type> class DBlock
{
public:
  /* Constructor. Arguments are the block size and
     (optionally) the pointer to the function which
     will be called if allocation failed; the message
     passed to this function is "Not enough memory!" */
  DBlock(int size, void (*err_function)(const char *) = nullptr) {
    first = nullptr;
    first_free = nullptr;
    block_size = size;
    error_function = err_function;
  }

  /* Destructor. Deallocates all items added so far */
  ~DBlock() {
    while (first) {
      block *next = first -> next;
      delete[] ((char*)first);
      first = next;
    }
  }

  /* Allocates one item */
  Type *New() {
    block_item *item;

    if (!first_free) {
      block *next = first;
      first = (block *) new char [sizeof(block) + (block_size-1)*sizeof(block_item)];
      if (!first) {
        if (error_function) (*error_function)("Not enough memory!");
        exit(1);
      }
      first_free = & (first -> data[0] );
      for (item=first_free; item<first_free+block_size-1; item++)
        item -> next_free = item + 1;
      item -> next_free = nullptr;
      first -> next = next;
    }

    item = first_free;
    first_free = item -> next_free;
    return (Type *) item;
  }

  /* Deletes an item allocated previously */
  void Delete(Type *t) {
    ((block_item *) t) -> next_free = first_free;
    first_free = (block_item *) t;
  }

  /***********************************************************************/

private:

  typedef union block_item_st {
    Type                        t;
    block_item_st        *next_free;
  } block_item;

  typedef struct block_st {
    struct block_st                        *next;
    block_item                                data[1];
  } block;

  int                        block_size;
  block                *first;
  block_item        *first_free;

  void        (*error_function)(const char *);
};


#endif



/* graph.h */
/*
        This software library implements the maxflow algorithm
        described in

                An Experimental Comparison of Min-Cut/Max-Flow Algorithms
                for Energy Minimization in Vision.
                Yuri Boykov and Vladimir Kolmogorov.
                In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI),
                September 2004

        This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
        at Siemens Corporate Research. To make it available for public use,
        it was later reimplemented by Vladimir Kolmogorov based on open publications.

        If you use this software for research purposes, you should cite
        the aforementioned paper in any resulting publication.

        ----------------------------------------------------------------

        For description, license, example usage, discussion of graph representation        and memory usage see README.TXT.
*/

#ifndef __MAXFLOW_GRAPH_H__
#define __MAXFLOW_GRAPH_H__

//#include "block.h"
#include <stdio.h>
/*
        Nodes, arcs and pointers to nodes are
        added in blocks for memory and time efficiency.
        Below are numbers of items in blocks
*/
#define MF_NODE_BLOCK_SIZE 512
#define MF_ARC_BLOCK_SIZE 1024
#define MF_NODEPTR_BLOCK_SIZE 128

template <std::size_t size>
struct Int_to_ptr;

template<> struct Int_to_ptr<sizeof(int)> {
  typedef int type;
};
#if INT_MAX != LONG_MAX
template<> struct Int_to_ptr<sizeof(long)> {
  typedef long type;
};
#else
template<> struct Int_to_ptr<sizeof(long long)> {
  typedef long long type;
};
#endif


class Graph
{
public:
  typedef enum {
    SOURCE        = 0,
    SINK        = 1
  } termtype; /* terminals */

  /* Type of edge weights.
     Can be changed to char, int, float, double, ... */
  typedef double captype;
  /* Type of total flow */
  typedef double flowtype;

  typedef void * node_id;

  /* interface functions */

  /* Constructor. Optional argument is the pointer to the
     function which will be called if an error occurs;
     an error message is passed to this function. If this
     argument is omitted, exit(1) will be called. */
  Graph(void (*err_function)(const char *) = nullptr);

  /* Destructor */
  ~Graph();

  /* Adds a node to the graph */
  node_id add_node();

  /* Adds a bidirectional edge between 'from' and 'to'
     with the weights 'cap' and 'rev_cap' */
  void add_edge(node_id from, node_id to, captype cap, captype rev_cap);

  /* Sets the weights of the edges 'SOURCE->i' and 'i->SINK'
     Can be called at most once for each node before any call to 'add_tweights'.
     Weights can be negative */
  void set_tweights(node_id i, captype cap_source, captype cap_sink);

  /* Adds new edges 'SOURCE->i' and 'i->SINK' with corresponding weights
     Can be called multiple times for each node.
     Weights can be negative */
  void add_tweights(node_id i, captype cap_source, captype cap_sink);

  /* After the maxflow is computed, this function returns to which
     segment the node 'i' belongs (Graph::SOURCE or Graph::SINK) */
  termtype what_segment(node_id i);

  /* Computes the maxflow. Can be called only once. */
  flowtype maxflow();

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/

private:
  /* internal variables and functions */

  struct arc_forward_st;
  struct arc_reverse_st;

  typedef Int_to_ptr< sizeof(void*) >::type INTEGER;
#define MF_IS_ODD(a) ((INTEGER)(a) & 1)
#define MF_MAKE_ODD(a)  ((arc_forward *) ((INTEGER)(a) | 1))
#define MF_MAKE_EVEN(a) ((arc_forward *) ((INTEGER)(a) & (~1)))
#define MF_MAKE_ODD_REV(a)  ((arc_reverse *) ((INTEGER)(a) | 1))
#define MF_MAKE_EVEN_REV(a) ((arc_reverse *) ((INTEGER)(a) & (~1)))
#define MF_POINTER_TO_INTEGER(ptr) ((INTEGER) ptr)



  /* node structure */
  typedef struct node_st {
    /*
            Usually i->first_out is the first outgoing
            arc, and (i+1)->first_out-1 is the last outgoing arc.
            However, it is not always possible, since
            arcs are allocated in blocks, so arcs corresponding
            to two consecutive nodes may be in different blocks.

            If outgoing arcs for i are last in the arc block,
            then a different mechanism is used. i->first_out
            is odd in this case; the first outgoing arc
            is (a+1), and the last outgoing arc is
            ((arc_forward *)(a->shift))-1, where
            a = (arc_forward *) (((char *)(i->first_out)) + 1);

            Similar mechanism is used for incoming arcs.
    */
    arc_forward_st        *first_out;        /* first outcoming arc */
    arc_reverse_st        *first_in;        /* first incoming arc */

    arc_forward_st        *parent;        /* describes node's parent
                                                                           if MF_IS_ODD(parent) then MF_MAKE_EVEN(parent) points to 'arc_reverse',
                                                                           otherwise parent points to 'arc_forward' */

    node_st                        *next;                /* pointer to the next active node
                                                                           (or to itself if it is the last node in the list) */

    int                                TS;                        /* timestamp showing when DIST was computed */
    int                                DIST;                /* distance to the terminal */
    short                        is_sink;        /* flag showing whether the node is in the source or in the sink tree */

    captype                        tr_cap;                /* if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
                                                                           otherwise         -tr_cap is residual capacity of the arc node->SINK */
  } node;

  /* arc structures */
#define MF_NEIGHBOR_NODE(i, shift) ((node *) ((char *)(i) + (shift)))
#define MF_NEIGHBOR_NODE_REV(i, shift) ((node *) ((char *)(i) - (shift)))
  typedef struct arc_forward_st {
    INTEGER                        shift;                /* node_to = MF_NEIGHBOR_NODE(node_from, shift) */
    captype                        r_cap;                /* residual capacity */
    captype                        r_rev_cap;        /* residual capacity of the reverse arc*/
  } arc_forward;

  typedef struct arc_reverse_st {
    arc_forward                *sister;        /* reverse arc */
  } arc_reverse;

  /* 'pointer to node' structure */
  typedef struct nodeptr_st {
    node_st                        *ptr;
    nodeptr_st                *next;
  } nodeptr;

  typedef struct node_block_st {
    node                                        *current;
    struct node_block_st        *next;
    node                                        nodes[MF_NODE_BLOCK_SIZE];
  } node_block;

  typedef struct arc_for_block_st {
    char                                        *start;                /* the actual start address of this block.
                                                                                           May be different from 'this' since 'this'
                                                                                           must be at an even address. */
    arc_forward                                *current;
    struct arc_for_block_st        *next;
    arc_forward
    arcs_for[MF_ARC_BLOCK_SIZE]; /* all arcs must be at even addresses */
    union {
      arc_forward                        dummy;
      node                                *LAST_NODE;        /* used in graph consruction */
    }                                                LAST_NODE;
  } arc_for_block;

  typedef struct arc_rev_block_st {
    char                                        *start;                /* the actual start address of this block.
                                                                                           May be different from 'this' since 'this'
                                                                                           must be at an even address. */
    arc_reverse                                *current;
    struct arc_rev_block_st        *next;
    arc_reverse
    arcs_rev[MF_ARC_BLOCK_SIZE]; /* all arcs must be at even addresses */
    union {
      arc_reverse                        dummy;
      node                                *LAST_NODE;        /* used in graph consruction */
    }                                                LAST_NODE;
  } arc_rev_block;

  node_block                        *node_block_first;
  arc_for_block                *arc_for_block_first;
  arc_rev_block                *arc_rev_block_first;
  DBlock<nodeptr>                *nodeptr_block;

  void        (*error_function)(const char
                          *);        /* this function is called if a error occurs,
                                                                                   with a corresponding error message
                                                                                   (or exit(1) is called if it's nullptr) */

  flowtype                        flow;                /* total flow */

  /***********************************************************************/

  node                                *queue_first[2], *queue_last[2];        /* list of active nodes */
  nodeptr                                *orphan_first, *orphan_last;                /* list of pointers to orphans */
  int                                        TIME;                                                                /* monotonically increasing global counter */

  /***********************************************************************/

  /* functions for processing active list */
  void set_active(node *i);
  node *next_active();

  void prepare_graph();
  void maxflow_init();
  void augment(node *s_start, node *t_start, captype *cap_middle,
               captype *rev_cap_middle);
  void process_source_orphan(node *i);
  void process_sink_orphan(node *i);
};
/* graph.cpp */


//#include <stdio.h>
//#include "graph.h"

inline Graph::Graph(void (*err_function)(const char *))
{
  error_function = err_function;
  node_block_first = nullptr;
  arc_for_block_first = nullptr;
  arc_rev_block_first = nullptr;
  orphan_first = nullptr;
  orphan_last = nullptr;
  flow = 0;
}

inline Graph::~Graph()
{
  while (node_block_first) {
    node_block *next = node_block_first -> next;
    delete node_block_first;
    node_block_first = next;
  }

  while (arc_for_block_first) {
    arc_for_block *next = arc_for_block_first -> next;
    delete[] arc_for_block_first -> start;
    arc_for_block_first = next;
  }

  while (arc_rev_block_first) {
    arc_rev_block *next = arc_rev_block_first -> next;
    delete[] arc_rev_block_first -> start;
    arc_rev_block_first = next;
  }
}

inline Graph::node_id Graph::add_node()
{
  node *i;

  if (!node_block_first
      || node_block_first->current+1 > &node_block_first->nodes[MF_NODE_BLOCK_SIZE-1]) {
    node_block *next = node_block_first;
    node_block_first = (node_block *) new node_block;
    if (!node_block_first) {
      if (error_function) (*error_function)("Not enough memory!");
      exit(1);
    }
    node_block_first -> current = & ( node_block_first -> nodes[0] );
    node_block_first -> next = next;
  }

  i = node_block_first -> current ++;
  i -> first_out = (arc_forward *) 0;
  i -> first_in = (arc_reverse *) 0;

  i -> tr_cap = 0;

  return (node_id) i;
}

inline void Graph::add_edge(node_id from, node_id to, captype cap,
                            captype rev_cap)
{
  arc_forward *a_for;
  arc_reverse *a_rev;

  if (!arc_for_block_first
      || arc_for_block_first->current+1 >
      &arc_for_block_first->arcs_for[MF_ARC_BLOCK_SIZE]) {
    arc_for_block *next = arc_for_block_first;
    char *ptr = new char[sizeof(arc_for_block)+1];
    if (!ptr) {
      if (error_function) (*error_function)("Not enough memory!");
      exit(1);
    }
    if (MF_IS_ODD(ptr)) arc_for_block_first = (arc_for_block *) (ptr + 1);
    else              arc_for_block_first = (arc_for_block *) ptr;
    arc_for_block_first -> start = ptr;
    arc_for_block_first -> current = & ( arc_for_block_first -> arcs_for[0] );
    arc_for_block_first -> next = next;
  }

  if (!arc_rev_block_first
      || arc_rev_block_first->current+1 >
      &arc_rev_block_first->arcs_rev[MF_ARC_BLOCK_SIZE]) {
    arc_rev_block *next = arc_rev_block_first;
    char *ptr = new char[sizeof(arc_rev_block)+1];
    if (!ptr) {
      if (error_function) (*error_function)("Not enough memory!");
      exit(1);
    }
    if (MF_IS_ODD(ptr)) arc_rev_block_first = (arc_rev_block *) (ptr + 1);
    else              arc_rev_block_first = (arc_rev_block *) ptr;
    arc_rev_block_first -> start = ptr;
    arc_rev_block_first -> current = & ( arc_rev_block_first -> arcs_rev[0] );
    arc_rev_block_first -> next = next;
  }

  a_for = arc_for_block_first -> current ++;
  a_rev = arc_rev_block_first -> current ++;

  a_rev -> sister = (arc_forward *) from;
  a_for -> shift  = MF_POINTER_TO_INTEGER(to);
  a_for -> r_cap = cap;
  a_for -> r_rev_cap = rev_cap;

  ((node *)from) -> first_out =
    (arc_forward *) (MF_POINTER_TO_INTEGER(((node *)from) -> first_out) + 1);
  ((node *)to) -> first_in =
    (arc_reverse *) (MF_POINTER_TO_INTEGER(((node *)to) -> first_in) + 1);
}

inline void Graph::set_tweights(node_id i, captype cap_source, captype cap_sink)
{
  flow += (cap_source < cap_sink) ? cap_source : cap_sink;
  ((node*)i) -> tr_cap = cap_source - cap_sink;
}

inline void Graph::add_tweights(node_id i, captype cap_source, captype cap_sink)
{
  captype delta = ((node*)i) -> tr_cap;
  if (delta > 0) cap_source += delta;
  else           cap_sink   -= delta;
  flow += (cap_source < cap_sink) ? cap_source : cap_sink;
  ((node*)i) -> tr_cap = cap_source - cap_sink;
}

/*
        Converts arcs added by 'add_edge()' calls
        to a forward star graph representation.

        Linear time algorithm.
        No or little additional memory is allocated
        during this process
        (it may be necessary to allocate additional
        arc blocks, since arcs corresponding to the
        same node must be contiguous, i.e. be in one
        arc block.)
*/
inline void Graph::prepare_graph()
{
  node *i;
  arc_for_block *ab_for, *ab_for_first;
  arc_rev_block *ab_rev, *ab_rev_first, *ab_rev_scan;
  arc_forward *a_for;
  arc_reverse *a_rev, *a_rev_scan, *a_rev_tmp=new arc_reverse;
  node_block *nb;
  bool for_flag = false, rev_flag = false;
  INTEGER k;

  if (!arc_rev_block_first) {
    node_id from = add_node(), to = add_node();
    add_edge(from, to, 1, 0);
  }

  /* FIRST STAGE */
  a_rev_tmp->sister = nullptr;
  for (a_rev=arc_rev_block_first->current;
       a_rev<&arc_rev_block_first->arcs_rev[MF_ARC_BLOCK_SIZE]; a_rev++) {
    a_rev -> sister = nullptr;
  }

  ab_for = ab_for_first = arc_for_block_first;
  ab_rev = ab_rev_first = ab_rev_scan = arc_rev_block_first;
  a_for = &ab_for->arcs_for[0];
  a_rev = a_rev_scan = &ab_rev->arcs_rev[0];

  for (nb=node_block_first; nb; nb=nb->next) {
    for (i=&nb->nodes[0]; i<nb->current; i++) {
      /* outgoing arcs */
      k = MF_POINTER_TO_INTEGER(i -> first_out);
      if (a_for + k > &ab_for->arcs_for[MF_ARC_BLOCK_SIZE]) {
        if (k > MF_ARC_BLOCK_SIZE) {
          if (error_function) (*error_function)("# of arcs per node exceeds block size!");
          exit(1);
        }
        if (for_flag) ab_for = nullptr;
        else          {
          ab_for = ab_for -> next;
          ab_rev_scan = ab_rev_scan -> next;
        }
        if (ab_for == nullptr) {
          arc_for_block *next = arc_for_block_first;
          char *ptr = new char[sizeof(arc_for_block)+1];
          if (!ptr) {
            if (error_function) (*error_function)("Not enough memory!");
            exit(1);
          }
          if (MF_IS_ODD(ptr)) arc_for_block_first = (arc_for_block *) (ptr + 1);
          else              arc_for_block_first = (arc_for_block *) ptr;
          arc_for_block_first -> start = ptr;
          arc_for_block_first -> current = & ( arc_for_block_first -> arcs_for[0] );
          arc_for_block_first -> next = next;
          ab_for = arc_for_block_first;
          for_flag = true;
        } else a_rev_scan = &ab_rev_scan->arcs_rev[0];
        a_for = &ab_for->arcs_for[0];
      }
      if (ab_rev_scan) {
        a_rev_scan += k;
        i -> parent = (arc_forward *) a_rev_scan;
      } else i -> parent = (arc_forward *) a_rev_tmp;
      a_for += k;
      i -> first_out = a_for;
      ab_for -> LAST_NODE.LAST_NODE = i;

      /* incoming arcs */
      k = MF_POINTER_TO_INTEGER(i -> first_in);
      if (a_rev + k > &ab_rev->arcs_rev[MF_ARC_BLOCK_SIZE]) {
        if (k > MF_ARC_BLOCK_SIZE) {
          if (error_function) (*error_function)("# of arcs per node exceeds block size!");
          exit(1);
        }
        if (rev_flag) ab_rev = nullptr;
        else          ab_rev = ab_rev -> next;
        if (ab_rev == nullptr) {
          arc_rev_block *next = arc_rev_block_first;
          char *ptr = new char[sizeof(arc_rev_block)+1];
          if (!ptr) {
            if (error_function) (*error_function)("Not enough memory!");
            exit(1);
          }
          if (MF_IS_ODD(ptr)) arc_rev_block_first = (arc_rev_block *) (ptr + 1);
          else              arc_rev_block_first = (arc_rev_block *) ptr;
          arc_rev_block_first -> start = ptr;
          arc_rev_block_first -> current = & ( arc_rev_block_first -> arcs_rev[0] );
          arc_rev_block_first -> next = next;
          ab_rev = arc_rev_block_first;
          rev_flag = true;
        }
        a_rev = &ab_rev->arcs_rev[0];
      }
      a_rev += k;
      i -> first_in = a_rev;
      ab_rev -> LAST_NODE.LAST_NODE = i;
    }
    /* i is the last node in block */
    i -> first_out = a_for;
    i -> first_in  = a_rev;
  }

  /* SECOND STAGE */
  for (ab_for=arc_for_block_first; ab_for; ab_for=ab_for->next) {
    ab_for -> current = ab_for -> LAST_NODE.LAST_NODE -> first_out;
  }

  for ( ab_for=ab_for_first, ab_rev=ab_rev_first;
        ab_for;
        ab_for=ab_for->next, ab_rev=ab_rev->next )
    for ( a_for=&ab_for->arcs_for[0], a_rev=&ab_rev->arcs_rev[0];
          a_for<&ab_for->arcs_for[MF_ARC_BLOCK_SIZE];
          a_for++, a_rev++ ) {
      arc_forward *af;
      arc_reverse *ar;
      node *from;
      INTEGER shift = 0, shift_new;
      captype r_cap=0, r_rev_cap=0, r_cap_new, r_rev_cap_new;

      if (!(from=(node *)(a_rev->sister))) continue;
      af = a_for;
      ar = a_rev;

      do {
        ar -> sister = nullptr;

        shift_new = ((char *)(af->shift)) - (char *)from;
        r_cap_new = af -> r_cap;
        r_rev_cap_new = af -> r_rev_cap;
        if (shift) {
          af -> shift = shift;
          af -> r_cap = r_cap;
          af -> r_rev_cap = r_rev_cap;
        }
        shift = shift_new;
        r_cap = r_cap_new;
        r_rev_cap = r_rev_cap_new;

        af = -- from -> first_out;
        if ((arc_reverse *)(from->parent) != a_rev_tmp) {
          from -> parent = (arc_forward *)(((arc_reverse *)(from -> parent)) - 1);
          ar = (arc_reverse *)(from -> parent);
        }
      } while ( (from=(node *)(ar->sister)) );

      af -> shift = shift;
      af -> r_cap = r_cap;
      af -> r_rev_cap = r_rev_cap;
    }

  for (ab_for=arc_for_block_first; ab_for; ab_for=ab_for->next) {
    i = ab_for -> LAST_NODE.LAST_NODE;
    a_for = i -> first_out;
    ab_for -> current -> shift     = a_for -> shift;
    ab_for -> current -> r_cap     = a_for -> r_cap;
    ab_for -> current -> r_rev_cap = a_for -> r_rev_cap;
    a_for -> shift = MF_POINTER_TO_INTEGER(ab_for -> current + 1);
    i -> first_out = (arc_forward *) (((char *)a_for) - 1);
  }

  /* THIRD STAGE */
  for (ab_rev=arc_rev_block_first; ab_rev; ab_rev=ab_rev->next) {
    ab_rev -> current = ab_rev -> LAST_NODE.LAST_NODE -> first_in;
  }

  for (nb=node_block_first; nb; nb=nb->next)
    for (i=&nb->nodes[0]; i<nb->current; i++) {
      arc_forward *a_for_first, *a_for_last;

      a_for_first = i -> first_out;
      if (MF_IS_ODD(a_for_first)) {
        a_for_first = (arc_forward *) (((char *)a_for_first) + 1);
        a_for_last = (arc_forward *) ((a_for_first ++) -> shift);
      } else a_for_last = (i + 1) -> first_out;

      for (a_for=a_for_first; a_for<a_for_last; a_for++) {
        node *to = MF_NEIGHBOR_NODE(i, a_for -> shift);
        a_rev = -- to -> first_in;
        a_rev -> sister = a_for;
      }
    }

  for (ab_rev=arc_rev_block_first; ab_rev; ab_rev=ab_rev->next) {
    i = ab_rev -> LAST_NODE.LAST_NODE;
    a_rev = i -> first_in;
    ab_rev -> current -> sister = a_rev -> sister;
    a_rev -> sister = (arc_forward *) (ab_rev -> current + 1);
    i -> first_in = (arc_reverse *) (((char *)a_rev) - 1);
  }
  delete a_rev_tmp;
}

/* maxflow.cpp */

//#include <stdio.h>
//#include "graph.h"

/*
        special constants for node->parent
*/
#define MF_TERMINAL ( (arc_forward *) 1 )                /* to terminal */
#define MF_ORPHAN   ( (arc_forward *) 2 )                /* orphan */

#define MF_INFINITE_D 1000000000                /* infinite distance to the terminal */

/***********************************************************************/

/*
        Functions for processing active list.
        i->next points to the next node in the list
        (or to i, if i is the last node in the list).
        If i->next is nullptr iff i is not in the list.

        There are two queues. Active nodes are added
        to the end of the second queue and read from
        the front of the first queue. If the first queue
        is empty, it is replaced by the second queue
        (and the second queue becomes empty).
*/

inline void Graph::set_active(node *i)
{
  if (!i->next) {
    /* it's not in the list yet */
    if (queue_last[1]) queue_last[1] -> next = i;
    else               queue_first[1]        = i;
    queue_last[1] = i;
    i -> next = i;
  }
}

/*
        Returns the next active node.
        If it is connected to the sink, it stays in the list,
        otherwise it is removed from the list
*/
inline Graph::node * Graph::next_active()
{
  node *i;

  while ( 1 ) {
    if (!(i=queue_first[0])) {
      queue_first[0] = i = queue_first[1];
      queue_last[0]  = queue_last[1];
      queue_first[1] = nullptr;
      queue_last[1]  = nullptr;
      if (!i) return nullptr;
    }

    /* remove it from the active list */
    if (i->next == i) queue_first[0] = queue_last[0] = nullptr;
    else              queue_first[0] = i -> next;
    i -> next = nullptr;

    /* a node in the list is active iff it has a parent */
    if (i->parent) return i;
  }
}

/***********************************************************************/

inline void Graph::maxflow_init()
{
  node *i;
  node_block *nb;

  queue_first[0] = queue_last[0] = nullptr;
  queue_first[1] = queue_last[1] = nullptr;
  orphan_first = nullptr;

  for (nb=node_block_first; nb; nb=nb->next)
    for (i=&nb->nodes[0]; i<nb->current; i++) {
      i -> next = nullptr;
      i -> TS = 0;
      if (i->tr_cap > 0) {
        /* i is connected to the source */
        i -> is_sink = 0;
        i -> parent = MF_TERMINAL;
        set_active(i);
        i -> TS = 0;
        i -> DIST = 1;
      } else if (i->tr_cap < 0) {
        /* i is connected to the sink */
        i -> is_sink = 1;
        i -> parent = MF_TERMINAL;
        set_active(i);
        i -> TS = 0;
        i -> DIST = 1;
      } else {
        i -> parent = nullptr;
      }
    }
  TIME = 0;
}

/***********************************************************************/

inline void Graph::augment(node *s_start, node *t_start, captype *cap_middle,
                           captype *rev_cap_middle)
{
  node *i;
  arc_forward *a;
  captype bottleneck;
  nodeptr *np;


  /* 1. Finding bottleneck capacity */
  /* 1a - the source tree */
  bottleneck = *cap_middle;
  for (i=s_start; ; ) {
    a = i -> parent;
    if (a == MF_TERMINAL) break;
    if (MF_IS_ODD(a)) {
      a = MF_MAKE_EVEN(a);
      if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
      i = MF_NEIGHBOR_NODE_REV(i, a -> shift);
    } else {
      if (bottleneck > a->r_rev_cap) bottleneck = a -> r_rev_cap;
      i = MF_NEIGHBOR_NODE(i, a -> shift);
    }
  }
  if (bottleneck > i->tr_cap) bottleneck = i -> tr_cap;
  /* 1b - the sink tree */
  for (i=t_start; ; ) {
    a = i -> parent;
    if (a == MF_TERMINAL) break;
    if (MF_IS_ODD(a)) {
      a = MF_MAKE_EVEN(a);
      if (bottleneck > a->r_rev_cap) bottleneck = a -> r_rev_cap;
      i = MF_NEIGHBOR_NODE_REV(i, a -> shift);
    } else {
      if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
      i = MF_NEIGHBOR_NODE(i, a -> shift);
    }
  }
  if (bottleneck > - i->tr_cap) bottleneck = - i -> tr_cap;


  /* 2. Augmenting */
  /* 2a - the source tree */
  *rev_cap_middle += bottleneck;
  *cap_middle -= bottleneck;
  for (i=s_start; ; ) {
    a = i -> parent;
    if (a == MF_TERMINAL) break;
    if (MF_IS_ODD(a)) {
      a = MF_MAKE_EVEN(a);
      a -> r_rev_cap += bottleneck;
      a -> r_cap -= bottleneck;
      if (!a->r_cap) {
        /* add i to the adoption list */
        i -> parent = MF_ORPHAN;
        np = nodeptr_block -> New();
        np -> ptr = i;
        np -> next = orphan_first;
        orphan_first = np;
      }
      i = MF_NEIGHBOR_NODE_REV(i, a -> shift);
    } else {
      a -> r_cap += bottleneck;
      a -> r_rev_cap -= bottleneck;
      if (!a->r_rev_cap) {
        /* add i to the adoption list */
        i -> parent = MF_ORPHAN;
        np = nodeptr_block -> New();
        np -> ptr = i;
        np -> next = orphan_first;
        orphan_first = np;
      }
      i = MF_NEIGHBOR_NODE(i, a -> shift);
    }
  }
  i -> tr_cap -= bottleneck;
  if (!i->tr_cap) {
    /* add i to the adoption list */
    i -> parent = MF_ORPHAN;
    np = nodeptr_block -> New();
    np -> ptr = i;
    np -> next = orphan_first;
    orphan_first = np;
  }
  /* 2b - the sink tree */
  for (i=t_start; ; ) {
    a = i -> parent;
    if (a == MF_TERMINAL) break;
    if (MF_IS_ODD(a)) {
      a = MF_MAKE_EVEN(a);
      a -> r_cap += bottleneck;
      a -> r_rev_cap -= bottleneck;
      if (!a->r_rev_cap) {
        /* add i to the adoption list */
        i -> parent = MF_ORPHAN;
        np = nodeptr_block -> New();
        np -> ptr = i;
        np -> next = orphan_first;
        orphan_first = np;
      }
      i = MF_NEIGHBOR_NODE_REV(i, a -> shift);
    } else {
      a -> r_rev_cap += bottleneck;
      a -> r_cap -= bottleneck;
      if (!a->r_cap) {
        /* add i to the adoption list */
        i -> parent = MF_ORPHAN;
        np = nodeptr_block -> New();
        np -> ptr = i;
        np -> next = orphan_first;
        orphan_first = np;
      }
      i = MF_NEIGHBOR_NODE(i, a -> shift);
    }
  }
  i -> tr_cap += bottleneck;
  if (!i->tr_cap) {
    /* add i to the adoption list */
    i -> parent = MF_ORPHAN;
    np = nodeptr_block -> New();
    np -> ptr = i;
    np -> next = orphan_first;
    orphan_first = np;
  }


  flow += bottleneck;
}

/***********************************************************************/

inline void Graph::process_source_orphan(node *i)
{
  node *j;
  arc_forward *a0_for, *a0_for_first, *a0_for_last;
  arc_reverse *a0_rev, *a0_rev_first, *a0_rev_last;
  arc_forward *a0_min = nullptr, *a;
  nodeptr *np;
  int d, d_min = MF_INFINITE_D;

  /* trying to find a new parent */
  a0_for_first = i -> first_out;
  if (MF_IS_ODD(a0_for_first)) {
    a0_for_first = (arc_forward *) (((char *)a0_for_first) + 1);
    a0_for_last = (arc_forward *) ((a0_for_first ++) -> shift);
  } else a0_for_last = (i + 1) -> first_out;
  a0_rev_first = i -> first_in;
  if (MF_IS_ODD(a0_rev_first)) {
    a0_rev_first = (arc_reverse *) (((char *)a0_rev_first) + 1);
    a0_rev_last  = (arc_reverse *) ((a0_rev_first ++) -> sister);
  } else a0_rev_last = (i + 1) -> first_in;


  for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++)
    if (a0_for->r_rev_cap) {
      j = MF_NEIGHBOR_NODE(i, a0_for -> shift);
      if (!j->is_sink && (a=j->parent)) {
        /* checking the origin of j */
        d = 0;
        while ( 1 ) {
          if (j->TS == TIME) {
            d += j -> DIST;
            break;
          }
          a = j -> parent;
          d ++;
          if (a==MF_TERMINAL) {
            j -> TS = TIME;
            j -> DIST = 1;
            break;
          }
          if (a==MF_ORPHAN) {
            d = MF_INFINITE_D;
            break;
          }
          if (MF_IS_ODD(a))
            j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
          else
            j = MF_NEIGHBOR_NODE(j, a -> shift);
        }
        if (d<MF_INFINITE_D) { /* j originates from the source - done */
          if (d<d_min) {
            a0_min = a0_for;
            d_min = d;
          }
          /* set marks along the path */
          for (j=MF_NEIGHBOR_NODE(i, a0_for->shift); j->TS!=TIME; ) {
            j -> TS = TIME;
            j -> DIST = d --;
            a = j->parent;
            if (MF_IS_ODD(a))
              j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
            else
              j = MF_NEIGHBOR_NODE(j, a -> shift);
          }
        }
      }
    }
  for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++) {
    a0_for = a0_rev -> sister;
    if (a0_for->r_cap) {
      j = MF_NEIGHBOR_NODE_REV(i, a0_for -> shift);
      if (!j->is_sink && (a=j->parent)) {
        /* checking the origin of j */
        d = 0;
        while ( 1 ) {
          if (j->TS == TIME) {
            d += j -> DIST;
            break;
          }
          a = j -> parent;
          d ++;
          if (a==MF_TERMINAL) {
            j -> TS = TIME;
            j -> DIST = 1;
            break;
          }
          if (a==MF_ORPHAN) {
            d = MF_INFINITE_D;
            break;
          }
          if (MF_IS_ODD(a))
            j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
          else
            j = MF_NEIGHBOR_NODE(j, a -> shift);
        }
        if (d<MF_INFINITE_D) { /* j originates from the source - done */
          if (d<d_min) {
            a0_min = MF_MAKE_ODD(a0_for);
            d_min = d;
          }
          /* set marks along the path */
          for (j=MF_NEIGHBOR_NODE_REV(i,a0_for->shift); j->TS!=TIME; ) {
            j -> TS = TIME;
            j -> DIST = d --;
            a = j->parent;
            if (MF_IS_ODD(a))
              j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
            else
              j = MF_NEIGHBOR_NODE(j, a -> shift);
          }
        }
      }
    }
  }

  if ( (i->parent = a0_min) ) {
    i -> TS = TIME;
    i -> DIST = d_min + 1;
  } else {
    /* no parent is found */
    i -> TS = 0;

    /* process neighbors */
    for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++) {
      j = MF_NEIGHBOR_NODE(i, a0_for -> shift);
      if (!j->is_sink && (a=j->parent)) {
        if (a0_for->r_rev_cap) set_active(j);
        if (a!=MF_TERMINAL && a!=MF_ORPHAN && MF_IS_ODD(a)
            && MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a)->shift)==i) {
          /* add j to the adoption list */
          j -> parent = MF_ORPHAN;
          np = nodeptr_block -> New();
          np -> ptr = j;
          if (orphan_last) orphan_last -> next = np;
          else             orphan_first        = np;
          orphan_last = np;
          np -> next = nullptr;
        }
      }
    }
    for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++) {
      a0_for = a0_rev -> sister;
      j = MF_NEIGHBOR_NODE_REV(i, a0_for -> shift);
      if (!j->is_sink && (a=j->parent)) {
        if (a0_for->r_cap) set_active(j);
        if (a!=MF_TERMINAL && a!=MF_ORPHAN && !MF_IS_ODD(a) && MF_NEIGHBOR_NODE(j, a->shift)==i) {
          /* add j to the adoption list */
          j -> parent = MF_ORPHAN;
          np = nodeptr_block -> New();
          np -> ptr = j;
          if (orphan_last) orphan_last -> next = np;
          else             orphan_first        = np;
          orphan_last = np;
          np -> next = nullptr;
        }
      }
    }
  }
}

inline void Graph::process_sink_orphan(node *i)
{
  node *j;
  arc_forward *a0_for, *a0_for_first, *a0_for_last;
  arc_reverse *a0_rev, *a0_rev_first, *a0_rev_last;
  arc_forward *a0_min = nullptr, *a;
  nodeptr *np;
  int d, d_min = MF_INFINITE_D;

  /* trying to find a new parent */
  a0_for_first = i -> first_out;
  if (MF_IS_ODD(a0_for_first)) {
    a0_for_first = (arc_forward *) (((char *)a0_for_first) + 1);
    a0_for_last = (arc_forward *) ((a0_for_first ++) -> shift);
  } else a0_for_last = (i + 1) -> first_out;
  a0_rev_first = i -> first_in;
  if (MF_IS_ODD(a0_rev_first)) {
    a0_rev_first = (arc_reverse *) (((char *)a0_rev_first) + 1);
    a0_rev_last  = (arc_reverse *) ((a0_rev_first ++) -> sister);
  } else a0_rev_last = (i + 1) -> first_in;


  for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++)
    if (a0_for->r_cap) {
      j = MF_NEIGHBOR_NODE(i, a0_for -> shift);
      if (j->is_sink && (a=j->parent)) {
        /* checking the origin of j */
        d = 0;
        while ( 1 ) {
          if (j->TS == TIME) {
            d += j -> DIST;
            break;
          }
          a = j -> parent;
          d ++;
          if (a==MF_TERMINAL) {
            j -> TS = TIME;
            j -> DIST = 1;
            break;
          }
          if (a==MF_ORPHAN) {
            d = MF_INFINITE_D;
            break;
          }
          if (MF_IS_ODD(a))
            j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
          else
            j = MF_NEIGHBOR_NODE(j, a -> shift);
        }
        if (d<MF_INFINITE_D) { /* j originates from the sink - done */
          if (d<d_min) {
            a0_min = a0_for;
            d_min = d;
          }
          /* set marks along the path */
          for (j=MF_NEIGHBOR_NODE(i, a0_for->shift); j->TS!=TIME; ) {
            j -> TS = TIME;
            j -> DIST = d --;
            a = j->parent;
            if (MF_IS_ODD(a))
              j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
            else
              j = MF_NEIGHBOR_NODE(j, a -> shift);
          }
        }
      }
    }
  for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++) {
    a0_for = a0_rev -> sister;
    if (a0_for->r_rev_cap) {
      j = MF_NEIGHBOR_NODE_REV(i, a0_for -> shift);
      if (j->is_sink && (a=j->parent)) {
        /* checking the origin of j */
        d = 0;
        while ( 1 ) {
          if (j->TS == TIME) {
            d += j -> DIST;
            break;
          }
          a = j -> parent;
          d ++;
          if (a==MF_TERMINAL) {
            j -> TS = TIME;
            j -> DIST = 1;
            break;
          }
          if (a==MF_ORPHAN) {
            d = MF_INFINITE_D;
            break;
          }
          if (MF_IS_ODD(a))
            j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
          else
            j = MF_NEIGHBOR_NODE(j, a -> shift);
        }
        if (d<MF_INFINITE_D) { /* j originates from the sink - done */
          if (d<d_min) {
            a0_min = MF_MAKE_ODD(a0_for);
            d_min = d;
          }
          /* set marks along the path */
          for (j=MF_NEIGHBOR_NODE_REV(i,a0_for->shift); j->TS!=TIME; ) {
            j -> TS = TIME;
            j -> DIST = d --;
            a = j->parent;
            if (MF_IS_ODD(a))
              j = MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a) -> shift);
            else
              j = MF_NEIGHBOR_NODE(j, a -> shift);
          }
        }
      }
    }
  }

  if ( (i->parent = a0_min) ) {
    i -> TS = TIME;
    i -> DIST = d_min + 1;
  } else {
    /* no parent is found */
    i -> TS = 0;

    /* process neighbors */
    for (a0_for=a0_for_first; a0_for<a0_for_last; a0_for++) {
      j = MF_NEIGHBOR_NODE(i, a0_for -> shift);
      if (j->is_sink && (a=j->parent)) {
        if (a0_for->r_cap) set_active(j);
        if (a!=MF_TERMINAL && a!=MF_ORPHAN && MF_IS_ODD(a)
            && MF_NEIGHBOR_NODE_REV(j, MF_MAKE_EVEN(a)->shift)==i) {
          /* add j to the adoption list */
          j -> parent = MF_ORPHAN;
          np = nodeptr_block -> New();
          np -> ptr = j;
          if (orphan_last) orphan_last -> next = np;
          else             orphan_first        = np;
          orphan_last = np;
          np -> next = nullptr;
        }
      }
    }
    for (a0_rev=a0_rev_first; a0_rev<a0_rev_last; a0_rev++) {
      a0_for = a0_rev -> sister;
      j = MF_NEIGHBOR_NODE_REV(i, a0_for -> shift);
      if (j->is_sink && (a=j->parent)) {
        if (a0_for->r_rev_cap) set_active(j);
        if (a!=MF_TERMINAL && a!=MF_ORPHAN && !MF_IS_ODD(a) && MF_NEIGHBOR_NODE(j, a->shift)==i) {
          /* add j to the adoption list */
          j -> parent = MF_ORPHAN;
          np = nodeptr_block -> New();
          np -> ptr = j;
          if (orphan_last) orphan_last -> next = np;
          else             orphan_first        = np;
          orphan_last = np;
          np -> next = nullptr;
        }
      }
    }
  }
}

/***********************************************************************/

inline Graph::flowtype Graph::maxflow()
{
  node *i, *j, *current_node = nullptr, *s_start, *t_start=nullptr;
  captype *cap_middle=nullptr, *rev_cap_middle=nullptr;
  arc_forward *a_for, *a_for_first, *a_for_last;
  arc_reverse *a_rev, *a_rev_first, *a_rev_last;
  nodeptr *np, *np_next;

  prepare_graph();
  maxflow_init();
  nodeptr_block = new DBlock<nodeptr>(MF_NODEPTR_BLOCK_SIZE, error_function);

  while ( 1 ) {
    if ( (i=current_node) ) {
      i -> next = nullptr; /* remove active flag */
      if (!i->parent) i = nullptr;
    }
    if (!i) {
      if (!(i = next_active())) break;
    }

    /* growth */
    s_start = nullptr;

    a_for_first = i -> first_out;
    if (MF_IS_ODD(a_for_first)) {
      a_for_first = (arc_forward *) (((char *)a_for_first) + 1);
      a_for_last = (arc_forward *) ((a_for_first ++) -> shift);
    } else a_for_last = (i + 1) -> first_out;
    a_rev_first = i -> first_in;
    if (MF_IS_ODD(a_rev_first)) {
      a_rev_first = (arc_reverse *) (((char *)a_rev_first) + 1);
      a_rev_last = (arc_reverse *) ((a_rev_first ++) -> sister);
    } else a_rev_last = (i + 1) -> first_in;

    if (!i->is_sink) {
      /* grow source tree */
      for (a_for=a_for_first; a_for<a_for_last; a_for++)
        if (a_for->r_cap) {
          j = MF_NEIGHBOR_NODE(i, a_for -> shift);
          if (!j->parent) {
            j -> is_sink = 0;
            j -> parent = MF_MAKE_ODD(a_for);
            j -> TS = i -> TS;
            j -> DIST = i -> DIST + 1;
            set_active(j);
          } else if (j->is_sink) {
            s_start = i;
            t_start = j;
            cap_middle     = & ( a_for -> r_cap );
            rev_cap_middle = & ( a_for -> r_rev_cap );
            break;
          } else if (j->TS <= i->TS &&
                     j->DIST > i->DIST) {
            /* heuristic - trying to make the distance from j to the source shorter */
            j -> parent = MF_MAKE_ODD(a_for);
            j -> TS = i -> TS;
            j -> DIST = i -> DIST + 1;
          }
        }
      if (!s_start)
        for (a_rev=a_rev_first; a_rev<a_rev_last; a_rev++) {
          a_for = a_rev -> sister;
          if (a_for->r_rev_cap) {
            j = MF_NEIGHBOR_NODE_REV(i, a_for -> shift);
            if (!j->parent) {
              j -> is_sink = 0;
              j -> parent = a_for;
              j -> TS = i -> TS;
              j -> DIST = i -> DIST + 1;
              set_active(j);
            } else if (j->is_sink) {
              s_start = i;
              t_start = j;
              cap_middle     = & ( a_for -> r_rev_cap );
              rev_cap_middle = & ( a_for -> r_cap );
              break;
            } else if (j->TS <= i->TS &&
                       j->DIST > i->DIST) {
              /* heuristic - trying to make the distance from j to the source shorter */
              j -> parent = a_for;
              j -> TS = i -> TS;
              j -> DIST = i -> DIST + 1;
            }
          }
        }
    } else {
      /* grow sink tree */
      for (a_for=a_for_first; a_for<a_for_last; a_for++)
        if (a_for->r_rev_cap) {
          j = MF_NEIGHBOR_NODE(i, a_for -> shift);
          if (!j->parent) {
            j -> is_sink = 1;
            j -> parent = MF_MAKE_ODD(a_for);
            j -> TS = i -> TS;
            j -> DIST = i -> DIST + 1;
            set_active(j);
          } else if (!j->is_sink) {
            s_start = j;
            t_start = i;
            cap_middle     = & ( a_for -> r_rev_cap );
            rev_cap_middle = & ( a_for -> r_cap );
            break;
          } else if (j->TS <= i->TS &&
                     j->DIST > i->DIST) {
            /* heuristic - trying to make the distance from j to the sink shorter */
            j -> parent = MF_MAKE_ODD(a_for);
            j -> TS = i -> TS;
            j -> DIST = i -> DIST + 1;
          }
        }
      for (a_rev=a_rev_first; a_rev<a_rev_last; a_rev++) {
        a_for = a_rev -> sister;
        if (a_for->r_cap) {
          j = MF_NEIGHBOR_NODE_REV(i, a_for -> shift);
          if (!j->parent) {
            j -> is_sink = 1;
            j -> parent = a_for;
            j -> TS = i -> TS;
            j -> DIST = i -> DIST + 1;
            set_active(j);
          } else if (!j->is_sink) {
            s_start = j;
            t_start = i;
            cap_middle     = & ( a_for -> r_cap );
            rev_cap_middle = & ( a_for -> r_rev_cap );
            break;
          } else if (j->TS <= i->TS &&
                     j->DIST > i->DIST) {
            /* heuristic - trying to make the distance from j to the sink shorter */
            j -> parent = a_for;
            j -> TS = i -> TS;
            j -> DIST = i -> DIST + 1;
          }
        }
      }
    }

    TIME ++;

    if (s_start) {
      i -> next = i; /* set active flag */
      current_node = i;

      /* augmentation */
      augment(s_start, t_start, cap_middle, rev_cap_middle);
      /* augmentation end */

      /* adoption */
      while ( (np=orphan_first) ) {
        np_next = np -> next;
        np -> next = nullptr;

        while ( (np=orphan_first) ) {
          orphan_first = np -> next;
          i = np -> ptr;
          nodeptr_block -> Delete(np);
          if (!orphan_first) orphan_last = nullptr;
          if (i->is_sink) process_sink_orphan(i);
          else            process_source_orphan(i);
        }

        orphan_first = np_next;
      }
      /* adoption end */
    } else current_node = nullptr;
  }

  delete nodeptr_block;

  return flow;
}

/***********************************************************************/

inline Graph::termtype Graph::what_segment(node_id i)
{
  if (((node*)i)->parent && !((node*)i)->is_sink) return SOURCE;
  return SINK;
}

#undef MF_NODE_BLOCK_SIZE
#undef MF_ARC_BLOCK_SIZE
#undef MF_NODEPTR_BLOCK_SIZE
#undef MF_IS_ODD
#undef MF_MAKE_ODD
#undef MF_MAKE_EVEN
#undef MF_MAKE_ODD_REV
#undef MF_MAKE_EVEN_REV
#undef MF_POINTER_TO_INTEGER
#undef MF_NEIGHBOR_NODE
#undef MF_NEIGHBOR_NODE_REV
#undef MF_TERMINAL
#undef MF_ORPHAN
#undef MF_INFINITE_D


#endif
