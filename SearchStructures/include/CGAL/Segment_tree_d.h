// Copyright (c) 1997  ETH Zurich (SwitzerlanBd).
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
// Author(s)     : Gabriele Neyer

#ifndef CGAL_SEGMENT_TREE_D_H
#define CGAL_SEGMENT_TREE_D_H

#include <iostream>
#include <iterator>
#include <algorithm>
#include <list>
#include <vector>
#include <CGAL/Tree_base.h>

// A d-dimensional Segment Tree or a multilayer tree consisting of a Segment
// and other trees that are derived public Tree_base<C_Data, C_Window, 
// C_Interface> can be constructed within this class.
// C_Data: container class which contains the d-dim data the tree holds.
// C_Window: Query window -- a d-dimensional interval
// C_Interface: Interface for the class with functions that allow to access the
//             data. cf. file _interface.h for the requirements.

namespace CGAL {

template <class C_Data, class C_Window, class C_Interface>
class Segment_tree_d;

template <class C_Data, class C_Window, class C_Interface>
struct Segment_tree_node
  : public Tree_node_base<Segment_tree_node<C_Data, C_Window, C_Interface> >
{
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef typename C_Interface::Key Key;
  typedef  C_Interface Interface;
  typedef Tree_base< C_Data,  C_Window> Tree_base_type;
  std::list< C_Data> objects;
  Key left_key;
  Key right_key;
  Tree_base<C_Data, C_Window> *sublayer;
public:
  friend class Segment_tree_d< C_Data,  C_Window,  C_Interface>;

  typedef Tree_node_base<Segment_tree_node<C_Data, C_Window, C_Interface> >  Base;

  Segment_tree_node()
    : sublayer(0)
  {}

  Segment_tree_node(const Key& p_left_key,
		    const Key& p_right_key)
    : left_key(p_left_key), right_key(p_right_key), sublayer(0)
  {}

  Segment_tree_node(Segment_tree_node * p_left,
		    Segment_tree_node * p_right,
		    const Key& p_left_key,
		    const Key& p_right_key)
    : Base(p_left,p_right), left_key(p_left_key), right_key(p_right_key), sublayer(0)
  {}
  
  ~Segment_tree_node(){
    objects.clear();
    if (sublayer != 0)//(tree_base_type *)
      delete sublayer;
  }
};


template <class C_Data, class C_Window, class C_Interface>
class Segment_tree_d: public Tree_base< C_Data,  C_Window>
{
private:
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef  typename C_Interface::Key Key;
  typedef  C_Interface Interface;
public:
  typedef Tree_base<C_Data, C_Window> tbt;
protected:
  Tree_base<C_Data, C_Window> *sublayer_tree; 

  // type of a vertex
  // struct Segment_tree_node;
  
  friend struct Segment_tree_node<C_Data,C_Window,C_Interface>;
  typedef Segment_tree_node<C_Data,C_Window,C_Interface> Segment_tree_node_t;
  typedef Segment_tree_node<C_Data,C_Window,C_Interface> *link_type;
  
  static std::allocator<Segment_tree_node_t> alloc;
  
  C_Interface m_interface;
  bool is_built;


  bool is_less_equal(const Key& x, const Key&  y) const
  {
    return (!m_interface.comp(y,x));
  }   
  
  static link_type left(link_type x) { 
    return x->left_link;
  }
  static link_type right(link_type x) {
    return x->right_link; 
   }
  static link_type parent(link_type x) {
    return x->parent_link;
  }

  link_type header;
  link_type node;
  link_type rightmost(){return right(header);}
  link_type leftmost(){return left(header);}
  link_type root() const{
    if(header!=0)
      return header->parent_link; 
    return 0;
  }
  
  // returns true, if the object lies inside of win
  bool is_inside( C_Window const &win,  C_Data const& object) const
  {
    if(is_less_equal(m_interface.get_left_win(win), m_interface.get_left(object)) 
       && is_less_equal(m_interface.get_right(object),
			m_interface.get_right_win(win)))
    {
      return sublayer_tree->is_inside(win,object);
    }
    return false;
  }

  // this tree is not a recursion anchor
  bool is_anchor() const
  { return false;}  

  void insert_segment(link_type v,  C_Data& element)
  {
    if ((is_less_equal(m_interface.get_left(element), v->left_key) && 
	 is_less_equal(v->right_key, m_interface.get_right(element)))
	|| left(v)==CGAL_TREE_BASE_NULL)
      v->objects.push_back( element);
    else
     {
       if (!is_less_equal((*left(v)).right_key, m_interface.get_left(element)))
	 insert_segment(left(v), element);
       if (!is_less_equal(m_interface.get_right(element), (*right(v)).left_key))
	 insert_segment(right(v), element);
     }	
   }
  
  // according to the list of elements at vertex v, a sublayer tree for these
  // elements is created.
   void build_next_dimension(link_type v)
   {
     if(left(v)!=CGAL_TREE_BASE_NULL)
     {
       build_next_dimension(left(v));
       build_next_dimension(right(v));
     }
     if(! v->objects.empty())
     {
       typename std::list< C_Data>::iterator sub_first = v->objects.begin();
       typename std::list< C_Data>::iterator sub_last = v->objects.end();

       Tree_base<C_Data, C_Window> *g = sublayer_tree->clone();
       g->make_tree(sub_first, sub_last);
       v->sublayer = g;
       if (!v->sublayer->is_anchor())
       {
	 v->objects.clear();
       }
     }
   }

  link_type  new_Segment_tree_node_t(link_type l, link_type r, const Key& kl, const Key& kr)
  {
    Segment_tree_node_t node(l,r,kl,kr);
    Segment_tree_node_t* node_ptr = alloc.allocate(1);
    alloc.construct(node_ptr, node);
    return node_ptr;
  }

  link_type  new_Segment_tree_node_t(const Key& kl, const Key& kr)
  {
    Segment_tree_node_t node(kl,kr);
    Segment_tree_node_t* node_ptr = alloc.allocate(1);
    alloc.construct(node_ptr, node);
    return node_ptr;
  }

  link_type  new_Segment_tree_node_t()
  {
    Segment_tree_node_t node;
    Segment_tree_node_t* node_ptr = alloc.allocate(1);
    alloc.construct(node_ptr, node);
    return node_ptr;
  }

  // the skeleton of the segment tree is constructed here.
   void build_segment_tree(int n, link_type& leftchild, link_type& rightchild,
		   link_type& prevchild, link_type& leftmostlink,
		   int& index, int last, const std::vector<Key>& keys)
   { 
     // only two elements ==> two leaves and a parent is constructed
     if (n==2)
     {
       link_type vright;
       link_type vleft = new_Segment_tree_node_t
	 (keys[index], keys[index+1]);
       index++;
       if(index+1>last)
       {
         vright = new_Segment_tree_node_t
	   (keys[index], keys[index]);
       }
       else
       {
	 vright = new_Segment_tree_node_t
	   (keys[index], keys[index+1]);
       }
       index++;
       link_type vparent = new_Segment_tree_node_t
	 (vleft, vright, vleft->left_key, vright->right_key);

       vleft->parent_link = vparent;
       vright->parent_link = vparent;
       leftchild = vleft;
       rightchild = vright;
       prevchild = vparent;
       if(leftmostlink == CGAL_TREE_BASE_NULL)
	 leftmostlink = leftchild;
     }
     else
       // only one element ==> one leaf is constructed
       if(n==1)
       {
	 link_type vright;
	 if(index+1 > last){
	   vright = new_Segment_tree_node_t
	     (keys[index], keys[index]);
	 }
	 else{
	   vright = new_Segment_tree_node_t
	     (keys[index], keys[index+1]);
	 }
	 index++;

	 prevchild = vright;
	 rightchild = vright;
       }
       else
       {
	 // recursiv call for the construction. the interval is devided.
	 build_segment_tree(n - (int)n/2, leftchild, rightchild, 
			 prevchild, leftmostlink, index, last, keys);
	 link_type vparent = new_Segment_tree_node_t
	   (prevchild, CGAL_TREE_BASE_NULL, prevchild->left_key, prevchild->left_key);
	 prevchild->parent_link   = vparent;
	 build_segment_tree((int)n/2, leftchild, rightchild, 
			 prevchild, leftmostlink, index, last, keys);
	 vparent->right_link = prevchild;
	 prevchild->parent_link = vparent;
	 vparent->right_key = prevchild->right_key;
	 prevchild = vparent;
       }
   }

  void delete_tree(link_type v)
  {
    if(v->left_link!=CGAL_TREE_BASE_NULL)
    { 
      delete_tree(left(v));
      delete_tree(right(v));
    }
    delete_node(v);
  }	    

  void delete_node(Segment_tree_node_t* node_ptr)
  {
    alloc.destroy(node_ptr);
    alloc.deallocate(node_ptr,1);
  }

  // all elements that contain win are inserted into result
  template <class A>
  inline  
  A enclosing_query( C_Window const &win,
		     A result,
		     link_type v)
   {
     if(is_less_equal(m_interface.get_right_win(win), v->left_key) 
	|| is_less_equal(v->right_key,m_interface.get_left_win(win)))
       return result;
     if (v->sublayer!=0 && (!v->sublayer->is_anchor())) //(tree_base_type *)
     {
       Tree_base<C_Data, C_Window> *T = v->sublayer;

       std::list< C_Data> tmp_result;
       std::back_insert_iterator<std::list< C_Data> > tmp_back_inserter = 
	 std::back_inserter(tmp_result);
       T->enclosing_query(win, tmp_back_inserter);
       typename std::list<  C_Data>::iterator tmp = tmp_result.begin();
       while(tmp!=tmp_result.end())
       {
	 if(is_less_equal(m_interface.get_left(*tmp), 
			  m_interface.get_left_win(win)))
	 {
	   if(is_less_equal(m_interface.get_right_win(win), 
			    m_interface.get_right(*tmp)))
	     if(is_less_equal(v->left_key, m_interface.get_left_win(win)))
	       *result++=(*tmp);
	 }
	 tmp++;
       }
     }
     else
     {
       if(! v->objects.empty())
       {
	 typename std::list< C_Data>::iterator j=v->objects.begin();
	 while (j!= v->objects.end())
	 {
	   if(is_less_equal(m_interface.get_left(*j), 
			    m_interface.get_left_win(win)))
	   {
	     if(is_less_equal(m_interface.get_right_win(win), 
			      m_interface.get_right(*j)))
	       if(is_less_equal(v->left_key, m_interface.get_left_win(win)))
		 *result++=(*j);
	   }
	   j++;
	 }
       }
     }
     if(left(v))
     {
       enclosing_query(win, result, left(v));
       enclosing_query(win, result, right(v));
     }
     return result;
   }


  // all elements that habe non empty intersection with win are put into result
  template <class A>
  inline 
  A window_query( C_Window const &win,
		  A result,
		  const link_type& v) // af: was not const
   {
     if(is_less_equal(m_interface.get_right_win(win), v->left_key) || 
	is_less_equal(v->right_key,m_interface.get_left_win(win)))
       return result;
     if (v->sublayer!=0 && (!v->sublayer->is_anchor())) //(tree_base_type *)
     {
       Tree_base<C_Data, C_Window> *T = v->sublayer;

       std::list< C_Data> tmp_result;
       std::back_insert_iterator<std::list< C_Data> > tmp_back_inserter = 
	 std::back_inserter(tmp_result);
       T->window_query(win, tmp_back_inserter);
       typename std::list< C_Data>::iterator tmp = tmp_result.begin();
       while(tmp!=tmp_result.end())
       {
	 if(m_interface.comp(m_interface.get_left(*tmp), 
			   m_interface.get_left_win(win)))
	 {
	   if(is_less_equal(v->left_key, m_interface.get_left_win(win))){
	     *result++=(*tmp);
	   }
	 }
	 else
	 {
	   if(is_less_equal(v->left_key,m_interface.get_left(*tmp))){
	     *result++=(*tmp);
	   }
	 }
	 tmp++;
       }
     }
     else
     {
       if(! v->objects.empty())
       {
	 typename std::list< C_Data>::iterator j=v->objects.begin();
	 while (j!= v->objects.end())
	 {
	   if(m_interface.comp(m_interface.get_left(*j), m_interface.get_left_win(win)))
	   {
	     if(is_less_equal(v->left_key, m_interface.get_left_win(win)))
	     {
	       *result++=(*j);
	     }
	   }
	   else
	     if(is_less_equal(v->left_key,m_interface.get_left(*j)))
	     {
	       *result++=(*j);
	     }
	   j++;
	 }
       }
     }
     if(left(v))
     {
       window_query(win, result, left(v));
       window_query(win, result, right(v));
     }
     return result;
   }
  
  bool is_valid(const link_type& v) const // af:was not const reference
  {
    if (v->sublayer != 0)//(tree_base_type *)
    {
      Tree_base<C_Data, C_Window> *T=v->sublayer;
      if(! T->is_valid())
	return false;
    }
    if(left(v)!=CGAL_TREE_BASE_NULL)
    {
      if(!is_valid(left(v)))
	return false;
      if(!is_valid(right(v)))
	return false;
    }
    if(! v->objects.empty())
    {
//      true falls das Object das Segment enthaelt, 
//	  der parent aber das Segment nicht enthaelt.
      typename std::list< C_Data>::iterator j=v->objects.begin();
      link_type parent_of_v = parent(v);
      while (j!= v->objects.end())
      {
	if(!is_less_equal(m_interface.get_left(*j), v->left_key))
	  return false;
	if(!is_less_equal( v->right_key, m_interface.get_right(*j)))
	  return false;
	if (parent_of_v != root())
	  if((is_less_equal(m_interface.get_left(*j), parent_of_v->left_key))&& 
	     (is_less_equal( parent_of_v->right_key, 
			     m_interface.get_right(*j))))
	    return false;
	j++;
      }
    }
    return true;
  } 



public:

  // construction of a tree
  Segment_tree_d(Segment_tree_d const &sub_tree, bool):
    sublayer_tree(sub_tree.sublayer_tree->clone()), is_built(false), header(CGAL_TREE_BASE_NULL)
  {}

  // construction of a tree, definition of the prototype of sublayer tree
  Segment_tree_d(Tree_base<C_Data, C_Window> const &sub_tree):
    sublayer_tree(sub_tree.clone()), is_built(false), header(CGAL_TREE_BASE_NULL)
  {}

  // destruction 
  ~Segment_tree_d()
  {
    link_type v=root();
    if(v!=CGAL_TREE_BASE_NULL)
      delete_tree(v);
    if (header!=0)  
      delete_node(header);
    if(sublayer_tree!=0)
      delete sublayer_tree;
  }
   
  // clone creates a prototype
  Tree_base<C_Data, C_Window> *clone() const { 
    return new Segment_tree_d(*this, true); }

  bool make_tree(const typename std::list< C_Data>::iterator& beg,
		const typename std::list< C_Data>::iterator& end,
                 typename tbt::lit * =0){ 
    return make_tree_impl(beg,end);
  }

  #ifdef stlvector
  bool make_tree(const typename std::vector< C_Data>::iterator& beg, 
                 const typename std::vector< C_Data>::iterator& end,
                 typename tbt::vbit * =0){ 
    return make_tree_impl(beg,end);
  }
  #endif
  #ifdef carray
  bool make_tree(const C_Data *beg, 
                 const C_Data *end){
     return make_tree_impl(beg,end);
   }
  #endif

  // the tree is build according to Data [first,last)
  template<class A>
  inline 
  bool make_tree_impl(const A& first, const A& last)
  {
    if(!is_built)
      is_built = true;
    else
      return false;

    A count = first;
    int n=0;
    std::vector<Key> keys(2* std::distance(first, last));
    while(count!=last)
    {
      if (m_interface.comp(m_interface.get_left(*count),
			 m_interface.get_right(*count)))
      { 
	keys[n++]=m_interface.get_left(*count);
	keys[n++]=m_interface.get_right(*count);
      }
      else
      {
	CGAL_Tree_warning_msg(m_interface.comp(m_interface.get_left(*count),
					     m_interface.get_right(*count)), 
				  "invalid segment ignored");
      }
      count++;
    }

    if(n==0)
    {
      is_built = false;
      return true;
    }
    std::sort(keys.begin(), keys.end(), m_interface.comp);
    std::vector<Key> keys2(2*n + 1);
    int m=0;
    int num=1;
    keys2[0]=keys[0];
    for(m=1;m<n;m++)
    {
      if(m_interface.comp(keys[m],keys2[num-1])|| 
	 m_interface.comp(keys2[num-1],keys[m]))
      {
	keys2[num++]=keys[m];
      }
    }

    link_type leftchild;
    link_type rightchild;
    link_type prevchild;
    link_type leftmostlink = CGAL_TREE_BASE_NULL;

    int start = 0;
    build_segment_tree(num-1, leftchild, rightchild, prevchild, 
		      leftmostlink, start, num-1, keys2);

    header = new_Segment_tree_node_t();
    header->right_link = rightchild;
    header->parent_link = prevchild;
    prevchild->parent_link = prevchild;
    header->left_link = leftmostlink;

    A current = first;
    link_type r = root();
    do
    {
      if (m_interface.comp(m_interface.get_left(*current),
			 m_interface.get_right(*current)))
	insert_segment(r, *current);
    }while(++current!=last);

    link_type v=root();
    build_next_dimension(v);
    return true;
  }


  std::back_insert_iterator< std::list< C_Data> > window_query
          ( C_Window const &win, 
            std::back_insert_iterator< std::list< C_Data> > out,
            typename tbt::lbit * =0){
    return window_query_impl(win,out);
  }


  std::back_insert_iterator< std::vector< C_Data> > window_query
          ( C_Window const &win, 
            std::back_insert_iterator< std::vector< C_Data> > out,
            typename tbt::vbit * =0){
    return window_query_impl(win,out);
  }
  #ifdef carray
  C_Data *window_query( C_Window const &win, C_Data *out){
    return window_query_impl(win,out);
   }
  #endif

#ifdef ostreamiterator
  std::ostream_iterator< C_Data>  window_query( C_Window const &win, 
                     std::ostream_iterator< C_Data> out,
                     typename tbt::oit *dummy=0){
    return window_query_impl(win,out);
  }
#endif



  // all elements that ly inside win are inserted into result
  template <class A>
  inline A window_query_impl( C_Window const &win, 
			     A result,typename tbt::lbit * =0)
  {
    if(is_less_equal(m_interface.get_right_win(win), 
		     m_interface.get_left_win(win)))
    { 
      CGAL_Tree_warning_msg(m_interface.comp(m_interface.get_right_win(win), 
					       m_interface.get_left_win(win)),
				"invalid window -- query ignored");
      return result;
    }
    link_type v = root();
    if(v!=CGAL_TREE_BASE_NULL)
      return window_query(win, result, v);  
    return result;
  }

  
  std::back_insert_iterator< std::list< C_Data> > enclosing_query( 
	       C_Window const &win, 
               std::back_insert_iterator< std::list< C_Data> > out,
               typename tbt::lbit * =0){
    return enclosing_query_impl(win,out);
  }

  std::back_insert_iterator< std::vector< C_Data> > enclosing_query( 
	      C_Window const &win, 
              std::back_insert_iterator< std::vector< C_Data> > out,
              typename tbt::vbit * =0){
    return enclosing_query_impl(win,out);
  }


  #ifdef carray
  C_Data *enclosing_query( C_Window const &win, C_Data *out){
    return enclosing_query_impl(win,out);
   }
  #endif

#ifdef ostreamiterator
  std::ostream_iterator< C_Data>  
  enclosing_query( C_Window const &win, 
		   std::ostream_iterator< C_Data> out,
		   typename tbt::oit *dummy=0){
    return enclosing_query_impl(win,out);
  }
#endif



  // all objects that enclose win are inserted into result
  template <class A>
  inline
  A enclosing_query_impl( C_Window const &win, 
			  A result, typename tbt::lbit * =0)
  {
    if(is_less_equal(m_interface.get_right_win(win), 
		     m_interface.get_left_win(win)))
    { 
      CGAL_Tree_warning_msg(m_interface.comp(m_interface.get_right_win(win), 
					       m_interface.get_left_win(win)),
				"invalid window -- query ignored");
      return result;
    }
    link_type v = root();
    if(v!=CGAL_TREE_BASE_NULL)
      return enclosing_query(win, result, v);
    return result;
  }

  bool is_valid() const
  {
    link_type v= root();
    if(v!=CGAL_TREE_BASE_NULL)
      return is_valid(v);
    return true;
  }
};

template <class C_Data, class C_Window, class C_Interface>
std::allocator<Segment_tree_node<C_Data,C_Window,C_Interface> > 
    Segment_tree_d<C_Data,C_Window,C_Interface>::alloc;

} //namespace CGAL

#endif
