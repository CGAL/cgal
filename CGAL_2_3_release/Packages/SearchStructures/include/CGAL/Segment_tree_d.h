// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 09
//
// file          : include/CGAL/Segment_tree_d.h
// package       : SearchStructures (2.54)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// source        : include/CGAL/Segment_tree_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer
//
// coordinator   : Peter Widmayer, ETH Zurich
//
//
// email         : cgal@cs.uu.nl
//
// ======================================================================

#ifndef CGAL_Segment_tree_d__
#define CGAL_Segment_tree_d__

#include <iostream>
#include <iterator>
#include <algorithm>
#include <list>
#include <vector>
#include <CGAL/Tree_base.h>

// A d-dimensional Segment Tree or a multilayer tree consisting of a Segment
// and other trees that are derived public Tree_base<C_Data, C_Window, 
// C_Interface> can be constructed within this class.
// C_Data: container class which contains the d-dimensional data the tree holds.
// C_Window: Query window -- a d-dimensional interval
// C_Interface: Interface for the class with functions that allow to access the 
//             data. cf. file _interface.h for the requirements.

CGAL_BEGIN_NAMESPACE

template <class C_Data, class C_Window, class C_Interface>
class Segment_tree_d;

template <class C_Data, class C_Window, class C_Interface>
struct segment_tree_node: public tree_node_base
{
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef typename C_Interface::Key Key;
  typedef  C_Interface Interface;
  typedef tree_base< C_Data,  C_Window> tree_base_type;
  typedef Segment_tree_d< C_Data,  C_Window,  C_Interface> sT_d;
  std::list< C_Data> objects;
  Key left_key;
  Key right_key;
  tree_base<C_Data, C_Window> *sublayer;
public:
  friend sT_d;
  
  segment_tree_node(){
    sublayer = 0; //(tree_base_type *)
		      left_link = TREE_BASE_NULL;
    right_link = TREE_BASE_NULL;
  }
  segment_tree_node(segment_tree_node * p_left,
		    segment_tree_node * p_right,
		    const Key p_left_key,
		    const Key p_right_key)
    {
      left_link =p_left;
      right_link =p_right;
      left_key = p_left_key;
      right_key = p_right_key;
      sublayer = 0; //(tree_base_type *) 
			} 
  
  ~segment_tree_node(){
    objects.erase(objects.begin(), objects.end());
    if (sublayer != 0)//(tree_base_type *)
			delete sublayer;
  }
};


template <class C_Data, class C_Window, class C_Interface>
class Segment_tree_d: public tree_base< C_Data,  C_Window>
{
private:
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef  typename C_Interface::Key Key;
  typedef  C_Interface Interface;
public:
  typedef tree_base<C_Data, C_Window> tbt;
protected:
  typedef Segment_tree_d< C_Data,  C_Window,  C_Interface> sT_d;
  tree_base<C_Data, C_Window> *sublayer_tree; 
  
  // type of a vertex
  // struct segment_tree_node;
  
  friend segment_tree_node<C_Data,C_Window,C_Interface>;
  typedef segment_tree_node<C_Data,C_Window,C_Interface> segment_tree_node_t;
  typedef segment_tree_node<C_Data,C_Window,C_Interface> *link_type;
  
  C_Interface interface;
  bool is_build;


  bool is_less_equal(const Key x, const Key  y)
  {
    return (!interface.comp(y,x));
  }   
  
  static link_type& left(link_type x) { 
    return CGAL__static_cast(link_type&, (*x).left_link);
  }
  static link_type& right(link_type x) {
     return CGAL__static_cast(link_type&, (*x).right_link); 
   }
  static link_type& parent(link_type x) {
    return CGAL__static_cast(link_type&, (*x).parent_link);
  }

  link_type header;
  link_type node;
  link_type rightmost(){return right(header);}
  link_type leftmost(){return left(header);}
  link_type root(){
    if(header!=0)
      return CGAL__static_cast(link_type&, header->parent_link); 
    // return parent(header);
    else return 0;
  }
  
  // returns true, if the object lies inside of win
  bool is_inside( C_Window const &win,  C_Data const& object)
  {
    if(is_less_equal(interface.get_left_win(win), interface.get_left(object)) 
       && is_less_equal(interface.get_right(object),
			interface.get_right_win(win)))
    {
      return sublayer_tree->is_inside(win,object);
    }
    else
      return false;
  }

  // this tree is not a recursion anchor
  bool is_anchor()
  { return false;}  

  void insert_segment(link_type v,  C_Data& element)
  {
    if ((is_less_equal(interface.get_left(element), (*v).left_key) && 
	 is_less_equal((*v).right_key, interface.get_right(element)))
	|| left(v)==TREE_BASE_NULL)
      (*v).objects.insert((*v).objects.end(), element);
    else
     {
       if (!is_less_equal((*left(v)).right_key, interface.get_left(element)))
	 insert_segment(left(v), element);
       if (!is_less_equal(interface.get_right(element), (*right(v)).left_key))
	 insert_segment(right(v), element);
     }	
   }
  
  // according to the list of elements at vertex v, a sublayer tree for these
  // elements is created.
   void build_next_dimension(link_type v)
   {
     if(left(v)!=TREE_BASE_NULL)
     {
       build_next_dimension(left(v));
       build_next_dimension(right(v));
     }
     if(v->objects.size()>0)
     {
       typename std::list< C_Data>::iterator sub_first = v->objects.begin();
       typename std::list< C_Data>::iterator sub_last = v->objects.end();

       tree_base<C_Data, C_Window> *g = sublayer_tree->clone();
       g->make_tree(sub_first, sub_last);
       v->sublayer = g;
       if (!v->sublayer->is_anchor())
       {
	 sub_first = v->objects.begin();
	 sub_last = v->objects.end();
	 v->objects.erase(sub_first, sub_last);
       }
     }
   }

  // the sceleton of the segment tree is constructed here.
   void build_segment_tree(int n, link_type& leftchild, link_type& rightchild,
		   link_type& prevchild, link_type& leftmostlink,
		   int& index, int last, Key *keys)
   { 
     // only two elements ==> two leaves and a parent is constructed
     if (n==2)
     {
       link_type vright;
       link_type vleft = new segment_tree_node_t
	 (TREE_BASE_NULL, TREE_BASE_NULL, keys[index], keys[index+1]);
       index++;
       if(index+1>last)
       {
         vright = new segment_tree_node_t
	   (TREE_BASE_NULL, TREE_BASE_NULL, keys[index], keys[index]);
       }
       else
       {
	 vright = new segment_tree_node_t
	   (TREE_BASE_NULL, TREE_BASE_NULL, keys[index], keys[index+1]);
       }
       index++;
       link_type vparent = new segment_tree_node_t
	 (vleft, vright, vleft->left_key, vright->right_key);

       vleft->parent_link = vparent;
       vright->parent_link = vparent;
       leftchild = vleft;
       rightchild = vright;
       prevchild = vparent;
       if(leftmostlink == TREE_BASE_NULL)
	 leftmostlink = leftchild;
     }
     else
       // only one element ==> one leaf is constructed
       if(n==1)
       {
	 link_type vright;
	 if(index+1 > last){
	   vright = new segment_tree_node_t
	     (TREE_BASE_NULL, TREE_BASE_NULL, keys[index], keys[index]);
	 }
	 else{
	   vright = new segment_tree_node_t
	     (TREE_BASE_NULL, TREE_BASE_NULL, keys[index], keys[index+1]);
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
	 link_type vparent = new segment_tree_node_t
	   (prevchild, TREE_BASE_NULL, prevchild->left_key, prevchild->left_key);
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
    if(v->left_link!=TREE_BASE_NULL)
    { 
      delete_tree(left(v));
      delete_tree(right(v));
    }
    delete v;
  }	    


  // all elements that contain win are inserted into result
  template <class A>
  inline  
  A enclosing_query( C_Window const &win,
		     A result,
				         		      link_type v)
   {
     if(is_less_equal(interface.get_right_win(win), (*v).left_key) 
	|| is_less_equal((*v).right_key,interface.get_left_win(win)))
       return result;
     if (v->sublayer!=0 && (!v->sublayer->is_anchor())) //(tree_base_type *)
     {
       tree_base<C_Data, C_Window> *T = v->sublayer;

       std::list< C_Data> tmp_result;
       std::back_insert_iterator<std::list< C_Data> > tmp_back_inserter = 
	 std::back_inserter(tmp_result);
       (*T).enclosing_query(win, tmp_back_inserter);
       typename std::list<  C_Data>::iterator tmp = tmp_result.begin();
       while(tmp!=tmp_result.end())
       {
	 if(is_less_equal(interface.get_left(*tmp), 
			  interface.get_left_win(win)))
	 {
	   if(is_less_equal(interface.get_right_win(win), 
			    interface.get_right(*tmp)))
	     if(is_less_equal((*v).left_key, interface.get_left_win(win)))
	       *result++=(*tmp);
	 }
	 tmp++;
       }
     }
     else
     {
       if(v->objects.size()>0)
       {
	 typename std::list< C_Data>::iterator j=v->objects.begin();
	 while (j!= v->objects.end())
	 {
	   if(is_less_equal(interface.get_left(*j), 
			    interface.get_left_win(win)))
	   {
	     if(is_less_equal(interface.get_right_win(win), 
			      interface.get_right(*j)))
	       if(is_less_equal((*v).left_key, interface.get_left_win(win)))
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
						   link_type& v)
   {
     if(is_less_equal(interface.get_right_win(win), (*v).left_key) || 
	is_less_equal((*v).right_key,interface.get_left_win(win)))
       return result;
     if (v->sublayer!=0 && (!v->sublayer->is_anchor())) //(tree_base_type *)
     {
       tree_base<C_Data, C_Window> *T = v->sublayer;

       std::list< C_Data> tmp_result;
       std::back_insert_iterator<std::list< C_Data> > tmp_back_inserter = 
	 std::back_inserter(tmp_result);
       (*T).window_query(win, tmp_back_inserter);
       typename std::list< C_Data>::iterator tmp = tmp_result.begin();
       while(tmp!=tmp_result.end())
       {
	 if(interface.comp(interface.get_left(*tmp), 
			   interface.get_left_win(win)))
	 {
	   if(is_less_equal((*v).left_key, interface.get_left_win(win))){
	     *result++=(*tmp);
	   }
	 }
	 else
	 {
	   if(is_less_equal((*v).left_key,interface.get_left(*tmp))){
	     *result++=(*tmp);
	   }
	 }
	 tmp++;
       }
     }
     else
     {
       if(v->objects.size()>0)
       {
	 typename std::list< C_Data>::iterator j=v->objects.begin();
	 while (j!= v->objects.end())
	 {
	   if(interface.comp(interface.get_left(*j), interface.get_left_win(win)))
	   {
	     if(is_less_equal((*v).left_key, interface.get_left_win(win)))
	     {
	       *result++=(*j);
	     }
	   }
	   else
	     if(is_less_equal((*v).left_key,interface.get_left(*j)))
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
  
  bool is_valid(link_type& v)
  {
    if (v->sublayer != 0)//(tree_base_type *)
    {
      tree_base<C_Data, C_Window> *T=v->sublayer;
      if(! (*T).is_valid())
	return false;
    }
    if(left(v)!=TREE_BASE_NULL)
    {
      if(!is_valid(left(v)))
	return false;
      if(!is_valid(right(v)))
	return false;
    }
    if(v->objects.size()>0)
    {
//      true falls das Object das Segment enthaelt, 
//	  der parent aber das Segmetn nicht enthaelt.
      typename std::list< C_Data>::iterator j=v->objects.begin();
      link_type parent_of_v = parent(v);
      while (j!= v->objects.end())
      {
	if(!is_less_equal(interface.get_left(*j), (*v).left_key))
	  return false;
	if(!is_less_equal( (*v).right_key, interface.get_right(*j)))
	  return false;
	if (parent_of_v != root())
	  if((is_less_equal(interface.get_left(*j),(*parent_of_v).left_key))&& 
	     (is_less_equal( (*parent_of_v).right_key, 
			     interface.get_right(*j))))
	    return false;
	j++;
      }
    }
    return true;
  } 



public:

  // construction of a tree
  Segment_tree_d(Segment_tree_d const &sub_tree, bool):
    sublayer_tree(sub_tree.sublayer_tree->clone()), is_build(false)
  {
    header = TREE_BASE_NULL;
  }

  // construction of a tree, definition of the prototype of sublayer tree
  Segment_tree_d(tree_base<C_Data, C_Window> const &sub_tree):
    sublayer_tree(sub_tree.clone()), is_build(false)
  {
    header = TREE_BASE_NULL;
  }

  // destruction 
  ~Segment_tree_d()
  {
    link_type v=root();
    if(v!=TREE_BASE_NULL)
      delete_tree(v);
    if (header!=0)  
      delete header;
    if(sublayer_tree!=0)
      delete sublayer_tree;
  }
   
  // clone creates a prototype
  tree_base<C_Data, C_Window> *clone() const { 
    return new Segment_tree_d(*this, true); }

 bool make_tree(typename std::list< C_Data>::iterator& beg, 
                 typename std::list< C_Data>::iterator& end,
                 typename tbt::lit *dummy=0){ 
    return make_tree_impl(beg,end);
  }

  #ifdef stlvector
  bool make_tree(typename std::vector< C_Data>::iterator& beg, 
                 typename std::vector< C_Data>::iterator& end,
                 typename tbt::vbit *dummy=0){ 
    return make_tree_impl(beg,end);
  }
  #endif
  #ifdef carray
  bool make_tree(C_Data *beg, 
                 C_Data *end){
     return make_tree_impl(beg,end);
   }
  #endif

  // the tree is build according to Data [first,last)
  template<class A>
  inline 
  bool make_tree_impl(A& first,A& last)
  {
    if(!is_build)
      is_build = true;
    else
      return false;

    A count = first;
    int n=0;
    Key *keys = new Key[2*count_elements__C(first, last) + 1];
    while(count!=last)
    {
      if (interface.comp(interface.get_left(*count),
			 interface.get_right(*count)))
      { 
	keys[n++]=interface.get_left(*count);
	keys[n++]=interface.get_right(*count);
      }
      else
      {
	CGAL_Tree_warning_msg(interface.comp(interface.get_left(*count),
						 interface.get_right(*count)), 
				  "invalid segment ignored");
      }
      count++;
    }

    if(n==0)
    {
      is_build = false;
      return true;
    }
    std::sort(&keys[0], &keys[n], interface.comp);
    Key *keys2 = new Key[2*n + 1];
    int m=0;
    int num=1;
    keys2[0]=keys[0];
    for(m=1;m<n;m++)
    {
      if(interface.comp(keys[m],keys2[num-1])|| 
	 interface.comp(keys2[num-1],keys[m]))
      {
	keys2[num++]=keys[m];
      }
    }

    delete[] keys;
    link_type leftchild;
    link_type rightchild;
    link_type prevchild;
    link_type leftmostlink = TREE_BASE_NULL;

    int *start = new int(0);
    build_segment_tree(num-1, leftchild, rightchild, prevchild, 
		      leftmostlink, *start, num-1, keys2);
    delete[] keys2;
    delete start;

    header = new segment_tree_node_t();
    header->right_link = rightchild;
    header->parent_link = prevchild;
    prevchild->parent_link = prevchild;
    header->left_link = leftmostlink;

    A current = first;
    link_type r = root();
    do
    {
      if (interface.comp(interface.get_left(*current),
			 interface.get_right(*current)))
	insert_segment(r, *current);
    }while(++current!=last);

    link_type v=root();
    build_next_dimension(v);
    return true;
  }


  std::back_insert_iterator< std::list< C_Data> > window_query
          ( C_Window const &win, 
            std::back_insert_iterator< std::list< C_Data> > out,
            typename tbt::lbit *dummy=0){
    return window_query_impl(win,out);
  }


  std::back_insert_iterator< std::vector< C_Data> > window_query
          ( C_Window const &win, 
            std::back_insert_iterator< std::vector< C_Data> > out,
            typename tbt::vbit *dummy=0){
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
			     A result,typename tbt::lbit *dummy=0)
  {
    if(is_less_equal(interface.get_right_win(win), 
		     interface.get_left_win(win)))
    { 
      CGAL_Tree_warning_msg(interface.comp(interface.get_right_win(win), 
					       interface.get_left_win(win)),
				"invalid window -- query ignored");
      return result;
    }
    link_type v = root();
    if(v!=TREE_BASE_NULL)
      return window_query(win, result, v);  
    return result;
  }

  
  std::back_insert_iterator< std::list< C_Data> > enclosing_query( 
	       C_Window const &win, 
               std::back_insert_iterator< std::list< C_Data> > out,
               typename tbt::lbit *dummy=0){
    return enclosing_query_impl(win,out);
  }

  std::back_insert_iterator< std::vector< C_Data> > enclosing_query( 
	      C_Window const &win, 
              std::back_insert_iterator< std::vector< C_Data> > out,
              typename tbt::vbit *dummy=0){
    return enclosing_query_impl(win,out);
  }


  #ifdef carray
  C_Data *enclosing_query( C_Window const &win, C_Data *out){
    return enclosing_query_impl(win,out);
   }
  #endif

#ifdef ostreamiterator
  std::ostream_iterator< C_Data>  enclosing_query( C_Window const &win, 
                             std::ostream_iterator< C_Data> out,
                             typename tbt::oit *dummy=0){
    return enclosing_query_impl(win,out);
  }
#endif



  // all objects that enclose win are inserted into result
  template <class A>
  inline
  A enclosing_query_impl( C_Window const &win, 
		     A result,typename tbt::lbit *dummy=0)
  {
    if(is_less_equal(interface.get_right_win(win), 
		     interface.get_left_win(win)))
    { 
      CGAL_Tree_warning_msg(interface.comp(interface.get_right_win(win), 
					       interface.get_left_win(win)),
				"invalid window -- query ignored");
      return result;
    }
    link_type v = root();
    if(v!=TREE_BASE_NULL)
      return enclosing_query(win, result, v);
    return result;
  }

  bool is_valid()
  {
    link_type v= root();
    if(v!=TREE_BASE_NULL)
      return is_valid(v);
    return true;
  }
};
CGAL_END_NAMESPACE
#endif




